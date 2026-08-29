// Created on: 1993-01-21
// Created by: Jacques GOUSSARD
// Copyright (c) 1993-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// --- Rust translation of OCCT IntPatch_Intersection ---
//
// This file is a line-by-line Rust translation of OCCT's
// IntPatch_Intersection.hxx / .cxx, adapted for rcad's data model.
// See the OCCT source for original algorithmic documentation.

use super::imp_imp_intersection::ImpImpIntersection;
use super::{
    classify_surface_type, GeomAbsSurfaceType, IntPatchLine, IntPatchPoint,
};
use glam::DVec3;
use rcad_kernel::geom::{Surface3, SurfaceEval};

// ============================================================================
// IntPatch_Intersection (IntPatch_Intersection.hxx)
//
// OCCT L240-254: member fields
// ============================================================================
pub struct IntPatchIntersection {
    // OCCT L240: done
    done: bool,
    // OCCT L241: empt
    empt: bool,
    // OCCT L242: tgte
    tgte: bool,
    // OCCT L243: oppo
    oppo: bool,
    // OCCT L244: NCollection_Sequence<IntPatch_Point> spnt
    spnt: Vec<IntPatchPoint>,
    // OCCT L245: NCollection_Sequence<handle<IntPatch_Line>> slin
    slin: Vec<IntPatchLine>,
    // OCCT L246-249
    my_tol_arc: f64,
    my_tol_tang: f64,
    my_uv_max_step: f64,
    my_fleche: f64,
    // OCCT L250-254
    my_is_start_pnt: bool,
    my_u1_start: f64,
    my_v1_start: f64,
    my_u2_start: f64,
    my_v2_start: f64,
}

impl IntPatchIntersection {
    // =========================================================================
    // OCCT L48-63: Constructors
    // =========================================================================

    /// OCCT L48-63: IntPatch_Intersection() 鈥?default constructor
    pub fn new() -> Self {
        Self {
            done: false,
            empt: true,
            tgte: false,
            oppo: false,
            spnt: Vec::new(),
            slin: Vec::new(),
            my_tol_arc: 0.0,
            my_tol_tang: 0.0,
            my_uv_max_step: 0.0,
            my_fleche: 0.0,
            my_is_start_pnt: false,
            my_u1_start: 0.0,
            my_v1_start: 0.0,
            my_u2_start: 0.0,
            my_v2_start: 0.0,
        }
    }

    // =========================================================================
    // OCCT L70-73: SetTolerances
    // =========================================================================

    /// OCCT L134-172: SetTolerances
    pub fn set_tolerances(&mut self, tol_arc: f64, tol_tang: f64, uv_max_step: f64, fleche: f64) {
        self.my_tol_arc = tol_arc;
        self.my_tol_tang = tol_tang;
        self.my_uv_max_step = uv_max_step;
        self.my_fleche = fleche;
        // OCCT L143-146
        if self.my_tol_arc < 1e-8 {
            self.my_tol_arc = 1e-8;
        }
        if self.my_tol_tang < 1e-8 {
            self.my_tol_tang = 1e-8;
        }
        // OCCT L147-151
        if self.my_tol_arc > 0.5 {
            self.my_tol_arc = 0.5;
        }
        if self.my_tol_tang > 0.5 {
            self.my_tol_tang = 0.5;
        }
        // OCCT L159-162
        if self.my_fleche < 1.0e-3 {
            self.my_fleche = 1e-3;
        }
        if self.my_fleche > 10.0 {
            self.my_fleche = 10.0;
        }
        // OCCT L169-171
        if self.my_uv_max_step > 0.5 {
            self.my_uv_max_step = 0.5;
        }
    }

    // =========================================================================
    // OCCT L86-94: Perform(S1, D1, S2, D2, TolArc, TolTang, ...)
    //
    // rcad: D1/D2 (TopolTool) replaced by Surface3 references.
    //       IntSurf_PntOn2S parameter (L115) not exposed in this overload.
    // =========================================================================

    /// OCCT L1066-1372: Perform
    ///
    /// rcad adaptation:
    ///   - Adaptor3d_Surface 鈫?&Surface3 (one parameter instead of surface+tool)
    ///   - Adaptor3d_TopolTool 鈫?the corrected FF UV rectangles uv1/uv2 plus
    ///     the faces' boundary pcurve arcs (Restriction mode) when available
    ///   - isGeomInt, theIsReqToKeepRLine, theIsReqToPostWLProc 鈫?default params
    pub fn perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        uv1: [f64; 4],
        uv2: [f64; 4],
        bnd1: &[super::so_on_bounds::BoundaryArc],
        bnd2: &[super::so_on_bounds::BoundaryArc],
        tol_arc: f64,
        tol_tang: f64,
    ) {
        // ===== OCCT L1076-1085: set tolerances
        self.my_tol_arc = tol_arc;
        self.my_tol_tang = tol_tang;
        if self.my_fleche <= 1e-12 {
            self.my_fleche = 0.01;
        }
        if self.my_uv_max_step <= 1e-12 {
            self.my_uv_max_step = 0.01;
        }

        // ===== OCCT L1087-1092: clear state
        self.done = false;
        self.spnt.clear();
        self.slin.clear();
        self.empt = true;
        self.tgte = false;
        self.oppo = false;

        // ===== OCCT L1094-1095: surface type classification
        let typs1 = classify_surface_type(s1);
        let typs2 = classify_surface_type(s2);

        // ===== OCCT L1098-1234: Cone/Torus special treatment =====
        let mut treat_as_bi_parametric = false;
        let mut b_geom_geom: i32 = 0;

        let has_cone_or_torus = typs1 == GeomAbsSurfaceType::Cone
            || typs2 == GeomAbsSurfaceType::Cone
            || typs1 == GeomAbsSurfaceType::Torus
            || typs2 == GeomAbsSurfaceType::Torus;

        if has_cone_or_torus {
            // OCCT L1104-1148: Cone special treatment
            let (ct_surf, geom_surf, ct_type) =
                if typs1 == GeomAbsSurfaceType::Cone || typs1 == GeomAbsSurfaceType::Torus {
                    (s1, s2, typs1)
                } else {
                    (s2, s1, typs2)
                };

            let mut b_to_check = false;
            let mut ct_axis_loc = DVec3::ZERO;
            let mut ct_axis_dir = DVec3::Z;

            if typs1 == GeomAbsSurfaceType::Cone || typs2 == GeomAbsSurfaceType::Cone {
                let cone_semi_angle = match ct_type {
                    GeomAbsSurfaceType::Cone => match ct_surf {
                        Surface3::Cone(c) => c.half_angle_rad,
                        _ => 0.0,
                    },
                    _ => match geom_surf {
                        Surface3::Cone(c) => c.half_angle_rad,
                        _ => 0.0,
                    },
                };
                let a1 = cone_semi_angle.abs();
                b_to_check = a1 < 0.02 || a1 > 1.55;

                if typs1 == typs2 {
                    let a2 = match geom_surf {
                        Surface3::Cone(c) => c.half_angle_rad.abs(),
                        _ => 0.0,
                    };
                    b_to_check = b_to_check || a2 < 0.02 || a2 > 1.55;

                    if a1 > 1.55 && a2 > 1.55 {
                        // Quasi-planes: if same domain, treat as canonic
                        match (ct_surf, geom_surf) {
                            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                                let apex1 = c1.apex_point();
                                let apex2 = c2.apex_point();
                                let dir1 = c1.axis_dir();
                                let dist = (apex1 - apex2).dot(dir1.cross(apex1 - apex2))
                                    / (dir1.length() * dir1.length());
                                if dir1.dot(c2.axis_dir()).abs() > (1.0 - 1e-12)
                                    && dist.abs() < 1e-7
                                {
                                    b_to_check = false;
                                }
                            }
                            _ => {}
                        }
                    }
                }

                treat_as_bi_parametric = b_to_check;
                if ct_type == GeomAbsSurfaceType::Cone {
                    if let Surface3::Cone(c) = ct_surf {
                        ct_axis_loc = c.apex;
                        ct_axis_dir = c.axis_dir();
                    }
                }
            }

            // OCCT L1150-1164: Torus special treatment
            if typs1 == GeomAbsSurfaceType::Torus || typs2 == GeomAbsSurfaceType::Torus {
                let tor_major = match ct_type {
                    GeomAbsSurfaceType::Torus => match ct_surf {
                        Surface3::Torus(t) => t.major_radius,
                        _ => 0.0,
                    },
                    _ => match geom_surf {
                        Surface3::Torus(t) => t.major_radius,
                        _ => 0.0,
                    },
                };
                let tor_minor = match ct_type {
                    GeomAbsSurfaceType::Torus => match ct_surf {
                        Surface3::Torus(t) => t.minor_radius,
                        _ => 0.0,
                    },
                    _ => match geom_surf {
                        Surface3::Torus(t) => t.minor_radius,
                        _ => 0.0,
                    },
                };
                b_to_check = tor_major > tor_minor;
                if typs1 == typs2 {
                    let tor2_minor = match geom_surf {
                        Surface3::Torus(t) => t.minor_radius,
                        _ => 0.0,
                    };
                    b_to_check = b_to_check && tor2_minor > 0.0;
                }
                if ct_type == GeomAbsSurfaceType::Torus {
                    if let Surface3::Torus(t) = ct_surf {
                        ct_axis_loc = t.center;
                        ct_axis_dir = t.axis.normalize_or_zero();
                    }
                }
            }

            // OCCT L1166-1233: Check axes for bGeomGeom
            if b_to_check {
                let ct_dir = ct_axis_dir;
                let ct_axis_line_origin = ct_axis_loc;

                let (gtype, gaxis_dir, gaxis_loc) = match geom_surf {
                    Surface3::Plane(p) => {
                        let gdir = p.normal;
                        // bGeomGeom = 1 if axes parallel, or if perpendicular & origin
                        let par = ct_dir.dot(gdir).abs();
                        let mut bg = 1i32;
                        if ct_type == GeomAbsSurfaceType::Cone {
                            bg = 1;
                            // Check if the cone semi-angle is very small
                            if let Surface3::Cone(c) = ct_surf {
                                if c.half_angle_rad.abs() < 0.02 {
                                    if par < 0.015 {
                                        bg = 0;
                                    }
                                }
                            }
                        } else {
                            // Torus-Plane
                            let normal = gdir;
                            let cross = ct_dir.cross(normal);
                            let perp = cross.length();
                            let surf_center = p.origin;
                            let dist_to_plane = (ct_axis_loc - surf_center).dot(normal).abs();
                            if perp < 1e-10  // parallel
                                || (perp > (1.0 - 1e-10) && dist_to_plane < 1e-7)
                            // normal + axis thru plane
                            {
                                bg = 1;
                            }
                        }
                        b_geom_geom = bg;
                        b_to_check = false;
                        (GeomAbsSurfaceType::Plane, gdir, p.origin)
                    }
                    Surface3::Sphere(s) => {
                        if (ct_axis_loc - s.center).cross(ct_dir).length().abs() < 1e-7 {
                            b_geom_geom = 1;
                        }
                        b_to_check = false;
                        (
                            GeomAbsSurfaceType::Sphere,
                            s.axis.normalize_or_zero(),
                            s.center,
                        )
                    }
                    Surface3::Cylinder(c) => (
                        GeomAbsSurfaceType::Cylinder,
                        c.axis.normalize_or_zero(),
                        c.origin,
                    ),
                    Surface3::Cone(c) => (GeomAbsSurfaceType::Cone, c.axis_dir(), c.apex),
                    Surface3::Torus(t) => (
                        GeomAbsSurfaceType::Torus,
                        t.axis.normalize_or_zero(),
                        t.center,
                    ),
                    _ => {
                        b_to_check = false;
                        (GeomAbsSurfaceType::OtherSurface, DVec3::Z, DVec3::ZERO)
                    }
                };

                // OCCT L1220-1227: remaining bToCheck cases (non-Plane/Sphere)
                if b_to_check {
                    let par = ct_dir.dot(gaxis_dir).abs();
                    let dist = (ct_axis_line_origin - gaxis_loc).cross(ct_dir).length();
                    if par > (1.0 - 1e-10) && dist < 1e-7 {
                        b_geom_geom = 1;
                    }
                }

                // OCCT L1229-1232: if bGeomGeom == 1, cancel TreatAsBiParametric
                if b_geom_geom == 1 {
                    treat_as_bi_parametric = false;
                }
            }
        }

        // OCCT L1242-1261: if TreatAsBiParametric, override surface types
        // Forces typs1/typs2 to BezierSurface to route through ImpPrm or PrmPrm.
        if treat_as_bi_parametric {
            // rcad: does not have ImpPrm/PrmPrm sub-algorithms in OCCT form yet.
            // When typsX is overridden to BezierSurface, OCCT routes to ImpPrm (one-side)
            // or PrmPrm (both sides). rcad: use marching for these paths.
        }

        // ===== OCCT L1264-1294: ts1/ts2 classification
        // Surface type definition: 1=analytic (geom), 0=parametric
        // OCCT L1264-1278:
        //   Plane, Cylinder, Sphere, Cone 鈫?ts=1
        //   Torus 鈫?ts = bGeomGeom (0 normally, 1 if coaxial with compatible analytic surface)
        //   BSpline/Bezier/etc 鈫?ts=0
        let ts1 = match typs1 {
            GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Cone => 1,
            GeomAbsSurfaceType::Torus => b_geom_geom,
            _ => 0,
        };
        let ts2 = match typs2 {
            GeomAbsSurfaceType::Plane
            | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Sphere
            | GeomAbsSurfaceType::Cone => 1,
            GeomAbsSurfaceType::Torus => b_geom_geom,
            _ => 0,
        };

        // ===== OCCT L1298-1339: dispatch
        // OCCT L1302-1324: Geom-Geom (ts1 == ts2 == 1)
        if ts1 == ts2 && ts1 == 1 {
            self.geom_geom_perform(s1, s2, uv1, uv2, bnd1, bnd2, typs1, typs2);
        }
        // OCCT L1326-1330: Geom-Param (ts1 != ts2)
        if ts1 != ts2 {
            self.geom_param_perform(s1, s2, uv1, uv2, ts1 == 0, typs1, typs2);
        }
        // OCCT L1332-1339: Param-Param (ts1 == ts2 == 0)
        if ts1 == ts2 && ts1 == 0 {
            self.param_param_perform(s1, s2, uv1, uv2, tol_arc, tol_tang, typs1, typs2);
        }

        // ===== OCCT L1346-1371: Post-process WLines
        // OCCT uses IntPatch_WLineTool::ComputePurgedWLine for each WLine
        // rcad: currently a no-op 鈥?marching in intss/ produces clean polylines
    }

    // =========================================================================
    // OCCT L139-170: Accessors
    // =========================================================================

    /// OCCT L139: IsDone
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT L142: IsEmpty
    pub fn is_empty(&self) -> bool {
        self.empt
    }
    /// OCCT L147: TangentFaces
    pub fn tangent_faces(&self) -> bool {
        self.tgte
    }
    /// OCCT L154: OppositeFaces
    pub fn opposite_faces(&self) -> bool {
        self.oppo
    }
    /// OCCT L157: NbPnts
    pub fn nb_points(&self) -> usize {
        self.spnt.len()
    }
    /// OCCT L161: Point(Index)
    pub fn point(&self, index: usize) -> &IntPatchPoint {
        &self.spnt[index]
    }
    /// OCCT L164: NbLines
    pub fn nb_lines(&self) -> usize {
        self.slin.len()
    }
    /// OCCT L168: Line(Index)
    pub fn line(&self, index: usize) -> &IntPatchLine {
        &self.slin[index]
    }
    /// OCCT L168: ChangeLine(Index) 鈥?mutable access for PutPointsOnLine.
    pub fn line_mut(&mut self, index: usize) -> &mut IntPatchLine {
        &mut self.slin[index]
    }
    /// OCCT L170: SequenceOfLine
    pub fn sequence_of_line(&self) -> &[IntPatchLine] {
        &self.slin
    }
    pub fn slin_mut(&mut self) -> &mut Vec<IntPatchLine> {
        &mut self.slin
    }

    // =========================================================================
    // OCCT L1659-1774: ParamParamPerfom (private)
    // =========================================================================
    fn param_param_perform(
        &mut self,
        _s1: &Surface3,
        _s2: &Surface3,
        _uv1: [f64; 4],
        _uv2: [f64; 4],
        _tol_arc: f64,
        _tol_tang: f64,
        _typs1: GeomAbsSurfaceType,
        _typs2: GeomAbsSurfaceType,
    ) {
        // OCCT L1669: IntPatch_PrmPrmIntersection interpp;
        // interpp.Perform(S1, D1, S2, D2, TolTang, TolArc, myFleche, myUVMaxStep, ListOfPnts);
        // rcad: IntPatch_PrmPrmIntersection (parametric-parametric walking) is
        // not yet ported to rcad-algo.  The analytic pairs (Geom-Geom) and the
        // analytic-parametric pairs (Geom-Param, ImpPrm) are the active paths
        // for the PaveFiller FF stage.
        self.done = false;
    }

    // =========================================================================
    // OCCT L1778-1905: GeomGeomPerfom (private)
    // =========================================================================
    fn geom_geom_perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        uv1: [f64; 4],
        uv2: [f64; 4],
        bnd1: &[super::so_on_bounds::BoundaryArc],
        bnd2: &[super::so_on_bounds::BoundaryArc],
        typs1: GeomAbsSurfaceType,
        typs2: GeomAbsSurfaceType,
    ) {
        // OCCT L1789-1797: IntPatch_ImpImpIntersection interii(...); if
        // (!interii.IsDone()) { done = false; ParamParamPerfom(...); return; }
        let mut imp_imp = ImpImpIntersection::new();
        imp_imp.perform(s1, s2, uv1, uv2, bnd1, bnd2, self.my_tol_arc, self.my_tol_tang);
        if !imp_imp.is_done() {
            self.done = false;
            // OCCT L1795: ParamParamPerfom(...) — rcad: parametric-parametric
            // marching is not ported (see param_param_perform).
            self.param_param_perform(s1, s2, uv1, uv2, self.my_tol_arc, self.my_tol_tang, typs1, typs2);
            return;
        }

        // OCCT L1799-1805: done = (GetStatus() == OK); empt = IsEmpty(); if (empt) return.
        self.done = imp_imp.status() == super::imp_imp_intersection::IntStatus::OK;
        self.empt = imp_imp.is_empty();
        if self.empt {
            return;
        }

        // OCCT L1807: const int aNbPointsInALine = 200;
        let a_nb_points_in_a_line = 200;

        // OCCT L1809-1813: tgte = interii.TangentFaces(); if (tgte) oppo = interii.OppositeFaces();
        self.tgte = imp_imp.tangent_faces();
        if self.tgte {
            self.oppo = imp_imp.opposite_faces();
        }

        // OCCT L1815-1838: IntPatch_ALineToWLine AToW(S1, S2, aNbPointsInALine);
        //   for each line: if (ArcType == IntPatch_Analytic) { isWLExist = true;
        //   AToW.MakeWLine(down_cast<ALine>(line), slin); } else {
        //   if (ArcType == IntPatch_Walking) WLine->EnablePurging(false);
        //   if ((ArcType != IntPatch_Restriction) || keepRLine) slin.Append(line); }
        let mut is_wl_exist = false;
        let a_to_w = super::a_line_to_w_line::ALineToWLine::new(s1, s2, a_nb_points_in_a_line);
        let keep_r_line = false;
        self.slin.clear();
        for i in 0..imp_imp.nb_lines() {
            if imp_imp.line(i).line_type == super::IntPatchIType::Analytic
                && imp_imp.line(i).a_curve.is_some()
            {
                if std::env::var("RCAD_FF_DEBUG").is_ok() {
                    eprintln!("[FF] geom_geom: ALine->WLine i={} n_vtx={}", i, imp_imp.line(i).nb_vertex());
                }
                is_wl_exist = true;
                a_to_w.make_wline(imp_imp.line_mut(i), &mut self.slin);
            } else {
                let line = imp_imp.line(i).clone();
                let mut line = line;
                if line.line_type == super::IntPatchIType::Walking {
                    line.is_purging_allowed = false;
                }
                if (line.line_type != super::IntPatchIType::Restriction) || keep_r_line {
                    self.slin.push(line);
                }
            }
        }

        // OCCT L1840-1843: copy spnt from interii.
        self.spnt.clear();
        for i in 0..imp_imp.nb_points() {
            self.spnt.push(imp_imp.point(i).clone());
        }

        // OCCT L1845-1848: if (typs1 == Cylinder && typs2 == Cylinder)
        //   IntPatch_WLineTool::JoinWLines(slin, spnt, S1, S2, TolTang).
        if typs1 == GeomAbsSurfaceType::Cylinder && typs2 == GeomAbsSurfaceType::Cylinder {
            super::w_line_tool::join_w_lines(&mut self.slin, s1, s2, self.my_tol_tang);
        }

        // OCCT L1850-1904: if (isWLExist) { build boxes, periodic array,
        //   critical points; IntPatch_WLineTool::ExtendTwoWLines(...); }
        if is_wl_exist {
            // OCCT L1852-1867: Bnd_Box2d aBx1/aBx2 from the surface UV
            // domains.  OCCT uses theS1->FirstUParameter()/FirstVParameter()
            // etc. — the GeomAdaptor_Surface was loaded (IntTools_FaceFace)
            // with the corrected face UV bounds, which is the rcad uv1/uv2.
            let mut a_bx1 = [uv1[0], uv1[2], uv1[1], uv1[3]];
            let mut a_bx2 = [uv2[0], uv2[2], uv2[1], uv2[3]];
            a_bx1 = enlarge_box2d(a_bx1, rcad_kernel::precision::PCONFUSION);
            a_bx2 = enlarge_box2d(a_bx2, rcad_kernel::precision::PCONFUSION);
            // OCCT L1869-1872: anArrOfPeriod.
            let an_arr_of_period = [
                period_of(s1, true),
                period_of(s1, false),
                period_of(s2, true),
                period_of(s2, false),
            ];
            // OCCT L1874-1894: aListOfCriticalPoints (cone apex / sphere poles).
            let mut a_list_of_critical_points: Vec<glam::DVec3> = Vec::new();
            match s1 {
                Surface3::Cone(c) => a_list_of_critical_points.push(c.apex_point()),
                Surface3::Sphere(sp) => {
                    a_list_of_critical_points.push(sp.point_at(0.0, std::f64::consts::FRAC_PI_2));
                    a_list_of_critical_points.push(sp.point_at(0.0, -std::f64::consts::FRAC_PI_2));
                }
                _ => {}
            }
            match s2 {
                Surface3::Cone(c) => a_list_of_critical_points.push(c.apex_point()),
                Surface3::Sphere(sp) => {
                    a_list_of_critical_points.push(sp.point_at(0.0, std::f64::consts::FRAC_PI_2));
                    a_list_of_critical_points.push(sp.point_at(0.0, -std::f64::consts::FRAC_PI_2));
                }
                _ => {}
            }
            super::w_line_tool::extend_two_w_lines(
                &mut self.slin,
                s1,
                s2,
                self.my_tol_tang,
                &an_arr_of_period,
                a_bx1,
                a_bx2,
                &a_list_of_critical_points,
            );
        }
    }

    // =========================================================================
    // OCCT L1909-2000: GeomParamPerfom (private)
    // =========================================================================
    fn geom_param_perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        uv1: [f64; 4],
        uv2: [f64; 4],
        is_not_analytical: bool,
        _typs1: GeomAbsSurfaceType,
        _typs2: GeomAbsSurfaceType,
    ) {
        // OCCT L1917: IntPatch_ImpPrmIntersection interip.
        let mut interip = super::imp_prm::ImpPrmIntersection::new();
        // OCCT L1918-1928: if (myIsStartPnt) SetStartPoint.
        if self.my_is_start_pnt {
            if is_not_analytical {
                interip.set_start_point(self.my_u1_start, self.my_v1_start);
            } else {
                interip.set_start_point(self.my_u2_start, self.my_v2_start);
            }
        }
        // OCCT L1930-1964: domains are always finite in the rcad FF path
        // (the corrected UV rectangles), so the else branch applies.
        interip.perform(
            s1,
            s2,
            uv1,
            uv2,
            self.my_tol_arc,
            self.my_tol_tang,
            self.my_fleche,
            self.my_uv_max_step,
        );

        // OCCT L1970-1999: copy the lines and points back.
        if interip.is_done() {
            self.done = true;
            self.empt = interip.is_empty();
            if !self.empt {
                let a_nb_lines = interip.nb_lines();
                for i in 1..=a_nb_lines {
                    let line = interip.line(i).clone();
                    if line.line_type != super::IntPatchIType::Walking {
                        self.slin.push(line);
                    }
                }
                for i in 1..=a_nb_lines {
                    let line = interip.line(i).clone();
                    if line.line_type == super::IntPatchIType::Walking {
                        self.slin.push(line);
                    }
                }
                for i in 1..=interip.nb_points() {
                    let p = interip.point(i);
                    self.spnt.push(super::IntPatchPoint {
                        p1: p.p3d,
                        p2: p.p3d,
                        u1: p.u1,
                        v1: p.v1,
                        u2: p.u2,
                        v2: p.v2,
                        tolerance: p.tolerance,
                    });
                }
            }
        }
    }
}

/// OCCT Bnd_Box2d::Enlarge(delta): grow the UV rectangle on all sides.
fn enlarge_box2d(mut b: [f64; 4], delta: f64) -> [f64; 4] {
    b[0] -= delta;
    b[1] -= delta;
    b[2] += delta;
    b[3] += delta;
    b
}

/// OCCT theS->UPeriod()/VPeriod() (0 when not periodic).
fn period_of(s: &Surface3, is_u: bool) -> f64 {
    if is_u {
        if s.is_u_periodic() { std::f64::consts::TAU } else { 0.0 }
    } else if s.is_v_periodic() {
        std::f64::consts::TAU
    } else {
        0.0
    }
}
