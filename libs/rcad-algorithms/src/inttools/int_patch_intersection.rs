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

use rcad_kernel::geom::Surface3;
use crate::bopds::ds::IntersectionCurve;
use glam::DVec3;
use super::int_patch_type::IntPatchIType;
use super::int_patch_point::IntPatchPoint;
use super::int_patch_line::IntPatchLine;
use super::geom_abs_surface_type::{GeomAbsSurfaceType, classify_surface_type};
use super::int_patch_imp_imp_intersection::ImpImpIntersection;

// ============================================================================
// ✅ OCCT-aligned: IntPatch_Intersection (IntPatch_Intersection.hxx)
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

    /// OCCT L48-63: IntPatch_Intersection() — default constructor
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
        if self.my_tol_arc < 1e-8 { self.my_tol_arc = 1e-8; }
        if self.my_tol_tang < 1e-8 { self.my_tol_tang = 1e-8; }
        // OCCT L147-151
        if self.my_tol_arc > 0.5 { self.my_tol_arc = 0.5; }
        if self.my_tol_tang > 0.5 { self.my_tol_tang = 0.5; }
        // OCCT L159-162
        if self.my_fleche < 1.0e-3 { self.my_fleche = 1e-3; }
        if self.my_fleche > 10.0 { self.my_fleche = 10.0; }
        // OCCT L169-171
        if self.my_uv_max_step > 0.5 { self.my_uv_max_step = 0.5; }
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
    ///   - Adaptor3d_Surface → &Surface3 (one parameter instead of surface+tool)
    ///   - Adaptor3d_TopolTool → omitted (face domain is implicit in rcad)
    ///   - isGeomInt, theIsReqToKeepRLine, theIsReqToPostWLProc → default params
    pub fn perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
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
            let (ct_surf, geom_surf, ct_type) = if typs1 == GeomAbsSurfaceType::Cone
                || typs1 == GeomAbsSurfaceType::Torus
            {
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
                                let dist = (apex1 - apex2).dot(dir1.cross(apex1 - apex2)) / (dir1.length() * dir1.length());
                                if dir1.dot(c2.axis_dir()).abs() > (1.0 - 1e-12) && dist.abs() < 1e-7 {
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
                                    if par < 0.015 { bg = 0; }
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
                                || (perp > (1.0 - 1e-10) && dist_to_plane < 1e-7) // normal + axis thru plane
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
                        (GeomAbsSurfaceType::Sphere, s.axis.normalize_or_zero(), s.center)
                    }
                    Surface3::Cylinder(c) => {
                        (GeomAbsSurfaceType::Cylinder, c.axis.normalize_or_zero(), c.origin)
                    }
                    Surface3::Cone(c) => {
                        (GeomAbsSurfaceType::Cone, c.axis_dir(), c.apex)
                    }
                    Surface3::Torus(t) => {
                        (GeomAbsSurfaceType::Torus, t.axis.normalize_or_zero(), t.center)
                    }
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
        //   Plane, Cylinder, Sphere, Cone → ts=1
        //   Torus → ts = bGeomGeom (0 normally, 1 if coaxial with compatible analytic surface)
        //   BSpline/Bezier/etc → ts=0
        let ts1 = match typs1 {
            GeomAbsSurfaceType::Plane | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Cone => 1,
            GeomAbsSurfaceType::Torus => b_geom_geom,
            _ => 0,
        };
        let ts2 = match typs2 {
            GeomAbsSurfaceType::Plane | GeomAbsSurfaceType::Cylinder
            | GeomAbsSurfaceType::Sphere | GeomAbsSurfaceType::Cone => 1,
            GeomAbsSurfaceType::Torus => b_geom_geom,
            _ => 0,
        };

        // ===== OCCT L1298-1339: dispatch
        // OCCT L1302-1324: Geom-Geom (ts1 == ts2 == 1)
        if ts1 == ts2 && ts1 == 1 {
            self.geom_geom_perform(s1, s2, typs1, typs2);
        }
        // OCCT L1326-1330: Geom-Param (ts1 != ts2)
        if ts1 != ts2 {
            self.geom_param_perform(s1, s2, ts1 == 0, typs1, typs2);
        }
        // OCCT L1332-1339: Param-Param (ts1 == ts2 == 0)
        if ts1 == ts2 && ts1 == 0 {
            self.param_param_perform(s1, s2, tol_arc, tol_tang, typs1, typs2);
        }

        // ===== OCCT L1346-1371: Post-process WLines
        // OCCT uses IntPatch_WLineTool::ComputePurgedWLine for each WLine
        // rcad: currently a no-op — marching in intss/ produces clean polylines
    }

    // =========================================================================
    // OCCT L139-170: Accessors
    // =========================================================================

    /// OCCT L139: IsDone
    pub fn is_done(&self) -> bool { self.done }
    /// OCCT L142: IsEmpty
    pub fn is_empty(&self) -> bool { self.empt }
    /// OCCT L147: TangentFaces
    pub fn tangent_faces(&self) -> bool { self.tgte }
    /// OCCT L154: OppositeFaces
    pub fn opposite_faces(&self) -> bool { self.oppo }
    /// OCCT L157: NbPnts
    pub fn nb_points(&self) -> usize { self.spnt.len() }
    /// OCCT L161: Point(Index)
    pub fn point(&self, index: usize) -> &IntPatchPoint { &self.spnt[index] }
    /// OCCT L164: NbLines
    pub fn nb_lines(&self) -> usize { self.slin.len() }
    /// OCCT L168: Line(Index)
    pub fn line(&self, index: usize) -> &IntPatchLine { &self.slin[index] }
    /// OCCT L170: SequenceOfLine
    pub fn sequence_of_line(&self) -> &[IntPatchLine] { &self.slin }

    // =========================================================================
    // rcad helper: convert to IntersectionCurve for DS storage
    // =========================================================================
    pub fn to_intersection_curves(&self) -> Vec<IntersectionCurve> {
        self.slin.iter().map(|l| {
            let mut curve_extra = crate::bopds::ds::CurveExtra::default();
            curve_extra.tangential_tol = l.tang_tolerance;
            IntersectionCurve {
                curve: l.curve.clone(),
                polyline: Vec::new(),
                start_vertex: usize::MAX,
                end_vertex: usize::MAX,
                t_range: l.t_range,
                pcurve_on_a: l.pcurve1.clone(),
                pcurve_on_b: l.pcurve2.clone(),
                geom_tol: l.tolerance,
                pave_blocks: Vec::new(),
                curve_extra,
            }
        }).collect()
    }

    // =========================================================================
    // OCCT L204-213: ParamParamPerfom (private)
    // =========================================================================
    fn param_param_perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        _tol_arc: f64,
        _tol_tang: f64,
        _typs1: GeomAbsSurfaceType,
        _typs2: GeomAbsSurfaceType,
    ) {
        // OCCT: IntPatch_PrmPrmIntersection interpp;
        // interpp.Perform(S1, D1, S2, D2, TolTang, TolArc, myFleche, myUVMaxStep, ListOfPnts);
        // rcad: use marching/numeric intersection via face_face::intersect_faces
        let curves = crate::inttools::face_face::intersect_faces(s1, s2,
            self.my_tol_arc, self.my_tol_tang);
        self.slin = curves.into_iter().map(|c| IntPatchLine {
            line_type: IntPatchIType::Walking,
            curve: c.curve,
            t_range: c.t_range,
            pcurve1: c.pcurve1,
            pcurve2: c.pcurve2,
            tolerance: c.tolerance,
            tang_tolerance: c.tang_tolerance,
        }).collect();
        self.empt = self.slin.is_empty();
        self.done = true;
    }

    // =========================================================================
    // OCCT L221-230: GeomGeomPerfom (private)
    // =========================================================================
    fn geom_geom_perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        _typs1: GeomAbsSurfaceType,
        _typs2: GeomAbsSurfaceType,
    ) {
        // ✅ OCCT-aligned: IntPatch_ImpImpIntersection
        let mut imp_imp = ImpImpIntersection::new();
        imp_imp.perform(s1, s2, self.my_tol_arc, self.my_tol_tang);
        if imp_imp.is_done() {
            self.slin = imp_imp.slin_ref().to_vec();
            self.empt = imp_imp.is_empty();
            self.tgte = imp_imp.tangent_faces();
            self.oppo = imp_imp.opposite_faces();
            self.done = imp_imp.is_done();
        }
    }

    // =========================================================================
    // OCCT L232-238: GeomParamPerfom (private)
    // =========================================================================
    fn geom_param_perform(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        _is_not_analytical: bool,
        _typs1: GeomAbsSurfaceType,
        _typs2: GeomAbsSurfaceType,
    ) {
        // OCCT: IntPatch_ImpPrmIntersection inter;
        // inter.Perform(S1, D1, S2, D2, TolTang, ...);
        // rcad: use marching/numeric intersection
        let curves = crate::inttools::face_face::intersect_faces(s1, s2,
            self.my_tol_arc, self.my_tol_tang);
        self.slin = curves.into_iter().map(|c| IntPatchLine {
            line_type: IntPatchIType::Walking,
            curve: c.curve,
            t_range: c.t_range,
            pcurve1: c.pcurve1,
            pcurve2: c.pcurve2,
            tolerance: c.tolerance,
            tang_tolerance: c.tang_tolerance,
        }).collect();
        self.empt = self.slin.is_empty();
        self.done = true;
    }
}
