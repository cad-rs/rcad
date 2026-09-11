//! OCCT IntPatch_RstInt::PutVertexOnLine (IntPatch_RstInt.cxx L433-1355) with
//! its support classes, translated 1:1:
//!   - Adaptor2d_Line2d restriction arcs and the Adaptor3d_TopolTool
//!     rectangular-domain instantiation (Adaptor3d_TopolTool.cxx Initialize(S)
//!     L57-160, Initialize(C) L241-258, Has3d L991, Identical L1030;
//!     Adaptor3d_HVertex.cxx Parameter L43-46)
//!   - IntPatch_PolyLine (IntPatch_PolyLine.cxx L22-220 + IntPatch_Polygo.lxx)
//!   - IntPatch_PolyArc (IntPatch_PolyArc.cxx L21-227)
//!   - IntPatch_CSFunction (IntPatch_CSFunction.cxx L25-125)
//!   - the RstInt file statics Recadre (L49-131) / Tol3d (L133-146) /
//!     CoincideOnArc (L147-177) / VerifyTgline (L178-219) / GetLinePoint2d
//!     (L220-268) / FindParameter (L269-421)
//!
//! IntPatch_CurvIntSurf is the IntImp_IntCS.gxx instantiation, already ported
//! as the generic [`int_imp::int_cs::IntCS`]; IntPatch_SearchPnt is
//! Intf_InterferencePolygon2d, already ported in
//! [`crate::geomalgo::intf_interference_polygon2d`].
//!
//! Architecture notes:
//! - The OCCT Domain (Adaptor3d_TopolTool) appears here only in its
//!   rectangular-domain form: the domains handed to
//!   IntPatch_ImpPrmIntersection are the corrected UV rectangles of the two
//!   faces (IntTools_TopolTool), whose restriction arcs are the four
//!   Adaptor2d_Line2d rectangle edges, with Has3d() == false.
//! - The OCCT Adaptor3d_HVertex handles stored by
//!   IntPatch_Point::SetVertex are identified by `(arc index, endpoint index)`
//!   within the RectDomain; IntPatchVertex carries these ids in `vtx_on_s1` /
//!   `vtx_on_s2` (OCCT VertexOnS1/S2 handles).
//! - The polyline of the line is snapshotted once (OCCT's IntPatch_PolyLine
//!   keeps its own reference and is never refreshed while vertices are
//!   appended), so the walk mutates the IntPatchLine freely.

use std::marker::PhantomData;

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve2d, Line2d, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::super::int_imp::int_cs::{IntCS, ZerCSAccessors};
use super::super::int_imp::{CurveTool3d, PSurfaceTool};
use super::super::intf_interference::IntfPolygon2d;
use super::super::intf_interference_polygon2d::InterferencePolygon2d;
use super::transitions::{self, Transition};
use super::IntPatchIType;
use super::IntPatchLine;

/// OCCT IntPatch_PolyLine INITDEFLE (PolyLine.cxx L22).
const INITDEFLE: f64 = rcad_kernel::precision::PCONFUSION * 100.0;

/// OCCT IntPatch_RstInt.cxx L132: `const double Confusion`.
const CONFUSION: f64 = rcad_kernel::precision::CONFUSION;

/// OCCT Adaptor3d_TopolTool::myInfinite (the restriction bound clamp).
const TOPO_TOOL_INFINITE: f64 = 1e10;

/// OCCT GeomAbs_SurfaceType codes used by the RstInt statics.
const SURF_PLANE: i32 = 0;
const SURF_CYLINDER: i32 = 1;
const SURF_CONE: i32 = 2;
const SURF_SPHERE: i32 = 3;
const SURF_TORUS: i32 = 4;

fn surface_type(s: &Surface3) -> i32 {
    match s {
        Surface3::Plane(_) => SURF_PLANE,
        Surface3::Cylinder(_) => SURF_CYLINDER,
        Surface3::Cone(_) => SURF_CONE,
        Surface3::Sphere(_) => SURF_SPHERE,
        Surface3::Torus(_) => SURF_TORUS,
        _ => 10,
    }
}

/// OCCT GeomAdaptor_Surface::UResolution for the polyline tolerance — the
/// analytic branches map R3d through the radius; every other (incl.
/// Revolution/Extrusion) type takes Precision::Parametric(R3d) = R3d/100.
fn u_resolution(s: &Surface3, r3d: f64) -> f64 {
    match s {
        Surface3::Plane(_) => r3d,
        _ => r3d * 0.01,
    }
}

fn v_resolution(s: &Surface3, r3d: f64) -> f64 {
    match s {
        Surface3::Plane(_) => r3d,
        _ => r3d * 0.01,
    }
}

// ===========================================================================
// Adaptor2d_Line2d — the 2D line restriction arc of a rectangular domain
// ===========================================================================

/// OCCT Adaptor2d_Line2d for the rectangle-edge restrictions: the 2D line
/// `origin + t * dir` restricted to `first..last` (the constructor form used
/// by Adaptor3d_TopolTool::Initialize(S), dir being one of X / Y / NX / NY).
#[derive(Debug, Clone)]
pub struct Line2dArc {
    pub origin: DVec2,
    pub dir: DVec2,
    pub first: f64,
    pub last: f64,
}

impl Line2dArc {
    /// OCCT Adaptor2d_Line2d::Value(u).
    pub fn value(&self, t: f64) -> DVec2 {
        self.origin + self.dir * t
    }

    /// OCCT Adaptor2d_Line2d::D1 — the (unit) line direction.
    pub fn d1(&self) -> DVec2 {
        self.dir
    }

    /// OCCT Adaptor2d_Line2d::Line() — the gp_Lin2d.
    pub fn line(&self) -> Line2d {
        Line2d {
            origin: self.origin,
            direction: self.dir,
        }
    }

    /// OCCT IntPatch_HInterTool::Bounds(arc, PFirst, PLast).
    pub fn bounds(&self) -> (f64, f64) {
        (self.first, self.last)
    }
}

/// OCCT Adaptor3d_HVertex of a restriction-arc endpoint (Adaptor3d_TopolTool
/// Initialize(C) L241-258): the 2D point plus the IsSame identity inside the
/// domain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DomainVertex {
    /// The vertex's 2D point (C->Value(bound)).
    pub p2d: DVec2,
    /// Identity: (restriction index, endpoint index) — the handle.
    pub id: (u16, u16),
}

/// OCCT Adaptor3d_TopolTool over a rectangular surface domain
/// (Adaptor3d_TopolTool.cxx Initialize(S) L57-160): up to four
/// Adaptor2d_Line2d restriction arcs in the OCCT construction order
/// [Vinf(+X), Usup(+Y), Vsup(-X), Uinf(-Y)].
#[derive(Debug, Clone)]
pub struct RectDomain {
    pub u_inf: f64,
    pub v_inf: f64,
    pub u_sup: f64,
    pub v_sup: f64,
    /// The restriction arcs in OCCT construction order.
    pub restrictions: Vec<Line2dArc>,
}

impl RectDomain {
    /// OCCT Adaptor3d_TopolTool::Initialize(S) (L57-160) for a finite
    /// rectangular domain.
    pub fn from_bounds(u_inf: f64, u_sup: f64, v_inf: f64, v_sup: f64) -> Self {
        let mut restrictions: Vec<Line2dArc> = Vec::new();
        let u_infinfinite = rcad_kernel::precision::is_negative_infinite_value(u_inf);
        let u_supinfinite = rcad_kernel::precision::is_positive_infinite_value(u_sup);
        let v_infinfinite = rcad_kernel::precision::is_negative_infinite_value(v_inf);
        let v_supinfinite = rcad_kernel::precision::is_positive_infinite_value(v_sup);

        // L70-83: V = Vinf, direction +X, parameter in U.
        if !v_infinfinite {
            let deltap = (u_sup - u_inf).min(2.0 * TOPO_TOOL_INFINITE);
            let (pinf, psup) = if u_inf >= -TOPO_TOOL_INFINITE {
                (u_inf, u_inf + deltap)
            } else if u_sup <= TOPO_TOOL_INFINITE {
                (u_sup - deltap, u_sup)
            } else {
                (-TOPO_TOOL_INFINITE, TOPO_TOOL_INFINITE)
            };
            restrictions.push(Line2dArc {
                origin: DVec2::new(0.0, v_inf),
                dir: DVec2::new(1.0, 0.0),
                first: pinf,
                last: psup,
            });
        }
        // L85-99: U = Usup, direction +Y, parameter in V.
        if !u_supinfinite {
            let deltap = (v_sup - v_inf).min(2.0 * TOPO_TOOL_INFINITE);
            let (pinf, psup) = if v_inf >= -TOPO_TOOL_INFINITE {
                (v_inf, v_inf + deltap)
            } else if v_sup <= TOPO_TOOL_INFINITE {
                (v_sup - deltap, v_sup)
            } else {
                (-TOPO_TOOL_INFINITE, TOPO_TOOL_INFINITE)
            };
            restrictions.push(Line2dArc {
                origin: DVec2::new(u_sup, 0.0),
                dir: DVec2::new(0.0, 1.0),
                first: pinf,
                last: psup,
            });
        }
        // L101-116: V = Vsup, direction -X, parameter in -U.
        if !v_supinfinite {
            let deltap = (u_sup - u_inf).min(2.0 * TOPO_TOOL_INFINITE);
            let (pinf, psup) = if -u_sup >= -TOPO_TOOL_INFINITE {
                (-u_sup, -u_sup + deltap)
            } else if -u_inf <= TOPO_TOOL_INFINITE {
                (-u_inf - deltap, -u_inf)
            } else {
                (-TOPO_TOOL_INFINITE, TOPO_TOOL_INFINITE)
            };
            restrictions.push(Line2dArc {
                origin: DVec2::new(0.0, v_sup),
                dir: DVec2::new(-1.0, 0.0),
                first: pinf,
                last: psup,
            });
        }
        // L118-133: U = Uinf, direction -Y, parameter in -V.
        if !u_infinfinite {
            let deltap = (v_sup - v_inf).min(2.0 * TOPO_TOOL_INFINITE);
            let (pinf, psup) = if -v_sup >= -TOPO_TOOL_INFINITE {
                (-v_sup, -v_sup + deltap)
            } else if -v_inf <= TOPO_TOOL_INFINITE {
                (-v_inf - deltap, -v_inf)
            } else {
                (-TOPO_TOOL_INFINITE, TOPO_TOOL_INFINITE)
            };
            restrictions.push(Line2dArc {
                origin: DVec2::new(u_inf, 0.0),
                dir: DVec2::new(0.0, -1.0),
                first: pinf,
                last: psup,
            });
        }
        RectDomain {
            u_inf,
            v_inf,
            u_sup,
            v_sup,
            restrictions,
        }
    }

    /// OCCT Has3d() (L991-993).
    pub fn has_3d(&self) -> bool {
        false
    }

    /// OCCT Initialize(C) L241-258 + the vertex iterator: the arc's
    /// finite-bound endpoints, in [first, last] order.
    pub fn vertices_of(&self, arc_index: usize) -> Vec<DomainVertex> {
        let a = &self.restrictions[arc_index];
        let mut out = Vec::new();
        if a.first > -TOPO_TOOL_INFINITE {
            out.push(DomainVertex {
                p2d: a.value(a.first),
                id: (arc_index as u16, 0),
            });
        }
        if a.last < TOPO_TOOL_INFINITE {
            out.push(DomainVertex {
                p2d: a.value(a.last),
                id: (arc_index as u16, 1),
            });
        }
        out
    }

    /// OCCT Adaptor3d_HVertex::Parameter(C) (Adaptor3d_HVertex.cxx L43-46):
    /// ElCLib::Parameter(C->Line(), myPnt).
    pub fn vertex_parameter(&self, vtx: &DomainVertex, arc_index: usize) -> f64 {
        let a = &self.restrictions[arc_index];
        (vtx.p2d - a.origin).dot(a.dir)
    }

    /// OCCT Identical(V1, V2) (L1030-1034) — handle identity.
    pub fn identical(v1: &DomainVertex, v2: &DomainVertex) -> bool {
        v1.id == v2.id
    }
}

// ===========================================================================
// IntPatch_PolyLine — the line polygon (IntPatch_PolyLine.cxx L22-220)
// ===========================================================================

/// OCCT IntPatch_PolyLine: the polygon of a Walking or Restriction line, on
/// the surface selected by `on_first`.
pub struct PolyLine {
    typ: IntPatchIType,
    /// Snapshot of the line's points (OCCT keeps the handle; the polyline is
    /// never refreshed while vertices are appended).
    pnts: Vec<super::WLinePnt>,
    on_first: bool,
    my_error: f64,
    my_box: BndBox2d,
}

impl PolyLine {
    /// OCCT SetWLine (L32-39) / SetRLine (L42-49) — both call Prepare().
    pub fn new(line: &IntPatchLine, on_first: bool) -> Self {
        let typ = line.line_type;
        let mut pl = PolyLine {
            typ,
            pnts: line.wline_pnts.clone(),
            on_first,
            my_error: INITDEFLE,
            my_box: BndBox2d::new(),
        };
        pl.prepare();
        pl
    }

    /// OCCT NbPoints (L83-87).
    pub fn nb_points(&self) -> usize {
        self.pnts.len()
    }

    /// OCCT Prepare (L51-81).
    fn prepare(&mut self) {
        self.my_box.set_void();
        let n = self.nb_points();
        let eps_2 = self.my_error * self.my_error;

        let mut p1 = DVec2::ZERO;
        let mut p2 = DVec2::ZERO;
        if n >= 3 {
            p1 = self.point(1);
            p2 = self.point(2);
        }
        for i in 1..=n {
            let p3 = self.point(i);
            if i >= 3 {
                let v13 = p3 - p1;
                let v12 = p2 - p1;
                let d13_2 = v13.length_squared();
                let mut d_2;
                if d13_2 > eps_2 {
                    d_2 = v13.perp_dot(v12) * v13.perp_dot(v12) / d13_2;
                } else {
                    d_2 = eps_2;
                }
                if d_2 > self.my_error * self.my_error {
                    // try to compute deflection more precisely using parabola
                    // interpolation
                    let v23 = p3 - p2;
                    let d12 = v12.length();
                    let d23 = v23.length();
                    // compute parameter of P2 (assume parameters of P1,P3 are 0,1)
                    let mut tm = d12 / (d12 + d23);
                    if tm > 0.1 && tm < 0.9 {
                        tm -= (tm - 0.5) * 0.6;
                        let tm1mtm = tm * (1.0 - tm);
                        // coefficients of parabola
                        let ax = (tm * v13.x - v12.x) / tm1mtm;
                        let bx = (v12.x - tm * tm * v13.x) / tm1mtm;
                        let cx = p1.x;
                        let ay = (tm * v13.y - v12.y) / tm1mtm;
                        let by = (v12.y - tm * tm * v13.y) / tm1mtm;
                        let cy = p1.y;
                        // equations of lines P1-P2 and P2-P3
                        let a1 = v12.y / d12;
                        let b1 = -v12.x / d12;
                        let c1 = (p2.x * p1.y - p1.x * p2.y) / d12;
                        let a2 = v23.y / d23;
                        let b2 = -v23.x / d23;
                        let c2 = (p3.x * p2.y - p2.x * p3.y) / d23;
                        // points on parabola with max deflection
                        let t1 = -0.5 * (a1 * bx + b1 * by) / (a1 * ax + b1 * ay);
                        let t2 = -0.5 * (a2 * bx + b2 * by) / (a2 * ax + b2 * ay);
                        let xt1 = ax * t1 * t1 + bx * t1 + cx;
                        let yt1 = ay * t1 * t1 + by * t1 + cy;
                        let xt2 = ax * t2 * t2 + bx * t2 + cx;
                        let yt2 = ay * t2 * t2 + by * t2 + cy;
                        // max deflection on segments P1-P2 and P2-P3
                        let mut d1 = (a1 * xt1 + b1 * yt1 + c1).abs();
                        let d2v = (a2 * xt2 + b2 * yt2 + c2).abs();
                        if d2v > d1 {
                            d1 = d2v;
                        }
                        // select min deflection from linear and parabolic ones
                        if d1 * d1 < d_2 {
                            d_2 = d1 * d1;
                        }
                    }
                    if d_2 > self.my_error * self.my_error {
                        self.my_error = d_2.sqrt();
                    }
                }
                p1 = p2;
                p2 = p3;
            }
            self.my_box.update_xy(p3.x, p3.y);
        }
        self.my_box.enlarge(self.my_error);
    }

    /// OCCT ResetError (L89-92).
    pub fn reset_error(&mut self) {
        self.my_error = INITDEFLE;
    }

    /// OCCT Point (L93-137): the 2D parameters on the selected surface, with
    /// the 1e-7 nudge at the polyline ends (Walking lines only).
    pub fn point(&self, index: usize) -> DVec2 {
        let n = self.nb_points();
        let w = &self.pnts[index - 1];
        let (mut x, mut y) = if self.on_first {
            (w.u1, w.v1)
        } else {
            (w.u2, w.v2)
        };
        let mut dx = 0.0;
        let mut dy = 0.0;
        if self.typ == IntPatchIType::Walking {
            if index == 1 {
                let w1 = &self.pnts[1];
                let (x1, y1) = if self.on_first {
                    (w1.u1, w1.v1)
                } else {
                    (w1.u2, w1.v2)
                };
                dx = 0.0000001 * (x - x1);
                dy = 0.0000001 * (y - y1);
            } else if index == n {
                let w1 = &self.pnts[index - 2];
                let (x1, y1) = if self.on_first {
                    (w1.u1, w1.v1)
                } else {
                    (w1.u2, w1.v2)
                };
                dx = 0.0000001 * (x - x1);
                dy = 0.0000001 * (y - y1);
            }
        }
        x += dx;
        y += dy;
        DVec2::new(x, y)
    }
}

impl IntfPolygon2d for PolyLine {
    fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }
    fn bounding_mut(&mut self) -> &mut BndBox2d {
        &mut self.my_box
    }
    fn closed(&self) -> bool {
        false
    }
    fn deflection_over_estimation(&self) -> f64 {
        self.my_error
    }
    fn nb_segments(&self) -> i32 {
        (self.nb_points() as i32) - 1
    }
    fn segment(&self, the_index: i32) -> (DVec2, DVec2) {
        let i = the_index as usize;
        (self.point(i), self.point(i + 1))
    }
}

// ===========================================================================
// IntPatch_PolyArc — the restriction-arc polygon (IntPatch_PolyArc.cxx)
// ===========================================================================

fn min_max(a1: f64, a2: f64) -> (f64, f64) {
    if a1 < a2 {
        (a1, a2)
    } else {
        (a2, a1)
    }
}

/// OCCT IntPatch_PolyArc (cxx L21-227): the polygon of a restriction arc.
pub struct PolyArc {
    brise: Vec<DVec2>,
    param: Vec<f64>,
    offsetx: f64,
    offsety: f64,
    my_error: f64,
    my_box: BndBox2d,
    ferme: bool,
}

impl PolyArc {
    /// OCCT IntPatch_PolyArc ctor (L28-193).
    pub fn new(
        line: &Line2dArc,
        nb_sample: usize,
        a_pdeb: f64,
        a_pfin: f64,
        box_other_polygon: &BndBox2d,
    ) -> Self {
        let mut pdeb = a_pdeb;
        let mut pfin = a_pfin;

        let mut index_inf = nb_sample + 1;
        let mut index_sup = 0usize;

        // OCCT L62: BoxOtherPolygon.Get(bxmin, bymin, bxmax, bymax) — a void
        // box yields the infinite default rectangle.
        let (bxmin, bymin, bxmax, bymax) = match box_other_polygon.get() {
            Some((xmin, ymin, xmax, ymax)) => (xmin, ymin, xmax, ymax),
            None => (
                -TOPO_TOOL_INFINITE,
                -TOPO_TOOL_INFINITE,
                TOPO_TOOL_INFINITE,
                TOPO_TOOL_INFINITE,
            ),
        };
        let mut r = (bxmax - bxmin) + (bymax - bymin);
        let bx0 = (bxmax + bxmin) * 0.5;
        let by0 = (bymax + bymin) * 0.5;

        r *= 0.8;
        let mut r2 = r * r * 49.0;
        let mut nbloop = 0usize;

        let mut brise = vec![DVec2::ZERO; nb_sample];
        let mut param = vec![0.0f64; nb_sample];
        let mut my_box = BndBox2d::new();
        let mut my_error = 0.0f64;

        loop {
            nbloop += 1;
            let pas = (pfin - pdeb) / (nb_sample as f64 - 1.0);
            param[0] = pdeb;
            let p2d = line.value(pdeb);
            let (mut xs, mut ys) = (p2d.x, p2d.y);
            brise[0] = p2d;

            my_box.set_void();
            my_box.update_xy(brise[0].x, brise[0].y);
            my_error = 0.0;

            for i in 2..=nb_sample {
                param[i - 1] = pdeb + (i as f64 - 1.0) * pas;
                let p2d = line.value(param[i - 1]);
                let (x, y) = (p2d.x, p2d.y);
                brise[i - 1] = p2d;
                let xxs = 0.5 * (xs + x);
                let yys = 0.5 * (ys + y);
                //------------------------------------------------------------
                //-- On recherche le debut et la fin de la zone significative
                //------------------------------------------------------------
                // MSV: (see cda 002 H2) if segment is too large (>>r) we have
                //      a risk to jump through BoxOtherPolygon, therefore we should
                //      check this condition if the first one is failure.
                let is_mid_pt_in_box = (bx0 - xxs).abs() + (by0 - yys).abs() < r;
                let mut is_seg_out = true;
                if !is_mid_pt_in_box {
                    let d = (x - xs) * (x - xs) + (y - ys) * (y - ys);
                    if d > r2 {
                        let (xmin, xmax) = min_max(xs, x);
                        let (ymin, ymax) = min_max(ys, y);
                        is_seg_out = xmax < bxmin || xmin > bxmax || ymax < bymin || ymin > bymax;
                    }
                }
                if is_mid_pt_in_box || !is_seg_out {
                    // MSV: take the previous and the next segments too, because of
                    //      we check only the middle point (see BUC60946)
                    if index_inf > i {
                        index_inf = i.saturating_sub(2).max(1);
                    }
                    if index_sup < i {
                        index_sup = (i + 1).min(nb_sample);
                    }
                }

                my_box.update_xy(brise[i - 1].x, brise[i - 1].y);
                let pm = line.value(param[i - 1] - pas * 0.5);
                let xm0 = pm.x - xxs;
                let ym0 = pm.y - yys;
                let xm = (xm0 * xm0 + ym0 * ym0).sqrt();
                my_error = my_error.max(xm);
                xs = x;
                ys = y;
            }
            if index_inf > index_sup {
                r += r;
                r2 = r * r * 49.0;
                //-- std::cout<<" Le rayon : "<<r<<" est insuffisant "<<std::endl;
            } else {
                //----------------------------------------------
                //-- Si le nombre de points significatifs est
                //-- insuffisant, on reechantillonne une fois
                //-- encore
                //----------------------------------------------
                if index_sup.saturating_sub(index_inf) < nb_sample / 2 {
                    nbloop = 10;
                    pdeb = param[index_inf - 1];
                    pfin = param[index_sup - 1];
                    index_inf = nb_sample + 1;
                    index_sup = 0;
                }
            }
            if !((index_inf > index_sup) && nbloop <= 10) {
                break;
            }
        }
        my_error *= 1.2;
        if my_error < 0.00000001 {
            my_error = 0.00000001;
        }
        my_box.enlarge(my_error);

        let ferme = line.value(a_pdeb).distance(line.value(a_pfin)) <= 1e-7;

        PolyArc {
            brise,
            param,
            offsetx: 0.0,
            offsety: 0.0,
            my_error,
            my_box,
            ferme,
        }
    }

    /// OCCT NbPoints (L199-202).
    pub fn nb_points(&self) -> usize {
        self.brise.len()
    }

    /// OCCT Point (L204-212).
    pub fn point(&self, index: usize) -> DVec2 {
        if self.offsetx == 0.0 && self.offsety == 0.0 {
            return self.brise[index - 1];
        }
        let p = self.brise[index - 1];
        DVec2::new(p.x + self.offsetx, p.y + self.offsety)
    }

    /// OCCT Parameter (L214-217).
    pub fn parameter(&self, index: usize) -> f64 {
        self.param[index - 1]
    }

    /// OCCT SetOffset (L219-227).
    pub fn set_offset(&mut self, ox: f64, oy: f64) {
        let Some((xmin, ymin, xmax, ymax)) = self.my_box.get() else {
            self.offsetx = ox;
            self.offsety = oy;
            return;
        };
        let g = self.my_box.get_gap();

        self.my_box.set_void();

        self.my_box
            .update(xmin - self.offsetx, ymin - self.offsety, xmax - self.offsetx, ymax - self.offsety);
        self.offsetx = ox;
        self.offsety = oy;
        self.my_box
            .update(xmin + ox, ymin + oy, xmax + ox, ymax + oy);
        self.my_box.set_gap(g);
    }
}

impl IntfPolygon2d for PolyArc {
    fn bounding(&self) -> &BndBox2d {
        &self.my_box
    }
    fn bounding_mut(&mut self) -> &mut BndBox2d {
        &mut self.my_box
    }
    fn closed(&self) -> bool {
        self.ferme
    }
    fn deflection_over_estimation(&self) -> f64 {
        self.my_error
    }
    fn nb_segments(&self) -> i32 {
        (self.nb_points() as i32) - 1
    }
    fn segment(&self, the_index: i32) -> (DVec2, DVec2) {
        let i = the_index as usize;
        (self.point(i), self.point(i + 1))
    }
}

// ===========================================================================
// IntPatch_CSFunction — the curve-on-surface / surface function
// (IntPatch_CSFunction.cxx L25-125)
// ===========================================================================

/// The OCCT surface adaptor loaded with corrected bounds — the ThePSurface
/// payload of the IntImp generics.
pub struct BoundedSurface<'a> {
    pub surf: &'a Surface3,
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

/// ThePSurfaceTool marker for [`BoundedSurface`].
pub struct BoundedSurfaceTool<'a> {
    _life: PhantomData<&'a ()>,
}

impl<'a> PSurfaceTool for BoundedSurfaceTool<'a> {
    type Surface = BoundedSurface<'a>;
    fn value(s: &BoundedSurface<'a>, u: f64, v: f64) -> DVec3 {
        SurfaceEval::point_at(s.surf, u, v)
    }
    fn d1(s: &BoundedSurface<'a>, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        SurfaceEval::derivatives(s.surf, u, v)
    }
    fn first_u_parameter(s: &BoundedSurface<'a>) -> f64 {
        s.u0
    }
    fn last_u_parameter(s: &BoundedSurface<'a>) -> f64 {
        s.u1
    }
    fn first_v_parameter(s: &BoundedSurface<'a>) -> f64 {
        s.v0
    }
    fn last_v_parameter(s: &BoundedSurface<'a>) -> f64 {
        s.v1
    }
    fn u_resolution(_s: &BoundedSurface<'a>, r3d: f64) -> f64 {
        r3d * 0.01
    }
    fn v_resolution(_s: &BoundedSurface<'a>, r3d: f64) -> f64 {
        r3d * 0.01
    }
}

/// The TheCurveTool (IntPatch_HCurve2dTool) for the Line2d restriction arc.
pub struct Line2dArcTool;

impl CurveTool3d for Line2dArcTool {
    type Curve = Line2dArc;
    fn value(c: &Line2dArc, u: f64) -> DVec3 {
        let p = c.value(u);
        DVec3::new(p.x, p.y, 0.0)
    }
    fn d1(c: &Line2dArc, u: f64) -> (DVec3, DVec3) {
        let p = c.value(u);
        (
            DVec3::new(p.x, p.y, 0.0),
            DVec3::new(c.dir.x, c.dir.y, 0.0),
        )
    }
    fn first_parameter(c: &Line2dArc) -> f64 {
        c.first
    }
    fn last_parameter(c: &Line2dArc) -> f64 {
        c.last
    }
    fn resolution(_c: &Line2dArc, r3d: f64) -> f64 {
        r3d
    }
}

/// OCCT IntPatch_CSFunction (IntPatch_CSFunction.cxx L25-125): the function
/// F(X1, X2, X3) = S1(X1, X2) - S2(C(X3)), C being a 2D curve on S2.
pub struct CSFunction<'a> {
    surface1: &'a BoundedSurface<'a>,
    surface2: &'a Surface3,
    curve: &'a Line2dArc,
    p: DVec3,
    f: f64,
}

impl<'a> CSFunction<'a> {
    /// OCCT IntPatch_CSFunction(S1, C, S2) (L37-46): S1 is the surface on
    /// which the intersection is searched; C is a curve on the surface S2.
    pub fn new(s1: &'a BoundedSurface<'a>, c: &'a Line2dArc, s2: &'a Surface3) -> Self {
        CSFunction {
            surface1: s1,
            surface2: s2,
            curve: c,
            p: DVec3::ZERO,
            f: 0.0,
        }
    }
}

impl FunctionSetWithDerivatives for CSFunction<'_> {
    /// OCCT NbVariables (L49-52) / NbEquations (L55-58).
    fn nb_variables(&self) -> usize {
        3
    }
    fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value (L60-79).
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let psurf = SurfaceEval::point_at(self.surface1.surf, x[0], x[1]);
        let p2d = self.curve.value(x[2]);
        let pcurv = SurfaceEval::point_at(self.surface2, p2d.x, p2d.y);

        f[0] = psurf.x - pcurv.x;
        f[1] = psurf.y - pcurv.y;
        f[2] = psurf.z - pcurv.z;
        self.f = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
        self.p = (psurf + pcurv) * 0.5;
        true
    }

    /// OCCT Derivatives (L81-105).
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        let (_psurf, d1u, d1v) = SurfaceEval::derivatives(self.surface1.surf, x[0], x[1]);
        let p2d = self.curve.value(x[2]);
        let d2d = self.curve.d1();
        let (_pcurv, du, dv) = SurfaceEval::derivatives(self.surface2, p2d.x, p2d.y);
        let d1w = d2d.x * du + d2d.y * dv;

        df[0][0] = d1u.x;
        df[0][1] = d1v.x;
        df[0][2] = -d1w.x;
        df[1][0] = d1u.y;
        df[1][1] = d1v.y;
        df[1][2] = -d1w.y;
        df[2][0] = d1u.z;
        df[2][1] = d1v.z;
        df[2][2] = -d1w.z;
        true
    }

    /// OCCT Values (L107-130).
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        let (psurf, d1u, d1v) = SurfaceEval::derivatives(self.surface1.surf, x[0], x[1]);
        let p2d = self.curve.value(x[2]);
        let d2d = self.curve.d1();
        let (pcurv, du, dv) = SurfaceEval::derivatives(self.surface2, p2d.x, p2d.y);
        let d1w = d2d.x * du + d2d.y * dv;

        df[0][0] = d1u.x;
        df[0][1] = d1v.x;
        df[0][2] = -d1w.x;
        df[1][0] = d1u.y;
        df[1][1] = d1v.y;
        df[1][2] = -d1w.y;
        df[2][0] = d1u.z;
        df[2][1] = d1v.z;
        df[2][2] = -d1w.z;
        f[0] = psurf.x - pcurv.x;
        f[1] = psurf.y - pcurv.y;
        f[2] = psurf.z - pcurv.z;
        self.f = f[0] * f[0] + f[1] * f[1] + f[2] * f[2];
        self.p = (psurf + pcurv) * 0.5;
        true
    }
}

impl<'a> ZerCSAccessors<BoundedSurface<'a>, Line2dArc> for CSFunction<'a> {
    fn root(&self) -> f64 {
        self.f
    }
    fn point(&self) -> DVec3 {
        self.p
    }
    fn auxillar_surface(&self) -> &BoundedSurface<'a> {
        self.surface1
    }
    fn auxillar_curve(&self) -> &Line2dArc {
        self.curve
    }
}

// ===========================================================================
// RstInt file statics (IntPatch_RstInt.cxx L49-421)
// ===========================================================================

/// OCCT Recadre (L49-131): re-frame the refined UVs into the WLine point's
/// period neighborhood (2π quadric types; the torus V shift included).
fn recadre(
    type_s1: i32,
    type_s2: i32,
    wlin: &IntPatchLine,
    param: usize,
    u1: &mut f64,
    v1: &mut f64,
    u2: &mut f64,
    v2: &mut f64,
) {
    let two_pi = std::f64::consts::PI + std::f64::consts::PI;
    let nbpnts = wlin.nb_points();
    let param = param.clamp(1, nbpnts);
    let w = wlin.point(param - 1);
    let (u1p, v1p, u2p, v2p) = (w.u1, w.v1, w.u2, w.v2);
    match type_s1 {
        SURF_CYLINDER | SURF_CONE | SURF_SPHERE | SURF_TORUS => {
            while *u1 < (u1p - 1.5 * two_pi) {
                *u1 += two_pi;
            }
            while *u1 > (u1p + 1.5 * two_pi) {
                *u1 -= two_pi;
            }
        }
        _ => {}
    }
    if type_s1 == SURF_TORUS {
        while *v1 < (v1p - 1.5 * two_pi) {
            *v1 += two_pi;
        }
        while *v1 > (v1p + 1.5 * two_pi) {
            *v1 -= two_pi;
        }
    }
    match type_s2 {
        SURF_CYLINDER | SURF_CONE | SURF_SPHERE | SURF_TORUS => {
            while *u2 < (u2p - 1.5 * two_pi) {
                *u2 += two_pi;
            }
            while *u2 > (u2p + 1.5 * two_pi) {
                *u2 -= two_pi;
            }
        }
        _ => {}
    }
    if type_s2 == SURF_TORUS {
        while *v2 < (v1p - 1.5 * two_pi) {
            *v2 += two_pi;
        }
        while *v2 > (v2p + 1.5 * two_pi) {
            *v2 -= two_pi;
        }
    }
}

/// OCCT Tol3d (L133-146) — Has3d()==false for the rectangular geometric
/// domain: max(tolDef, Confusion).
fn tol3d(tol_def: f64) -> f64 {
    if tol_def < CONFUSION {
        CONFUSION
    } else {
        tol_def
    }
}

/// OCCT CoincideOnArc (L147-177).
fn coincide_on_arc(
    ptsommet: DVec3,
    a: &Line2dArc,
    arc_index: usize,
    surf: &Surface3,
    toler: f64,
    domain: &RectDomain,
) -> Option<DomainVertex> {
    let mut distmin = f64::MAX;
    let tolarc = toler.max(tol3d(0.0));

    let mut result = None;
    for vtx in domain.vertices_of(arc_index) {
        let p2d = vtx.p2d;
        let point = SurfaceEval::point_at(surf, p2d.x, p2d.y);
        let dist = point.distance(ptsommet);
        let tol = tolarc.max(tol3d(0.0));

        if dist <= tol && dist <= distmin {
            // the best coincidence
            distmin = dist;
            result = Some(vtx);
        }
    }
    if distmin < f64::MAX {
        result
    } else {
        None
    }
}

/// OCCT VerifyTgline (L178-219).
fn verify_tgline(wlin: &IntPatchLine, param: usize, tol: f64, tgl: &mut DVec3) {
    if tgl.x.abs() < tol && tgl.y.abs() < tol && tgl.z.abs() < tol {
        //-- On construit une tangente plus grande
        //-- (Eviter des points tres proches ds Walking)
        let nbpt = wlin.nb_points();
        let mut forward = (nbpt - param) >= (param - 1);
        for _ in 0..2 {
            if forward {
                for i in param + 1..=nbpt {
                    let t = wlin.point(i - 1).p3d - wlin.point(param - 1).p3d;
                    if t.x.abs() >= tol || t.y.abs() >= tol || t.z.abs() >= tol {
                        *tgl = t;
                        return;
                    }
                }
            } else {
                for i in (1..param).rev() {
                    let t = wlin.point(param - 1).p3d - wlin.point(i - 1).p3d;
                    if t.x.abs() >= tol || t.y.abs() >= tol || t.z.abs() >= tol {
                        *tgl = t;
                        return;
                    }
                }
            }
            forward = !forward;
        }
    }
}

/// OCCT GetLinePoint2d (L220-268).
fn get_line_point_2d(l: &IntPatchLine, param: f64, on_first: bool, u: &mut f64, v: &mut f64) {
    let nbptlin = l.nb_points();

    let par_trunc = param.trunc();
    let mut irang = par_trunc as usize;
    let par;
    if irang == nbptlin {
        irang -= 1;
        par = 1.0;
    } else {
        par = (param - par_trunc).abs();
    }

    let pnt = |i: usize| -> (f64, f64) {
        let w = l.point(i - 1);
        if on_first {
            (w.u1, w.v1)
        } else {
            (w.u2, w.v2)
        }
    };
    let (us1, vs1) = pnt(irang);
    let (us2, vs2) = pnt(irang + 1);

    *u = (1.0 - par) * us1 + par * us2;
    *v = (1.0 - par) * vs1 + par * vs2;
}

/// OCCT FindParameter (L269-421) — Walking branch; the Restriction branch
/// needs IntPatch_HInterTool::Project (not translated) and returns false —
/// the OCCT failure path.
fn find_parameter(
    l: &IntPatchLine,
    _tol: f64,
    ptsom: DVec3,
    _ptsom2d: DVec2,
    param: &mut f64,
    tgl: &mut DVec3,
    param_approche: usize,
) -> bool {
    // MSV 28.03.2002: find parameter on WLine in 2d space
    *tgl = DVec3::ZERO;
    if l.line_type != IntPatchIType::Walking {
        // OCCT L292-327: the Restriction branch — projects Ptsom2d on the
        // RLine's ArcOnS1/S2 (IntPatch_HInterTool::Project, not translated).
        // Returning false keeps the OCCT failure path.
        return false;
    }
    let tol2 = _tol * _tol;

    let nbpt = l.nb_points();
    let mut param_search_inf = 1usize;
    let mut param_search_sup = nbpt;

    if param_approche >= 3 {
        param_search_inf = param_approche - 2;
    }
    if param_approche + 2 < param_search_sup {
        param_search_sup = param_approche + 2;
    }

    let mut inf = [0usize; 3];
    let mut sup = [0usize; 3];
    // first search inside close bounding around ParamApproche;
    // then search to the nearest end of line;
    // and then search to the farthest end of line.
    inf[0] = param_search_inf;
    sup[0] = param_search_sup;
    if param_search_inf - 1 < nbpt - param_search_sup {
        inf[1] = 1;
        sup[1] = param_search_inf;
        inf[2] = param_search_sup;
        sup[2] = nbpt;
    } else {
        inf[1] = param_search_sup;
        sup[1] = nbpt;
        inf[2] = 1;
        sup[2] = param_search_inf;
    }

    for is in 0..3 {
        let mut p1 = l.point(inf[is] - 1).p3d;
        let mut v1 = p1 - ptsom;
        let mut norm1 = v1.length_squared();
        let mut normmin = tol2;
        let mut ibest = 0usize;
        if norm1 <= normmin {
            normmin = norm1;
            ibest = inf[is];
        }
        let mut found = false;
        for i in inf[is] + 1..=sup[is] {
            let p2 = l.point(i - 1).p3d;
            let v2 = p2 - ptsom;
            let norm2 = v2.length_squared();
            if v1.dot(v2) < 0.0 {
                *param = (i - 1) as f64 + 1.0 / (1.0 + (norm2 / norm1).sqrt());
                *tgl = p2 - p1;
                found = true;
                break;
            } else if norm2 < normmin {
                normmin = norm2;
                ibest = i;
            }
            v1 = v2;
            p1 = p2;
            norm1 = norm2;
        }
        if !found && ibest != 0 {
            *param = ibest as f64;
            found = true;
        }
        if found {
            return true;
        }
    }
    false
}

// ===========================================================================
// IntPatch_RstInt::PutVertexOnLine (IntPatch_RstInt.cxx L433-1355)
// ===========================================================================

/// Store the domain-vertex identity (OCCT Sommet.SetVertex(OnFirst, vtxarc)).
fn set_domain_vertex(v: &mut super::IntPatchVertex, on_first: bool, id: (u16, u16)) {
    if on_first {
        v.on_dom_s1 = true;
        v.is_vertex_on_s1 = true;
        v.vtx_on_s1 = Some(id);
    } else {
        v.on_dom_s2 = true;
        v.is_vertex_on_s2 = true;
        v.vtx_on_s2 = Some(id);
    }
}

/// Is the stored arc the same restriction (a 2D line matching the domain arc
/// frame within PConfusion)?
fn same_arc(a: &Curve2d, arc: &Line2dArc) -> bool {
    match a {
        Curve2d::Line(l) => {
            (l.origin - arc.origin).length() <= rcad_kernel::precision::PCONFUSION
                && (l.direction - arc.dir).length() <= rcad_kernel::precision::PCONFUSION
        }
        _ => false,
    }
}

fn diff_arc(a: &Curve2d, arc: &Line2dArc) -> bool {
    !same_arc(a, arc)
}

/// OCCT PutVertexOnLine(L, Surf, Domain, OtherSurf, OnFirst, Tol) — place
/// vertices on the Walking/RLine `line` where it crosses the domain
/// restrictions of `surf`; `other_surf` is the second intersection surface.
/// `surf_bounds` / `other_bounds` are the corrected UV bounds of the two
/// surfaces (the loaded adaptors' domains).
pub fn put_vertex_on_line(
    line: &mut IntPatchLine,
    surf: &Surface3,
    surf_bounds: [f64; 4],
    other_surf: &Surface3,
    other_bounds: [f64; 4],
    on_first: bool,
    tol: f64,
) {
    let mut param_arc: f64 = 0.0;

    // L474-477: the polyline tolerance from the surface resolutions.
    let mut tol_plin = u_resolution(surf, CONFUSION);
    tol_plin = tol_plin.max(v_resolution(surf, CONFUSION));
    tol_plin = tol_plin.min(CONFUSION);
    let _ = tol_plin; // consumed by PolyLine::new via INITDEFLE + Prepare

    let domain = RectDomain::from_bounds(
        surf_bounds[0],
        surf_bounds[1],
        surf_bounds[2],
        surf_bounds[3],
    );

    let typ_l = line.line_type;
    if typ_l != IntPatchIType::Walking && typ_l != IntPatchIType::Restriction {
        panic!("Standard_DomainError");
    }
    let mut plin = PolyLine::new(line, on_first);
    if !domain.has_3d() {
        // don't use computed deflection in the mode of pure geometric intersection
        plin.reset_error();
    }

    let surface_is_u_closed = SurfaceEval::is_u_closed(surf) || SurfaceEval::is_u_periodic(surf);
    let surface_is_v_closed = SurfaceEval::is_v_closed(surf) || SurfaceEval::is_v_periodic(surf);
    let o_surface_is_u_closed =
        SurfaceEval::is_u_closed(other_surf) || SurfaceEval::is_u_periodic(other_surf);
    let o_surface_is_v_closed =
        SurfaceEval::is_v_closed(other_surf) || SurfaceEval::is_v_periodic(other_surf);
    let possibly_closed =
        surface_is_u_closed || surface_is_v_closed || o_surface_is_u_closed || o_surface_is_v_closed;
    let (mut tol_u_closed, mut tol_v_closed, mut tol_ou_closed, mut tol_ov_closed) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    if possibly_closed {
        if surface_is_u_closed {
            tol_u_closed = (surf_bounds[1] - surf_bounds[0]) * 0.01;
        }
        if surface_is_v_closed {
            tol_v_closed = (surf_bounds[3] - surf_bounds[2]) * 0.01;
        }
        if o_surface_is_u_closed {
            tol_ou_closed = (other_bounds[1] - other_bounds[0]) * 0.01;
        }
        if o_surface_is_v_closed {
            tol_ov_closed = (other_bounds[3] - other_bounds[2]) * 0.01;
        }
    }

    //------------------------------------------------------------------------
    //-- On traite le cas ou la surface est periodique                      --
    //-- il faut dans ce cas considerer la restriction                      --
    //--                                la restriction decalee de +-2PI     --
    //------------------------------------------------------------------------
    let (surf1, surf2) = if on_first {
        (surf, other_surf)
    } else {
        (other_surf, surf)
    };
    let (surf1_bounds, surf2_bounds) = if on_first {
        (surf_bounds, other_bounds)
    } else {
        (other_bounds, surf_bounds)
    };
    let type_s1 = surface_type(surf1);
    let type_s2 = surface_type(surf2);
    let mut surface_is_periodic = false;
    let mut surface_is_bi_periodic = false;
    let surfacetype = if on_first { type_s1 } else { type_s2 };
    if surfacetype == SURF_CYLINDER
        || surfacetype == SURF_CONE
        || surfacetype == SURF_SPHERE
        || surfacetype == SURF_TORUS
    {
        surface_is_periodic = true;
        if surfacetype == SURF_TORUS {
            surface_is_bi_periodic = true;
        }
    }

    let mut numero_edge = 0usize;
    for arc_index in 0..domain.restrictions.len() {
        numero_edge += 1;
        let arc = domain.restrictions[arc_index].clone();

        // MSV Oct 15, 2001: use tolerance of this edge if possible
        let edge_tol = tol3d(tol);

        let (mut p_first, mut p_last) = arc.bounds();
        if rcad_kernel::precision::is_negative_infinite_value(p_first) {
            p_first = -TOPO_TOOL_INFINITE;
        }
        if rcad_kernel::precision::is_positive_infinite_value(p_last) {
            p_last = TOPO_TOOL_INFINITE;
        }

        let bplin = plin.bounding().clone();
        let mut offset_v = 0.0f64;
        let mut offset_u = 0.0f64;

        // OCCT L597-629: the GeomAbs_Line case (the rectangle restriction
        // arcs are always 2D lines) — rectangular axis-aligned domains only.
        let nb_echant: usize = {
            let a_lin = arc.line();
            let a_dir = a_lin.direction;
            let is_along_u = a_dir.x.abs() > a_dir.y.abs();
            if surface_is_periodic && !is_along_u {
                // Shift along U-direction
                if let Some((xmin, _, _, _)) = bplin.get() {
                    let a_xmin = xmin;
                    let a_new_location =
                        in_period(a_lin.origin.x, a_xmin, a_xmin + std::f64::consts::PI * 2.0);
                    offset_u = a_new_location - a_lin.origin.x;
                }
            } else if surface_is_bi_periodic && is_along_u {
                // Shift along V-direction
                if let Some((_, ymin, _, _)) = bplin.get() {
                    let a_ymin = ymin;
                    let a_new_location =
                        in_period(a_lin.origin.y, a_ymin, a_ymin + std::f64::consts::PI * 2.0);
                    offset_v = a_new_location - a_lin.origin.y;
                }
            }
            10
        };

        let mut brise = PolyArc::new(&arc, nb_echant, p_first, p_last, &bplin);

        let mut indice_offset_bi_periodic = 0i32;
        let mut indice_offset_periodic = 0i32;
        let a_ref_ou = offset_u;
        let a_ref_ov = offset_v;

        loop {
            if indice_offset_bi_periodic == 1 {
                offset_v = a_ref_ov - std::f64::consts::PI - std::f64::consts::PI;
            } else if indice_offset_bi_periodic == 2 {
                offset_v = a_ref_ov + std::f64::consts::PI + std::f64::consts::PI;
            }

            loop {
                if indice_offset_periodic == 1 {
                    offset_u = a_ref_ou - std::f64::consts::PI - std::f64::consts::PI;
                } else if indice_offset_periodic == 2 {
                    offset_u = a_ref_ou + std::f64::consts::PI + std::f64::consts::PI;
                }

                brise.set_offset(offset_u, offset_v);

                let commun = InterferencePolygon2d::new_polygon_polygon(&plin, &brise);
                let mut locpt: Vec<DVec3> = Vec::new();
                let mut locpt2: Vec<DVec2> = Vec::new();

                // We do not need in putting vertex into tangent zone(s).
                // Therefore, only section points are interested by us.
                // Boundary of WLine (its first/last points) will be
                // marked by some vertex later. See bug #29494.
                let a_nb_section_pts = commun.interf.nb_section_points();
                let mut nb_treated = 0usize;
                for i in 1..=a_nb_section_pts {
                    let sp = commun.interf.pnt_value(i);
                    let a_w1 = sp.param_on_first();
                    let a_w2 = sp.param_on_second();

                    let mut u = 0.0f64;
                    let mut v = 0.0f64;
                    get_line_point_2d(line, a_w1 + 1.0, !on_first, &mut u, &mut v);

                    let par_trunc = a_w2.trunc();
                    let mut irang = par_trunc as usize + 1;
                    let par;
                    if irang == brise.nb_points() {
                        irang -= 1;
                        par = 1.0;
                    } else {
                        par = (a_w2 - par_trunc).abs();
                    }

                    let w = (1.0 - par) * brise.parameter(irang) + par * brise.parameter(irang + 1);

                    //------------------------------------------------------------------------
                    //-- On a trouve un point 2d approche Ua,Va  intersection de la ligne
                    //-- de cheminement et de la restriction.
                    //--
                    //-- On injecte ce point ds les intersections Courbe-Surface
                    //--
                    let other_bounded = BoundedSurface {
                        surf: other_surf,
                        u0: other_bounds[0],
                        u1: other_bounds[1],
                        v0: other_bounds[2],
                        v1: other_bounds[3],
                    };
                    // MSV: extend UV bounds to not miss solution near the boundary
                    let marg_coef = 0.004;
                    let mut refined = false;
                    let mut ptsommet = DVec3::ZERO;
                    let mut u2 = 0.0f64;
                    let mut v2 = 0.0f64;
                    {
                        let mut thefunc = CSFunction::new(&other_bounded, &arc, surf);
                        let mut int_cs =
                            IntCS::<BoundedSurface, Line2dArc, BoundedSurfaceTool, Line2dArcTool, CSFunction>::with_margin(
                                u, v, w, &mut thefunc, edge_tol, marg_coef,
                            );
                        if int_cs.is_done() && !int_cs.is_empty() {
                            ptsommet = int_cs.point();
                            let (uu, vv) = int_cs.parameter_on_surface();
                            u2 = uu;
                            v2 = vv;
                            let an_old_pnt = SurfaceEval::point_at(other_surf, u, v);
                            let a_new_pnt = SurfaceEval::point_at(other_surf, u2, v2);
                            // if (anOldPnt.SquareDistance(aNewPnt) < Precision::SquareConfusion())
                            let a_tol_conf = CONFUSION.max(edge_tol);

                            if an_old_pnt.distance_squared(a_new_pnt) < a_tol_conf * a_tol_conf {
                                u2 = u;
                                v2 = v;
                            }
                            param_arc = int_cs.parameter_on_curve();
                            refined = true;
                        }
                    }

                    if refined {
                        let mut duplicate = false;
                        for j in 0..locpt.len() {
                            if ptsommet.distance(locpt[j]) <= edge_tol {
                                if possibly_closed {
                                    let (lu, lv) = (locpt2[j].x, locpt2[j].y);
                                    if (o_surface_is_u_closed && (lu - u2).abs() > tol_ou_closed)
                                        || (o_surface_is_v_closed && (lv - v2).abs() > tol_ov_closed)
                                    {
                                        continue;
                                    }
                                }
                                duplicate = true;
                                break;
                            }
                        }

                        if !duplicate {
                            let param_approx_on_line = (a_w1 as i64 + 1) as usize;

                            let p2d = arc.value(param_arc);
                            let d2d = arc.d1();
                            let mut u1 = p2d.x;
                            let mut v1 = p2d.y;
                            if typ_l == IntPatchIType::Walking && surface_is_periodic {
                                recadre(
                                    type_s1,
                                    type_s2,
                                    line,
                                    param_approx_on_line,
                                    &mut u1,
                                    &mut v1,
                                    &mut u2,
                                    &mut v2,
                                );
                            }
                            locpt.push(ptsommet);
                            locpt2.push(DVec2::new(u2, v2));

                            let mut tgl = DVec3::ZERO;
                            let mut param_line: f64 = 0.0;
                            let mut found = find_parameter(
                                line,
                                edge_tol,
                                ptsommet,
                                DVec2::new(u2, v2),
                                &mut param_line,
                                &mut tgl,
                                param_approx_on_line,
                            );

                            if typ_l == IntPatchIType::Walking && found && possibly_closed {
                                // check in 2d
                                if surface_is_u_closed || surface_is_v_closed {
                                    let mut cu = 0.0;
                                    let mut cv = 0.0;
                                    get_line_point_2d(line, param_line, on_first, &mut cu, &mut cv);
                                    if (surface_is_u_closed && (cu - u1).abs() > tol_u_closed)
                                        || (surface_is_v_closed && (cv - v1).abs() > tol_v_closed)
                                    {
                                        found = false;
                                    }
                                }
                                if found && (o_surface_is_u_closed || o_surface_is_v_closed) {
                                    let mut cu = 0.0;
                                    let mut cv = 0.0;
                                    get_line_point_2d(line, param_line, !on_first, &mut cu, &mut cv);
                                    if (o_surface_is_u_closed && (cu - u2).abs() > tol_ou_closed)
                                        || (o_surface_is_v_closed && (cv - v2).abs() > tol_ov_closed)
                                    {
                                        found = false;
                                    }
                                }
                            }
                            if !found {
                                continue;
                            }

                            let mut vtx_on_arc =
                                coincide_on_arc(ptsommet, &arc, arc_index, surf, edge_tol, &domain);
                            let mut vtx_tol;
                            if vtx_on_arc.is_some() {
                                vtx_tol = tol3d(edge_tol);
                                if edge_tol > vtx_tol {
                                    vtx_tol = edge_tol;
                                }
                            } else {
                                vtx_tol = edge_tol;
                            }

                            //-- It is necessary to test that the point does not already exist
                            //--   - It can be already a point on arc
                            //--        BUT on a different arc
                            // MSV 27.03.2002: find the nearest point; add check in 2d
                            let mut ivtx = 0usize;
                            let mut ptline: Option<super::IntPatchVertex> = None;
                            let mut dmin = f64::MAX;
                            for j in 1..=line.nb_vertex() {
                                let rptline = line.vertex(j);
                                let a_point_on_rst_still_exist = (on_first
                                    && rptline.on_dom_s1
                                    && rptline.arc_on_s1.as_ref().is_some_and(|a| same_arc(a, &arc)))
                                    || (!on_first
                                        && rptline.on_dom_s2
                                        && rptline.arc_on_s2.as_ref().is_some_and(|a| same_arc(a, &arc)));
                                if !a_point_on_rst_still_exist {
                                    if possibly_closed {
                                        if surface_is_u_closed || surface_is_v_closed {
                                            let (ru, rv) = if on_first {
                                                (rptline.u1, rptline.v1)
                                            } else {
                                                (rptline.u2, rptline.v2)
                                            };
                                            if (surface_is_u_closed && (ru - u1).abs() > tol_u_closed)
                                                || (surface_is_v_closed && (rv - v1).abs() > tol_v_closed)
                                            {
                                                continue;
                                            }
                                        }
                                        if o_surface_is_u_closed || o_surface_is_v_closed {
                                            let (ru, rv) = if on_first {
                                                (rptline.u2, rptline.v2)
                                            } else {
                                                (rptline.u1, rptline.v1)
                                            };
                                            if (o_surface_is_u_closed && (ru - u2).abs() > tol_ou_closed)
                                                || (o_surface_is_v_closed && (rv - v2).abs() > tol_ov_closed)
                                            {
                                                continue;
                                            }
                                        }
                                    }
                                    let dist = ptsommet.distance(rptline.p3d);
                                    let dt = vtx_tol.max(rptline.tolerance);
                                    if dist < dmin {
                                        if dist <= dt {
                                            ptline = Some(rptline.clone());
                                            ivtx = j;
                                            if surfacetype == SURF_CONE {
                                                ivtx = 0;
                                            }
                                        } else {
                                            // cancel previous solution because this point is better
                                            // but its tolerance is not large enough
                                            ivtx = 0;
                                        }
                                        dmin = dist;
                                    }
                                }
                            }
                            if ivtx != 0 {
                                let pl = ptline.as_ref().unwrap();
                                if pl.tolerance > vtx_tol {
                                    vtx_tol = pl.tolerance;
                                    if vtx_on_arc.is_none() {
                                        // now we should repeat attempt to coincide on a bound of arc
                                        vtx_on_arc = coincide_on_arc(
                                            ptsommet,
                                            &arc,
                                            arc_index,
                                            surf,
                                            vtx_tol,
                                            &domain,
                                        );
                                        if vtx_on_arc.is_some() {
                                            let t = tol3d(edge_tol);
                                            if t > vtx_tol {
                                                vtx_tol = t;
                                            }
                                        }
                                    }
                                }
                            }

                            if typ_l == IntPatchIType::Walking {
                                verify_tgline(line, param_line as usize, edge_tol, &mut tgl);
                            }

                            let (_ptbid, d1u, d1v) = SurfaceEval::derivatives(surf, u1, v1);
                            let tgrst = d2d.x * d1u + d2d.y * d1v;

                            let normsurf = d1u.cross(d1v);
                            let mut transline = Transition::new();
                            let mut transarc = Transition::new();
                            if normsurf.length() < f64::MIN_POSITIVE {
                                transline.set_value_in_out(true, transitions::TypeTrans::Undecided);
                                transarc.set_value_in_out(true, transitions::TypeTrans::Undecided);
                            } else {
                                transitions::make_transition(
                                    tgl,
                                    tgrst,
                                    normsurf,
                                    &mut transline,
                                    &mut transarc,
                                );
                            }

                            nb_treated += 1;
                            if ivtx == 0 {
                                let mut sommet = super::IntPatchVertex::default();
                                sommet.set_value(ptsommet, vtx_tol, false); // pour tangence
                                if on_first {
                                    sommet.set_parameters(u1, v1, u2, v2);
                                } else {
                                    sommet.set_parameters(u2, v2, u1, v1);
                                }

                                if let Some(v) = &vtx_on_arc {
                                    set_domain_vertex(&mut sommet, on_first, v.id);
                                }

                                //---------------------------------------------------------
                                //-- lbr : On remplace le point d indice paramline sur la -
                                //-- ligne par le vertex .                                -
                                //---------------------------------------------------------
                                sommet.set_parameter(param_line); // sur ligne d intersection
                                sommet.set_arc(
                                    on_first,
                                    Curve2d::Line(arc.line()),
                                    param_arc,
                                    transline,
                                    transarc,
                                );

                                if typ_l == IntPatchIType::Walking {
                                    line.add_vertex(sommet);
                                } else {
                                    line.add_vertex(sommet);
                                }
                            } else {
                                // CAS DE FIGURE : en appelant s1 la surf sur laquelle on
                                //   connait les pts sur restriction, et s2 celle sur laquelle
                                //   on les cherche. Le point trouve verifie necessairement
                                //   IsOnDomS1 = True.
                                let pl = ptline.as_mut().unwrap();
                                let on_different_rst = (on_first
                                    && pl.on_dom_s1
                                    && pl.arc_on_s1.as_ref().is_some_and(|a| diff_arc(a, &arc)))
                                    || (!on_first
                                        && pl.on_dom_s2
                                        && pl.arc_on_s2.as_ref().is_some_and(|a| diff_arc(a, &arc)));
                                pl.set_tolerance(vtx_tol);
                                if (!pl.is_vertex_on_s1 && on_first)
                                    || (!pl.is_vertex_on_s2 && !on_first)
                                    || on_different_rst
                                {
                                    if (!pl.on_dom_s2 && !on_first)
                                        || (!pl.on_dom_s1 && on_first)
                                        || on_different_rst
                                    {
                                        pl.set_arc(
                                            on_first,
                                            Curve2d::Line(arc.line()),
                                            param_arc,
                                            transline,
                                            transarc,
                                        );
                                        if let Some(v) = &vtx_on_arc {
                                            set_domain_vertex(pl, on_first, v.id);
                                        }
                                        if typ_l == IntPatchIType::Walking {
                                            if on_different_rst {
                                                let pl = pl.clone();
                                                line.add_vertex(pl);
                                            } else {
                                                let pl = pl.clone();
                                                line.replace_vertex(ivtx, pl);
                                            }
                                        } else if on_different_rst {
                                            let pl = pl.clone();
                                            line.add_vertex(pl);
                                        } else {
                                            let pl = pl.clone();
                                            line.replace_vertex(ivtx, pl);
                                        }
                                    } else if (on_first && pl.is_vertex_on_s2)
                                        || (!on_first && pl.is_vertex_on_s1)
                                    {
                                        let mut sommet = pl.clone();
                                        sommet.set_arc(
                                            on_first,
                                            Curve2d::Line(arc.line()),
                                            param_arc,
                                            transline,
                                            transarc,
                                        );
                                        if let Some(v) = &vtx_on_arc {
                                            set_domain_vertex(&mut sommet, on_first, v.id);
                                        }
                                        line.add_vertex(sommet);
                                    }
                                } else {
                                    let vtxref_id = if on_first {
                                        pl.vtx_on_s1
                                    } else {
                                        pl.vtx_on_s2
                                    };
                                    if (on_first && !pl.on_dom_s2) || (!on_first && !pl.on_dom_s1)
                                    {
                                        pl.set_arc(
                                            on_first,
                                            Curve2d::Line(arc.line()),
                                            param_arc,
                                            transline,
                                            transarc,
                                        );
                                        if let Some(v) = &vtx_on_arc {
                                            set_domain_vertex(pl, on_first, v.id);
                                        }
                                        let pl = pl.clone();
                                        line.replace_vertex(ivtx, pl);

                                        for k in 1..=line.nb_vertex() {
                                            if k != ivtx {
                                                let mut ptline_k = line.vertex(k).clone();
                                                if (on_first && ptline_k.is_vertex_on_s1)
                                                    || (!on_first && ptline_k.is_vertex_on_s2)
                                                {
                                                    let kid = if on_first {
                                                        ptline_k.vtx_on_s1
                                                    } else {
                                                        ptline_k.vtx_on_s2
                                                    };
                                                    if vtxref_id.is_some() && kid == vtxref_id {
                                                        if ptline_k.tolerance < vtx_tol {
                                                            ptline_k.set_tolerance(vtx_tol);
                                                        }
                                                        ptline_k.set_arc(
                                                            on_first,
                                                            Curve2d::Line(arc.line()),
                                                            param_arc,
                                                            transline,
                                                            transarc,
                                                        );
                                                        if let Some(v) = &vtx_on_arc {
                                                            set_domain_vertex(
                                                                &mut ptline_k,
                                                                on_first,
                                                                v.id,
                                                            );
                                                        }
                                                        line.replace_vertex(k, ptline_k);
                                                    }
                                                }
                                            }
                                        }
                                    } else if (on_first && pl.is_vertex_on_s2)
                                        || (!on_first && pl.is_vertex_on_s1)
                                    {
                                        //                on doit avoir vtxons2 = vtxarc... pas de verif...
                                        let mut sommet = pl.clone();
                                        sommet.set_arc(
                                            on_first,
                                            Curve2d::Line(arc.line()),
                                            param_arc,
                                            transline,
                                            transarc,
                                        );
                                        for k in 1..=line.nb_vertex() {
                                            if k != ivtx {
                                                let mut ptline_k = line.vertex(k).clone();
                                                if (on_first && ptline_k.is_vertex_on_s1)
                                                    || (!on_first && ptline_k.is_vertex_on_s2)
                                                {
                                                    let kid = if on_first {
                                                        ptline_k.vtx_on_s1
                                                    } else {
                                                        ptline_k.vtx_on_s2
                                                    };
                                                    if vtxref_id.is_some() && kid == vtxref_id {
                                                        if ptline_k.tolerance < vtx_tol {
                                                            ptline_k.set_tolerance(vtx_tol);
                                                        }
                                                        ptline_k.set_arc(
                                                            on_first,
                                                            Curve2d::Line(arc.line()),
                                                            param_arc,
                                                            transline,
                                                            transarc,
                                                        );
                                                        line.replace_vertex(k, ptline_k);
                                                        line.add_vertex(sommet.clone());
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if nb_treated == 2 && typ_l == IntPatchIType::Walking {
                    // We processed a tangent zone, and both ends have been treated.
                    // So mark WLine as having arc
                    if on_first {
                        line.set_arc_on_s1(Curve2d::Line(arc.line()));
                    } else {
                        line.set_arc_on_s2(Curve2d::Line(arc.line()));
                    }
                }

                indice_offset_periodic += 1;
                if !(surface_is_periodic && indice_offset_periodic <= 2) {
                    break;
                }
            }

            indice_offset_bi_periodic += 1;
            if !(surface_is_bi_periodic && indice_offset_bi_periodic <= 2) {
                break;
            }
        }
        let _ = numero_edge;
    }

    //--------------------------------------------------------------------------------
    //-- On reprend la ligne et on recale les parametres des vertex.
    //--
    if typ_l == IntPatchIType::Walking {
        let two_pi = std::f64::consts::PI + std::f64::consts::PI;
        let mut pu1 = 0.0f64;
        let mut pv1 = 0.0f64;
        let mut pu2 = 0.0f64;
        let mut pv2 = 0.0f64;
        match type_s1 {
            SURF_CYLINDER | SURF_CONE | SURF_SPHERE => pu1 = two_pi,
            SURF_TORUS => {
                pu1 = two_pi;
                pv1 = two_pi;
            }
            _ => {
                if SurfaceEval::is_u_periodic(surf1) {
                    pu1 = two_pi;
                } else if SurfaceEval::is_u_closed(surf1) {
                    pu1 = surf1_bounds[1] - surf1_bounds[0];
                }
                if SurfaceEval::is_v_periodic(surf1) {
                    pv1 = two_pi;
                } else if SurfaceEval::is_v_closed(surf1) {
                    pv1 = surf1_bounds[3] - surf1_bounds[2];
                }
            }
        }
        match type_s2 {
            SURF_CYLINDER | SURF_CONE | SURF_SPHERE => pu2 = two_pi,
            SURF_TORUS => {
                pu2 = two_pi;
                pv2 = two_pi;
            }
            _ => {
                if SurfaceEval::is_u_periodic(surf2) {
                    pu2 = two_pi;
                } else if SurfaceEval::is_u_closed(surf2) {
                    pu2 = surf2_bounds[1] - surf2_bounds[0];
                }
                if SurfaceEval::is_v_periodic(surf2) {
                    pv2 = two_pi;
                } else if SurfaceEval::is_v_closed(surf2) {
                    pv2 = surf2_bounds[3] - surf2_bounds[2];
                }
            }
        }

        line.compute_vertex_parameters_wline(tol, [pu1, pv1, pu2, pv2]);
    } else {
        line.compute_vertex_parameters_rline();
    }
}

/// OCCT ElCLib::InPeriod(U, UFirst, ULast).
fn in_period(u: f64, u_first: f64, u_last: f64) -> f64 {
    let period = u_last - u_first;
    if period <= 0.0 {
        return u;
    }
    let mut r = (u - u_first) % period;
    if r < 0.0 {
        r += period;
    }
    u_first + r
}
