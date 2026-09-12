//! OCCT BRepLib_MakeFace (TKTopAlgo/BRepLib/BRepLib_MakeFace.cxx L1-925,
//! whole file) — 1:1 translation.
//!
//! Architecture adaptations (rcad kernel):
//! - `BRep_Builder::MakeFace(F, S, L, tol)` -> the kernel pool
//!   `BRep::add_tface_tol` (surface, no wires, tolerance); the OCCT
//!   location argument is stored into `TFaceData::surface_location`.
//! - `BRep_Builder::NaturalRestriction(F, b)` -> the
//!   `TFaceData::natural_restriction` field.
//! - `BRep_Builder::Add(aFace, aWire)` -> the outer wire slot for the first
//!   wire and the inner-wire list for the following ones.
//! - `Geom_Surface::Bounds` -> `SurfaceEval::default_domain`.
//! - `BRepLib::UpdateTolerances(myShape)` / `BRepLib::SameParameter(myShape,
//!   tol, WithPCurve)` (the shape-level statics) -> GAP leaves recorded as
//!   no-ops (the shape-level forms are pending; the edge-level SameParameter
//!   home is brep_lib.rs).
//! - `Geom_Surface::UIso/VIso` -> GAP leaf raising exactly where OCCT calls
//!   the untranslated kernel iso extraction (BRepLib_MakeFace.cxx
//!   L575/L580/L585/L590); the bounded-surface edge assembly below stays 1:1
//!   and becomes reachable when the kernel leaf lands.  The consumer anchor
//!   MakeFace(gp_Pln) passes through the all-infinite flags and never
//!   reaches the iso calls.

use rcad_kernel::core::precision::{
    is_infinite_value, is_negative_infinite_value, is_positive_infinite_value, COMPUTATIONAL,
    CONFUSION,
};
use rcad_kernel::geom::{
    BSplineCurve3, BezierCurve3, Circle3, Curve2d, Curve3, Line2d, Line3, Surface3, SurfaceEval,
    TrimmedSurface,
};
use rcad_kernel::topods::State;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, Shape, TShape};

use glam::DVec2;

use crate::topalgo::brep_lib::make_wire::{
    builder_add_edge_vertex, set_shape_closed, shape_closed, wire_edges,
};
use crate::topalgo::brep_lib_find_surface::BRepLibFindSurface;
use crate::topalgo::brep_top_adaptor::fclass2d::FClass2d;
use crate::topalgo::shape_source::FaceShapeSource;

// =========================================================================
// OCCT BRepLib_FaceError (BRepLib_MakeFace.hxx L24-30).
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepLibFaceError {
    NoFace,
    NotPlanar,
    FaceDone,
    ParametersOutOfRange,
}

// =========================================================================
// OCCT BRepLib_MakeFace (BRepLib_MakeFace.hxx L36-160).
// =========================================================================
pub struct BRepLibMakeFace {
    /// hxx: TopoDS_Shape myShape.
    my_shape: Shape,
    /// hxx: BRepLib_FaceError myError.
    my_error: BRepLibFaceError,
    /// OCCT BRepLib_MakeShape Done()/NotDone() flag.
    done: bool,
}

impl Default for BRepLibMakeFace {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepLibMakeFace {
    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace() (L56-59).
    pub fn new() -> Self {
        BRepLibMakeFace {
            my_shape: Shape::null(),
            my_error: BRepLibFaceError::NoFace,
            done: false,
        }
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const TopoDS_Face& F)
    /// (L63-66).
    pub fn new_with_face(brep: &mut BRep, bb: &mut BRepBuilder, f: &Shape) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_face(brep, bb, f);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const gp_Pln& P) (L70-74).
    pub fn new_with_plane(brep: &mut BRep, bb: &mut BRepBuilder, p: &rcad_kernel::geom::Plane) -> Self {
        // L72: Geom_Plane GP = new Geom_Plane(P).
        let gp = Surface3::Plane(p.clone());
        // L73: Init(GP, true, Precision::Confusion()).
        let mut r = BRepLibMakeFace::new();
        r.init_bound(brep, bb, &gp, true, CONFUSION);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const gp_Cylinder& C)
    /// (L78-82) — the kernel surface payload carries the Geom_ surface.
    pub fn new_with_cylinder(brep: &mut BRep, bb: &mut BRepBuilder, c: &Surface3) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_bound(brep, bb, c, true, CONFUSION);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const gp_Cone& C) (L86-90).
    pub fn new_with_cone(brep: &mut BRep, bb: &mut BRepBuilder, c: &Surface3) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_bound(brep, bb, c, true, CONFUSION);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const gp_Sphere& S) (L94-98).
    pub fn new_with_sphere(brep: &mut BRep, bb: &mut BRepBuilder, s: &Surface3) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_bound(brep, bb, s, true, CONFUSION);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const gp_Torus& T) (L102-106).
    pub fn new_with_torus(brep: &mut BRep, bb: &mut BRepBuilder, t: &Surface3) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_bound(brep, bb, t, true, CONFUSION);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const handle(Geom_Surface)& S,
    /// const double TolDegen) (L110-113).
    pub fn new_with_surface_tol(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        s: &Surface3,
        tol_degen: f64,
    ) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_bound(brep, bb, s, true, tol_degen);
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const TopoDS_Wire& W, const
    /// bool OnlyPlane) (L189-258).
    pub fn new_with_wire(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        w: &Shape,
        only_plane: bool,
    ) -> Self {
        let mut r = BRepLibMakeFace::new();
        // L193: Find a surface through the wire — BRepLib_FindSurface FS(W,
        // -1, OnlyPlane, true) (the OnlyClosed argument is true).
        let fs = BRepLibFindSurface::new_closed(brep, w, -1.0, only_plane, true);
        // L194-198.
        if !fs.found() {
            r.my_error = BRepLibFaceError::NotPlanar;
            return r;
        }

        // L200-203: build the face and add the wire.
        r.my_error = BRepLibFaceError::FaceDone;
        r.done = true;

        // L204: double tol = max(1.2 * FS.ToleranceReached(), FS.Tolerance()).
        let tol = (1.2 * fs.tolerance_reached()).max(fs.tolerance());

        // L206: B.MakeFace(TopoDS::Face(myShape), FS.Surface(), FS.Location(),
        // tol).
        let fsurf = fs.surface().expect("BRepLib_FindSurface::Surface");
        let f = brep.add_tface_tol(
            Some(fsurf),
            Shape::null(),
            Vec::new(),
            None,
            None,
            Vec::new(),
            false,
            tol,
        );
        r.my_shape = f;
        brep.face_mut(r.my_shape.clone()).surface_location = fs.location();

        // L208-246: aW — get rid of degenerative edges in the input wire
        // when OnlyPlane.
        let a_w = if only_plane {
            let mut aw = bb.make_wire(brep);
            let mut has_degenerated = false;
            // L215-218: aWForw = W; aWForw.Orientation(TopAbs_FORWARD).
            let mut a_w_forw = w.clone();
            a_w_forw.orientation = Orientation::Forward;
            // L219-231.
            for e in wire_edges(brep, &a_w_forw) {
                if brep.edge(e.clone()).degenerated {
                    // L223-226.
                    has_degenerated = true;
                } else {
                    // L229: aB.Add(aW, aE).
                    bb.add_to_wire(brep, aw.clone(), e.clone());
                }
            }

            if has_degenerated {
                // L233-237: return to original orient.
                aw.orientation = w.orientation;
                set_shape_closed(brep, &aw, shape_closed(brep, w));
                aw
            } else {
                // L239-241.
                w.clone()
            }
        } else {
            // L243-246.
            w.clone()
        };

        // L248: Add(aW).
        r.add_wire(brep, bb, &a_w);
        // L250: BRepLib::UpdateTolerances(myShape) — GAP leaf, see the module
        // header.
        brep_lib_update_tolerances_gap(&r.my_shape);
        // L252: BRepLib::SameParameter(myShape, tol, true) — GAP leaf, see
        // the module header.
        brep_lib_same_parameter_shape_gap(&r.my_shape, tol, true);
        // L254-257.
        if shape_closed(brep, &a_w) {
            r.check_inside(brep);
        }
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const gp_Pln& P, const
    /// TopoDS_Wire& W, const bool Inside) (L262-271).
    pub fn new_with_plane_and_wire(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        p: &rcad_kernel::geom::Plane,
        w: &Shape,
        inside: bool,
    ) -> Self {
        // L264: Geom_Plane Pl = new Geom_Plane(P).
        let pl = Surface3::Plane(p.clone());
        let mut r = BRepLibMakeFace::new();
        // L265: Init(Pl, false, Precision::Confusion()).
        r.init_bound(brep, bb, &pl, false, CONFUSION);
        // L266: Add(W).
        r.add_wire(brep, bb, w);
        // L267-270.
        if inside && shape_closed(brep, w) {
            r.check_inside(brep);
        }
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const handle(Geom_Surface)& S,
    /// const TopoDS_Wire& W, const bool Inside) (L327-337) — the general
    /// surface form; the cylinder/cone/sphere/torus+wire ctors (L275-323)
    /// follow the same three steps with their surface payloads.
    pub fn new_with_surface_and_wire(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        s: &Surface3,
        w: &Shape,
        inside: bool,
    ) -> Self {
        let mut r = BRepLibMakeFace::new();
        // L331: Init(S, false, Precision::Confusion()).
        r.init_bound(brep, bb, s, false, CONFUSION);
        // L332: Add(W).
        r.add_wire(brep, bb, w);
        // L333-336.
        if inside && shape_closed(brep, w) {
            r.check_inside(brep);
        }
        r
    }

    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const TopoDS_Face& F, const
    /// TopoDS_Wire& W) (L341-345).
    pub fn new_with_face_and_wire(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        f: &Shape,
        w: &Shape,
    ) -> Self {
        let mut r = BRepLibMakeFace::new();
        r.init_face(brep, bb, f);
        r.add_wire(brep, bb, w);
        r
    }

    /// OCCT BRepLib_MakeFace::Init(const TopoDS_Face& F) (L349-362).
    pub fn init_face(&mut self, brep: &mut BRep, _bb: &mut BRepBuilder, f: &Shape) {
        // L352: copy the face — myShape = F.EmptyCopied().
        self.my_shape = brep.empty_copied(f);
        // L353.
        self.my_error = BRepLibFaceError::FaceDone;
        self.done = true;

        // L355-361: BRep_Builder B; TopoDS_Iterator It(F); while (It.More())
        // B.Add(myShape, It.Value()) — the kernel empty_copied already
        // carries the child slots (outer wire / inner wires / internal
        // vertices), so the add loop is the copy itself.
    }

    /// OCCT BRepLib_MakeFace::Init(const handle(Geom_Surface)& S, const bool
    /// Bound, const double TolDegen) (L366-384).
    pub fn init_bound(
        &mut self,
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        s: &Surface3,
        bound: bool,
        tol_degen: f64,
    ) {
        // L370.
        self.my_error = BRepLibFaceError::FaceDone;
        self.done = true;
        if bound {
            // L373-375: S->Bounds(UMin, UMax, VMin, VMax); Init(...).
            let [umin, umax, vmin, vmax] = s.default_domain();
            self.init_bounds(brep, bb, s, umin, umax, vmin, vmax, tol_degen);
        } else {
            // L379-381: BRep_Builder B; B.MakeFace(myShape, S,
            // Precision::Confusion()).
            self.my_shape = brep.add_tface_tol(
                Some(s.clone()),
                Shape::null(),
                Vec::new(),
                None,
                None,
                Vec::new(),
                false,
                CONFUSION,
            );
        }
        // L382-383: B.NaturalRestriction(TopoDS::Face(myShape), true).
        brep.face_mut(self.my_shape.clone()).natural_restriction = true;
    }

    /// OCCT BRepLib_MakeFace::IsDegenerated (L392-459) — checks whether the
    /// passed curve is degenerated with the passed tolerance value.
    pub fn is_degenerated(
        the_curve: Option<&Curve3>,
        the_max_tol: f64,
        the_act_tol: &mut f64,
    ) -> bool {
        // L396-398.
        let a_confusion = CONFUSION;
        *the_act_tol = a_confusion;
        let curve = match the_curve {
            Some(c) => c,
            None => return false,
        };
        match curve {
            // L401-410: GeomAbs_Circle.
            Curve3::Circle(Circle3 { radius, .. }) => {
                if *radius > the_max_tol {
                    return false;
                }
                *the_act_tol = (*radius).max(a_confusion);
                true
            }
            // L411-433: GeomAbs_BSplineCurve.
            Curve3::BSpline(BSplineCurve3 { control_points, .. }) => {
                let nb_poles = control_points.len();
                let mut a_max_pole_dist2 = 0.0f64;
                let a_max_tol2 = the_max_tol * the_max_tol;
                let p1 = control_points[0];
                for p2 in control_points.iter().take(nb_poles).skip(1) {
                    let a_pole_dist2 = p1.distance_squared(*p2);
                    if a_pole_dist2 > a_max_tol2 {
                        return false;
                    }
                    if a_pole_dist2 > a_max_pole_dist2 {
                        a_max_pole_dist2 = a_pole_dist2;
                    }
                }
                *the_act_tol = (1.000001 * a_max_pole_dist2.sqrt()).max(a_confusion);
                true
            }
            // L434-456: GeomAbs_BezierCurve.
            Curve3::Bezier(BezierCurve3 { control_points, .. }) => {
                let nb_poles = control_points.len();
                let mut a_max_pole_dist2 = 0.0f64;
                let a_max_tol2 = the_max_tol * the_max_tol;
                let p1 = control_points[0];
                for p2 in control_points.iter().take(nb_poles).skip(1) {
                    let a_pole_dist2 = p1.distance_squared(*p2);
                    if a_pole_dist2 > a_max_tol2 {
                        return false;
                    }
                    if a_pole_dist2 > a_max_pole_dist2 {
                        a_max_pole_dist2 = a_pole_dist2;
                    }
                }
                *the_act_tol = (1.000001 * a_max_pole_dist2.sqrt()).max(a_confusion);
                true
            }
            // L458.
            _ => false,
        }
    }

    /// OCCT BRepLib_MakeFace::Init(const handle(Geom_Surface)& SS, const
    /// double Um, UM, Vm, VM, const double TolDegen) (L463-869).
    #[allow(clippy::too_many_arguments)]
    pub fn init_bounds(
        &mut self,
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        ss: &Surface3,
        um: f64,
        um_m: f64,
        vm: f64,
        vm_m: f64,
        tol_degen: f64,
    ) {
        // L470.
        self.my_error = BRepLibFaceError::FaceDone;
        self.done = true;

        // L472-475.
        let mut u_min = um;
        let mut u_max = um_m;
        let mut v_min = vm;
        let mut v_max = vm_m;

        // L477: double umin, umax, vmin, vmax, T.
        let mut t;

        // L479-485: S = SS, BS = SS; the RectangularTrimmed unwrap.
        let mut s = ss.clone();
        let bs: Surface3 = match ss {
            Surface3::Trimmed(tr) => (*tr.basis).clone(),
            other => other.clone(),
        };

        // L487.
        let offset_surface = matches!(bs, Surface3::Offset(_));

        // L491: constexpr double epsilon = Precision::PConfusion().
        let epsilon = CONFUSION * 0.01;

        // L493: BS->Bounds(umin, umax, vmin, vmax).
        let [umin, umax, vmin, vmax] = bs.default_domain();

        // L495-518: the offset-surface re-trimming.
        if offset_surface {
            if let Surface3::Offset(os) = &bs {
                let base: Surface3 = (*os.basis).clone();
                if matches!(base, Surface3::LinearExtrusion(_)) {
                    // L500-510.
                    if is_infinite_value(umin) || is_infinite_value(umax) {
                        s = Surface3::Trimmed(TrimmedSurface::new(
                            bs.clone(),
                            u_min,
                            u_max,
                            v_min,
                            v_max,
                        ));
                    } else {
                        // L508: the V-only trim keeps the natural U bounds.
                        s = Surface3::Trimmed(TrimmedSurface::new(
                            bs.clone(),
                            umin,
                            umax,
                            v_min,
                            v_max,
                        ));
                    }
                } else if matches!(base, Surface3::Revolution(_)) {
                    // L511-517.
                    if is_infinite_value(vmin) || is_infinite_value(vmax) {
                        s = Surface3::Trimmed(TrimmedSurface::new(
                            bs.clone(),
                            umin,
                            umax,
                            v_min,
                            v_max,
                        ));
                    }
                }
            }
        }

        // L520-534.
        if s.is_u_periodic() {
            elclib_adjust_periodic(umin, umax, epsilon, &mut u_min, &mut u_max);
        } else if u_min > u_max {
            // L524-534.
            t = u_min;
            u_min = u_max;
            u_max = t;
            if (umin - u_min > epsilon) || (u_max - umax > epsilon) {
                self.my_error = BRepLibFaceError::ParametersOutOfRange;
                return;
            }
        }

        // L536-550.
        if s.is_v_periodic() {
            elclib_adjust_periodic(vmin, vmax, epsilon, &mut v_min, &mut v_max);
        } else if v_min > v_max {
            // L540-550.
            t = v_min;
            v_min = v_max;
            v_max = t;
            if (vmin - v_min > epsilon) || (vmax - vmax > epsilon) {
                self.my_error = BRepLibFaceError::ParametersOutOfRange;
                return;
            }
        }

        // L552-556: compute infinite flags.
        let umininf = is_negative_infinite_value(u_min);
        let umaxinf = is_positive_infinite_value(u_max);
        let vmininf = is_negative_infinite_value(v_min);
        let vmaxinf = is_positive_infinite_value(v_max);

        // L558-563: closed flag.
        let uclosed = s.is_u_closed()
            && (u_min - umin).abs() < epsilon
            && (u_max - umax).abs() < epsilon;

        let vclosed = s.is_v_closed()
            && (v_min - vmin).abs() < epsilon
            && (v_max - vmax).abs() < epsilon;

        // L565-571: compute 3d curves and degenerate flag.
        let max_tol = tol_degen;
        let mut c_umin: Option<Curve3> = None;
        let mut c_umax: Option<Curve3> = None;
        let mut c_vmin: Option<Curve3> = None;
        let mut c_vmax: Option<Curve3> = None;
        let d_umin;
        let d_umax;
        let d_vmin;
        let d_vmax;
        let mut umin_tol = CONFUSION;
        let mut umax_tol = CONFUSION;
        let mut vmin_tol = CONFUSION;
        let mut vmax_tol = CONFUSION;

        if !umininf {
            // L575-576.
            c_umin = Some(surface_u_iso(&s, u_min));
            d_umin = Self::is_degenerated(c_umin.as_ref(), max_tol, &mut umin_tol);
        } else {
            d_umin = false;
        }
        if !umaxinf {
            // L580-581.
            c_umax = Some(surface_u_iso(&s, u_max));
            d_umax = Self::is_degenerated(c_umax.as_ref(), max_tol, &mut umax_tol);
        } else {
            d_umax = false;
        }
        if !vmininf {
            // L585-586.
            c_vmin = Some(surface_v_iso(&s, v_min));
            d_vmin = Self::is_degenerated(c_vmin.as_ref(), max_tol, &mut vmin_tol);
        } else {
            d_vmin = false;
        }
        if !vmaxinf {
            // L590-591.
            c_vmax = Some(surface_v_iso(&s, v_max));
            d_vmax = Self::is_degenerated(c_vmax.as_ref(), max_tol, &mut vmax_tol);
        } else {
            d_vmax = false;
        }

        // L594-620: compute vertices — B.MakeVertex(V, S->Value(...), tol).
        let mut v00 = Shape::null();
        let mut v10 = Shape::null();
        let mut v11 = Shape::null();
        let mut v01 = Shape::null();

        if !umininf {
            if !vmininf {
                // L603.
                let p = s.point_at(u_min, v_min);
                v00 = brep.add_tvertex(p);
                bb.update_vertex_point(brep, v00.clone(), p, umin_tol.max(vmin_tol));
            }
            if !vmaxinf {
                // L607.
                let p = s.point_at(u_min, v_max);
                v01 = brep.add_tvertex(p);
                bb.update_vertex_point(brep, v01.clone(), p, umin_tol.max(vmax_tol));
            }
        }
        if !umaxinf {
            if !vmininf {
                // L614.
                let p = s.point_at(u_max, v_min);
                v10 = brep.add_tvertex(p);
                bb.update_vertex_point(brep, v10.clone(), p, umax_tol.max(vmin_tol));
            }
            if !vmaxinf {
                // L618.
                let p = s.point_at(u_max, v_max);
                v11 = brep.add_tvertex(p);
                bb.update_vertex_point(brep, v11.clone(), p, umax_tol.max(vmax_tol));
            }
        }

        // L622-632.
        if uclosed {
            v10 = v00.clone();
            v11 = v01.clone();
        }
        if vclosed {
            v01 = v00.clone();
            v11 = v10.clone();
        }

        // L634-649.
        if d_umin {
            v00 = v01.clone();
        }
        if d_umax {
            v10 = v11.clone();
        }
        if d_vmin {
            v00 = v10.clone();
        }
        if d_vmax {
            v01 = v11.clone();
        }

        // L651-668: make the lines — new Geom2d_Line(gp_Pnt2d(...),
        // gp_Dir2d(...)).
        let l_umin = if !umininf {
            Some(Curve2d::Line(Line2d {
                origin: DVec2::new(u_min, 0.0),
                direction: DVec2::new(0.0, 1.0),
            }))
        } else {
            None
        };
        let l_umax = if !umaxinf {
            Some(Curve2d::Line(Line2d {
                origin: DVec2::new(u_max, 0.0),
                direction: DVec2::new(0.0, 1.0),
            }))
        } else {
            None
        };
        let l_vmin = if !vmininf {
            Some(Curve2d::Line(Line2d {
                origin: DVec2::new(0.0, v_min),
                direction: DVec2::new(1.0, 0.0),
            }))
        } else {
            None
        };
        let l_vmax = if !vmaxinf {
            Some(Curve2d::Line(Line2d {
                origin: DVec2::new(0.0, v_max),
                direction: DVec2::new(1.0, 0.0),
            }))
        } else {
            None
        };

        // L670-672: make the face — B.MakeFace(F, S, Precision::Confusion()).
        self.my_shape = brep.add_tface_tol(
            Some(s.clone()),
            Shape::null(),
            Vec::new(),
            None,
            None,
            Vec::new(),
            false,
            CONFUSION,
        );
        let f = self.my_shape.clone();

        // L675: make the edges.
        let mut e_umin = Shape::null();
        let mut e_umax = Shape::null();
        let mut e_vmin = Shape::null();
        let mut e_vmax = Shape::null();

        if !umininf {
            // L679-686: B.MakeEdge(eumin, Cumin, uminTol) / B.MakeEdge(eumin).
            e_umin = brep.add_tedge(c_umin.clone(), Shape::null(), Shape::null(), [0.0, 0.0]);
            if uclosed {
                // L689: B.UpdateEdge(eumin, Lumax, Lumin, F, max(uminTol,
                // umaxTol)).
                bb.update_edge_pcurve(
                    brep,
                    e_umin.clone(),
                    l_umax.clone().expect("Lumax"),
                    f.clone(),
                    umin_tol.max(umax_tol),
                );
                bb.update_edge_pcurve(
                    brep,
                    e_umin.clone(),
                    l_umin.clone().expect("Lumin"),
                    f.clone(),
                    umin_tol.max(umax_tol),
                );
            } else {
                // L693: B.UpdateEdge(eumin, Lumin, F, uminTol).
                bb.update_edge_pcurve(
                    brep,
                    e_umin.clone(),
                    l_umin.clone().expect("Lumin"),
                    f.clone(),
                    umin_tol,
                );
            }
            // L695: B.Degenerated(eumin, Dumin).
            bb.set_edge_degenerated(brep, e_umin.clone(), d_umin);
            if !vmininf {
                // L698-699: V00.Orientation(TopAbs_FORWARD); B.Add(eumin, V00).
                v00.orientation = Orientation::Forward;
                builder_add_edge_vertex(brep, &e_umin, &v00);
            }
            if !vmaxinf {
                // L703-704.
                v01.orientation = Orientation::Reversed;
                builder_add_edge_vertex(brep, &e_umin, &v01);
            }
            // L706: B.Range(eumin, VMin, VMax).
            bb.set_edge_range(brep, e_umin.clone(), v_min, v_max);
        }

        if !umaxinf {
            if uclosed {
                // L713: eumax = eumin.
                e_umax = e_umin.clone();
            } else {
                // L717-726.
                e_umax = brep.add_tedge(c_umax.clone(), Shape::null(), Shape::null(), [0.0, 0.0]);
                bb.update_edge_pcurve(
                    brep,
                    e_umax.clone(),
                    l_umax.clone().expect("Lumax"),
                    f.clone(),
                    umax_tol,
                );
                bb.set_edge_degenerated(brep, e_umax.clone(), d_umax);
                if !vmininf {
                    // L729-730.
                    v10.orientation = Orientation::Forward;
                    builder_add_edge_vertex(brep, &e_umax, &v10);
                }
                if !vmaxinf {
                    // L734-735.
                    v11.orientation = Orientation::Reversed;
                    builder_add_edge_vertex(brep, &e_umax, &v11);
                }
                // L737: B.Range(eumax, VMin, VMax).
                bb.set_edge_range(brep, e_umax.clone(), v_min, v_max);
            }
        }

        if !vmininf {
            // L743-759.
            e_vmin = brep.add_tedge(c_vmin.clone(), Shape::null(), Shape::null(), [0.0, 0.0]);
            if vclosed {
                // L753: B.UpdateEdge(evmin, Lvmin, Lvmax, F, max(vminTol,
                // vmaxTol)).
                bb.update_edge_pcurve(
                    brep,
                    e_vmin.clone(),
                    l_vmin.clone().expect("Lvmin"),
                    f.clone(),
                    vmin_tol.max(vmax_tol),
                );
                bb.update_edge_pcurve(
                    brep,
                    e_vmin.clone(),
                    l_vmax.clone().expect("Lvmax"),
                    f.clone(),
                    vmin_tol.max(vmax_tol),
                );
            } else {
                // L757: B.UpdateEdge(evmin, Lvmin, F, vminTol).
                bb.update_edge_pcurve(
                    brep,
                    e_vmin.clone(),
                    l_vmin.clone().expect("Lvmin"),
                    f.clone(),
                    vmin_tol,
                );
            }
            // L759.
            bb.set_edge_degenerated(brep, e_vmin.clone(), d_vmin);
            if !umininf {
                // L762-763.
                v00.orientation = Orientation::Forward;
                builder_add_edge_vertex(brep, &e_vmin, &v00);
            }
            if !umaxinf {
                // L767-768.
                v10.orientation = Orientation::Reversed;
                builder_add_edge_vertex(brep, &e_vmin, &v10);
            }
            // L770: B.Range(evmin, UMin, UMax).
            bb.set_edge_range(brep, e_vmin.clone(), u_min, u_max);
        }

        if !vmaxinf {
            if vclosed {
                // L777: evmax = evmin.
                e_vmax = e_vmin.clone();
            } else {
                // L781-790.
                e_vmax = brep.add_tedge(c_vmax.clone(), Shape::null(), Shape::null(), [0.0, 0.0]);
                bb.update_edge_pcurve(
                    brep,
                    e_vmax.clone(),
                    l_vmax.clone().expect("Lvmax"),
                    f.clone(),
                    vmax_tol,
                );
                bb.set_edge_degenerated(brep, e_vmax.clone(), d_vmax);
                if !umininf {
                    // L793-794.
                    v01.orientation = Orientation::Forward;
                    builder_add_edge_vertex(brep, &e_vmax, &v01);
                }
                if !umaxinf {
                    // L798-799.
                    v11.orientation = Orientation::Reversed;
                    builder_add_edge_vertex(brep, &e_vmax, &v11);
                }
                // L801: B.Range(evmax, UMin, UMax).
                bb.set_edge_range(brep, e_vmax.clone(), u_min, u_max);
            }
        }

        // L806-807: eumin.Orientation(TopAbs_REVERSED);
        // evmax.Orientation(TopAbs_REVERSED).
        e_umin.orientation = Orientation::Reversed;
        e_vmax.orientation = Orientation::Reversed;

        if !umininf && !umaxinf && vmininf && vmaxinf {
            // L811-821: two wires in u.
            let w = bb.make_wire(brep);
            bb.add_to_wire(brep, w.clone(), e_umin.clone());
            add_face_wire(brep, bb, &f, &w);
            let w = bb.make_wire(brep);
            bb.add_to_wire(brep, w.clone(), e_umax.clone());
            add_face_wire(brep, bb, &f, &w);
            set_shape_closed(brep, &f, uclosed);
        } else if umininf && umaxinf && !vmininf && !vmaxinf {
            // L823-833: two wires in v.
            let w = bb.make_wire(brep);
            bb.add_to_wire(brep, w.clone(), e_vmin.clone());
            add_face_wire(brep, bb, &f, &w);
            let w = bb.make_wire(brep);
            bb.add_to_wire(brep, w.clone(), e_vmax.clone());
            add_face_wire(brep, bb, &f, &w);
            set_shape_closed(brep, &f, vclosed);
        } else if !umininf || !umaxinf || !vmininf || !vmaxinf {
            // L835-858: one wire.
            let w = bb.make_wire(brep);
            if !umininf {
                bb.add_to_wire(brep, w.clone(), e_umin.clone());
            }
            if !vmininf {
                bb.add_to_wire(brep, w.clone(), e_vmin.clone());
            }
            if !umaxinf {
                bb.add_to_wire(brep, w.clone(), e_umax.clone());
            }
            if !vmaxinf {
                bb.add_to_wire(brep, w.clone(), e_vmax.clone());
            }
            add_face_wire(brep, bb, &f, &w);
            set_shape_closed(brep, &w, !umininf && !umaxinf && !vmininf && !vmaxinf);
            set_shape_closed(brep, &f, uclosed && vclosed);
        }

        // L860-866.
        if offset_surface {
            // BRepLib::SameParameter(F, Precision::Confusion(), true) — GAP
            // leaf, see the module header.
            brep_lib_same_parameter_shape_gap(&f, CONFUSION, true);
        }

        // L868: Done().
        self.done();
    }

    /// OCCT BRepLib_MakeFace::Add(const TopoDS_Wire& W) (L873-879).
    pub fn add_wire(&mut self, brep: &mut BRep, bb: &mut BRepBuilder, w: &Shape) {
        // L875-876: BRep_Builder B; B.Add(myShape, W).
        add_face_wire(brep, bb, &self.my_shape.clone(), w);
        // L877: B.NaturalRestriction(TopoDS::Face(myShape), false).
        brep.face_mut(self.my_shape.clone()).natural_restriction = false;
        // L878: Done().
        self.done();
    }

    /// OCCT BRepLib_MakeFace::Face() (L883-886).
    pub fn face(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepLib_MakeFace::Error() (L897-900).
    pub fn error(&self) -> BRepLibFaceError {
        self.my_error
    }

    /// OCCT BRepLib_MakeFace::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BRepLib_MakeFace::CheckInside() (L907-924) — reverses the current
    /// face if not a bounded area.
    pub fn check_inside(&mut self, brep: &mut BRep) {
        // L910: TopoDS_Face F = TopoDS::Face(myShape).
        let f = self.my_shape.clone();
        // L911: BRepTopAdaptor_FClass2d FClass(F, 0.).
        let surf = match brep.face(f.clone()).surface.clone() {
            Some(s) => s,
            None => return,
        };
        let src = FaceShapeSource::new(&f, surf, &[glam::DAffine3::IDENTITY]);
        let fclass = FClass2d::new(&src, 0, 0.0);
        // L912-923.
        if fclass.perform_infinite_point(&src) == State::In {
            // L914-922: S = myShape.EmptyCopied(); the children are added
            // Reversed().
            let s = brep.empty_copied(&self.my_shape);
            let (outer, inners, internals) = {
                let fd = brep.face(self.my_shape.clone());
                (
                    fd.outer_wire.clone(),
                    fd.inner_wires.clone(),
                    fd.internal_vertices.clone(),
                )
            };
            let flip = |x: &Shape| -> Shape {
                let mut c = x.clone();
                c.orientation = match c.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    other => other,
                };
                c
            };
            {
                let fd = brep.face_mut(s.clone());
                fd.outer_wire = flip(&outer);
                fd.inner_wires = inners.iter().map(flip).collect();
                fd.internal_vertices = internals.iter().map(flip).collect();
            }
            // L923: myShape = S.
            self.my_shape = s;
        }
    }

    /// OCCT BRepLib_MakeShape::Done().
    fn done(&mut self) {
        self.done = true;
    }
}

// =========================================================================
// OCCT `BRep_Builder::Add(TopoDS_Face, TopoDS_Wire)` — the module-level
// helper shared by BRepLib_MakeFace::Add and the Init natural-bound assembly.
// =========================================================================

/// OCCT `BRep_Builder::Add(TopoDS_Face, TopoDS_Wire)` semantics: the FIRST
/// wire added to a face becomes its outer wire, later wires are inner wires;
/// both are appended to the face's sub-shape list (TopoDS_Iterator
/// enumerates them in insertion order — BRep_Builder.cxx / TopoDS_Builder).
/// The slot test uses the sub-shape TYPE: a face created without a wire
/// carries the null placeholder (a Vertex) in its outer-wire slot, whose
/// `Shape::index` is also usize::MAX for pool-external builders — judging by
/// `is_null()` would misclassify them.
fn add_face_wire(brep: &mut BRep, bb: &mut BRepBuilder, f: &Shape, w: &Shape) {
    let has_outer = matches!(&*brep.face(f.clone()).outer_wire.data, TShape::Wire(_));
    if has_outer {
        bb.add_to_face(brep, f.clone(), w.clone());
    } else {
        let fd = brep.face_mut(f.clone());
        fd.outer_wire = w.clone();
        fd.my_shapes.push(w.clone());
    }
}

// =========================================================================
// OCCT ElCLib::AdjustPeriodic (ElCLib.cxx L115-146) — local pure-math
// helper (the kernel copy is pub(crate) to the kernel crate).
// =========================================================================
fn elclib_adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    // OCCT L121.
    if is_infinite_value(u_first) || is_infinite_value(u_last) {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    let a_period = u_last - u_first;
    // OCCT L128-135: the Epsilon(ULast) overflow guard.
    if a_period < u_last.abs() * COMPUTATIONAL {
        *u1 = u_first;
        *u2 = u_last;
        return;
    }
    *u1 -= ((*u1 - u_first) / a_period).floor() * a_period;
    if u_last - *u1 < preci {
        *u1 -= a_period;
    }
    *u2 -= ((*u2 - *u1) / a_period).floor() * a_period;
    if *u2 - *u1 < preci {
        *u2 += a_period;
    }
}

/// OCCT Geom_Surface::UIso — the Geom_Plane arm (Geom_Plane.cxx L96-100):
/// `Geom_Line(ElSLib::PlaneUIso(pos, U))`, i.e. the line through the plane
/// point at (U, 0) along the V direction (ElSLib.cxx PlaneUIso).
/// `Surface3::Trimmed` mirrors Geom_RectangularTrimmedSurface::UIso, which
/// delegates to the basis iso; the trimmed span rides on the edge Range the
/// caller sets (BRepLib_MakeFace.cxx L706/L737).  The other surface types
/// keep the OCCT-failure GAP (the E3-T kernel gap, blend-surface UIso note).
fn surface_u_iso(s: &Surface3, u: f64) -> Curve3 {
    match s {
        Surface3::Plane(p) => Curve3::Line(Line3::new(p.origin + u * p.u_dir, p.v_dir)),
        Surface3::Trimmed(t) => surface_u_iso(t.basis.as_ref(), u),
        _ => unimplemented!("GAP: Geom_Surface::UIso pending (TKGeomBase/Geom)"),
    }
}

/// OCCT Geom_Surface::VIso — the Geom_Plane arm (Geom_Plane.cxx L103-107):
/// `Geom_Line(ElSLib::PlaneVIso(pos, V))`, the line through (0, V) along the
/// U direction.  See [`surface_u_iso`] for the Trimmed delegation and the
/// remaining GAP.
fn surface_v_iso(s: &Surface3, v: f64) -> Curve3 {
    match s {
        Surface3::Plane(p) => Curve3::Line(Line3::new(p.origin + v * p.v_dir, p.u_dir)),
        Surface3::Trimmed(t) => surface_v_iso(t.basis.as_ref(), v),
        _ => unimplemented!("GAP: Geom_Surface::VIso pending (TKGeomBase/Geom)"),
    }
}

/// OCCT BRepLib::UpdateTolerances(myShape) (BRepLib.cxx) — GAP leaf recorded
/// as a no-op: the shape-level tolerance update is pending; consumers of the
/// plane ctor do not depend on the adjustment.
fn brep_lib_update_tolerances_gap(_shape: &Shape) {}

/// OCCT BRepLib::SameParameter(myShape, tol, WithPCurve) (BRepLib.cxx
/// L301-456 is the edge-level body) — GAP leaf recorded as a no-op at the
/// shape-level entry; the edge-level home is topalgo/brep_lib/brep_lib.rs.
fn brep_lib_same_parameter_shape_gap(_shape: &Shape, _tol: f64, _with_pcurve: bool) {}
