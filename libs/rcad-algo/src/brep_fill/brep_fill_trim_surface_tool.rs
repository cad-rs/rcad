//! OCCT BRepFill_TrimSurfaceTool — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepFill/
//!         BRepFill_TrimSurfaceTool.hxx (L33-77) +
//!         BRepFill_TrimSurfaceTool.cxx (L62-487).
//!
//! First consumer: BRepFill_Evolved::ElementaryPerform (the evolved batch).
//!
//! Architecture differences (referenced from the affected functions):
//! 1. `handle<Geom2d_Curve>` maps to the kernel [`rcad_kernel::geom::Curve2d`]
//!    enum; `Geom2d_TrimmedCurve` maps to `Curve2d::Trimmed(TrimmedCurve2)`.
//! 2. `Geom2dAdaptor_Curve(C)` (natural-domain form) maps to the
//!    `brep_fill_trim_edge_tool::Geom2dAdaptorCurve` re-host constructed on
//!    the curve default domain.
//! 3. `Geom2dInt_GInter` maps to
//!    `geomalgo::geom2d_int::TheIntPCurvePCurveOfGInter` + the bounded
//!    `IntRes2d_Domain` pair (the trim_edge_tool.rs precedent).
//! 4. `GeomAdaptor_Surface` (the EvalParameters degenerate branch) maps to a
//!    direct match on the kernel `Surface3` enum (the quadric Position
//!    encodings live on the surface types).
//! 5. `BRepIntCurveSurface_Inter::Init(Face, Line, Tol)` maps to the rcad
//!    `topalgo::brep_int_curve_surface::Inter` load + init_curve pair (the
//!    rcad Inter takes the explicit curve/surface form).
//! 6. OCCT `EvalParameters` mutates the caller's `Bis` handle in the
//!    elongate branch (`TBis->SetTrim(...)` through the shared handle —
//!    cxx L222-225); rcad mirrors the aliasing with a `&mut Curve2d`
//!    parameter (documented side effect, OCCT source as written).
//! 7. GAP carriers (dependencies outside the module, plan §0.6):
//!    `BRepFill_MultiLine` / `BRepFill_ApproxSeewing` (TKBool/BRepFill —
//!    separate classes, own batch; consumed only by
//!    [`BRepFillTrimSurfaceTool::project`]); the ProjLib_ComputeApprox
//!    approximation branch of GeomProjLib::Curve2d (TKGeomBase/ProjLib).
//!    The OCCT_DEBUG block (cxx L87-93 / L446-453) is compiled out.

use glam::{DVec2, DVec3};

use rcad_kernel::base::proj_lib::{
    project_pln_circ, project_pln_elips, project_pln_hypr, project_pln_lin, project_pln_parab,
};
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, Line3, Plane, Surface3, SurfaceEval, TrimmedCurve2,
    TrimmedCurve3,
};
use rcad_kernel::math::gp::Ax3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ShapeType, TShape};

use crate::brep_algo::tool::{brep_tool_curve, brep_tool_pnt, brep_tool_uv_points, explorer};
use crate::brep_fill::brep_fill_trim_edge_tool::{
    Geom2dAPIProjectPointOnCurve, Geom2dAdaptorCurve,
};
use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Curve2dType, TheIntPCurvePCurveOfGInter};
use crate::geomalgo::int_res2d::Domain as Res2dDomain;
use crate::topalgo::brep_int_curve_surface::inter::Inter as BRepIntCurveSurfaceInter;

/// OCCT static BRepFill_Precision() — the 1.e-6 used by EvalParameters
/// (cxx L178) and IsOnFace (cxx L417/L424).
const TOL_1E6: f64 = 1.0e-6;

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp re-hosts (the brep_algo::tool reductions)
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Degenerated(E).
fn brep_tool_degenerated(e: &Shape) -> bool {
    matches!(e.data.as_ref(), TShape::Edge(ed) if ed.degenerated)
}

/// OCCT TopExp::FirstVertex(E) — the first stored vertex (the returned
/// orientation is irrelevant for the point read).
fn top_exp_first_vertex(e: &Shape) -> Shape {
    match e.data.as_ref() {
        TShape::Edge(ed) => {
            if !ed.first.is_null() {
                ed.first.clone()
            } else {
                ed.last.clone()
            }
        }
        _ => Shape::null(),
    }
}

/// OCCT BRep_Tool::Surface(F, L) — the face surface (the location is
/// carried by the surface encoding in rcad; architecture difference).
fn brep_tool_surface(f: &Shape) -> Option<Surface3> {
    match f.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT `C = BRep_Tool::Curve(E, L, f, l)` + the
/// `CT = new Geom_TrimmedCurve(C, f, l)` construction +
/// `CT->Transform(L.Transformation())` — the trimmed copy of the edge 3d
/// curve (the rcad curve encodings carry the location — the
/// offset_wire.rs architecture difference #5 form).
fn edge_trimmed_curve(e: &Shape) -> Option<Curve3> {
    let (c, f, l) = brep_tool_curve(e)?;
    Some(Curve3::Trimmed(TrimmedCurve3::new(c, f, l)))
}

/// OCCT GeomProjLib::Curve2d(C, Plane) — the analytic dispatch of the
/// projection of a 3d curve on a plane (ProjLib::Project(Pl, C) +
/// ProjLib::MakePCurveOfType).  A trimmed curve projects its basis and
/// keeps the parameter range (the GeomProjLib trimmed-curve form).
pub(super) fn geom_proj_lib_curve2d(c: &Curve3, plane: &Plane) -> Curve2d {
    match c {
        Curve3::Line(l) => Curve2d::Line(project_pln_lin(plane, l)),
        Curve3::Circle(ci) => Curve2d::Circle(project_pln_circ(plane, ci)),
        Curve3::Ellipse(e) => Curve2d::Ellipse(project_pln_elips(plane, e)),
        Curve3::Parabola(p) => Curve2d::Parabola(project_pln_parab(plane, p)),
        Curve3::Hyperbola(h) => Curve2d::Hyperbola(project_pln_hypr(plane, h)),
        Curve3::Trimmed(t) => {
            // OCCT: the projection of a trimmed curve is the projected
            // basis restricted to the same parameter range.
            let base = geom_proj_lib_curve2d(&t.curve, plane);
            Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(base),
                t_min: t.first,
                t_max: t.last,
            })
        }
        // OCCT: the non-analytic branch builds a ProjLib_ComputedCurve
        // through ProjLib_ComputeApprox (GeomProjLib.cxx Curve2d) —
        // GAP (plan §0.6).
        _ => panic!(
            "GAP: GeomProjLib::Curve2d non-analytic branch needs ProjLib_ComputeApprox \
             (TKGeomBase/ProjLib not translated) — see file header"
        ),
    }
}

/// OCCT GeomAdaptor_Surface::Value(U, V) — the surface evaluation (the
/// Surface3 dispatch of architecture difference #4).
fn surface_value(s: &Surface3, u: f64, v: f64) -> DVec3 {
    match s {
        Surface3::Plane(p) => SurfaceEval::point_at(p, u, v),
        Surface3::Cylinder(p) => SurfaceEval::point_at(p, u, v),
        Surface3::Sphere(p) => SurfaceEval::point_at(p, u, v),
        Surface3::Cone(p) => SurfaceEval::point_at(p, u, v),
        Surface3::Torus(p) => SurfaceEval::point_at(p, u, v),
        Surface3::Revolution(p) => SurfaceEval::point_at(p, u, v),
        _ => panic!("GeomAdaptor_Surface::Value: unsupported surface type"),
    }
}

// ---------------------------------------------------------------------------
// GAP carriers (architecture difference #7)
// ---------------------------------------------------------------------------

/// GAP: BRepFill_MultiLine (TKBool/BRepFill — own class, own batch;
/// consumed only by [`BRepFillTrimSurfaceTool::project`]).
pub struct BRepFillMultiLine;

impl BRepFillMultiLine {
    /// OCCT BRepFill_MultiLine::BRepFill_MultiLine(F1, F2, E1, E2, Inv1,
    /// Inv2, Bis) (BRepFill_MultiLine.cxx).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        _the_face1: &Shape,
        _the_face2: &Shape,
        _the_edge1: &Shape,
        _the_edge2: &Shape,
        _the_inv1: bool,
        _the_inv2: bool,
        _the_bis: &Curve2d,
    ) -> Self {
        panic!("GAP: BRepFill_MultiLine (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_MultiLine::Continuity().
    pub fn continuity(&self) -> i32 {
        panic!("GAP: BRepFill_MultiLine (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_MultiLine::IsParticularCase().
    pub fn is_particular_case(&self) -> bool {
        panic!("GAP: BRepFill_MultiLine (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_MultiLine::Curves(Curve, PCurve1, PCurve2).
    pub fn curves(&self) -> (Option<Curve3>, Option<Curve2d>, Option<Curve2d>) {
        panic!("GAP: BRepFill_MultiLine (TKBool/BRepFill not translated) — see file header")
    }
}

/// GAP: BRepFill_ApproxSeewing (TKBool/BRepFill — own class, own batch;
/// consumed only by [`BRepFillTrimSurfaceTool::project`]).
pub struct BRepFillApproxSeewing;

impl BRepFillApproxSeewing {
    /// OCCT BRepFill_ApproxSeewing::BRepFill_ApproxSeewing(ML)
    /// (BRepFill_ApproxSeewing.cxx L45-49).
    pub fn new(_the_ml: &BRepFillMultiLine) -> Self {
        panic!("GAP: BRepFill_ApproxSeewing (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_ApproxSeewing::Curve().
    pub fn curve(&self) -> Option<Curve3> {
        panic!("GAP: BRepFill_ApproxSeewing (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_ApproxSeewing::CurveOnF1().
    pub fn curve_on_f1(&self) -> Option<Curve2d> {
        panic!("GAP: BRepFill_ApproxSeewing (TKBool/BRepFill not translated) — see file header")
    }

    /// OCCT BRepFill_ApproxSeewing::CurveOnF2().
    pub fn curve_on_f2(&self) -> Option<Curve2d> {
        panic!("GAP: BRepFill_ApproxSeewing (TKBool/BRepFill not translated) — see file header")
    }
}

// ---------------------------------------------------------------------------
// Static helpers (BRepFill_TrimSurfaceTool.cxx)
// ---------------------------------------------------------------------------

/// OCCT static Bubble (L101-119): order the sequence of points by
/// increasing x.
fn bubble(seq: &mut Vec<DVec3>) {
    let mut invert = true;
    let nb_points = seq.len() as i32;
    while invert {
        invert = false;
        for i in 1..nb_points {
            // OCCT L110-112: P1 = Seq.Value(i); P2 = Seq.Value(i + 1).
            let p1 = seq[(i - 1) as usize];
            let p2 = seq[i as usize];
            if p2.x < p1.x {
                seq.swap((i - 1) as usize, i as usize);
                invert = true;
            }
        }
    }
}

/// OCCT static EvalPhase (L123-153): the revolution phase (0 or M_PI) —
/// the point of parameter (0, V) on the surface is compared to the axis
/// X direction.
fn eval_phase(edge: &Shape, face: &Shape, gas: &Surface3, axis: &Ax3) -> f64 {
    // OCCT L128-132: BRep_Tool::UVPoints(Edge, Face, PE1, PE2); VDeg = PE1.Y().
    let (pe1, _pe2) = brep_tool_uv_points(edge, face);
    let v_deg = pe1.y;
    // OCCT L133-142: the V of the first face edge that is not Edge.
    let mut v = 0.0;
    let exp = explorer(face, ShapeType::Edge, ShapeType::Shape);
    for current in &exp {
        if !current.is_same(edge) {
            let (pf1, pf2) = brep_tool_uv_points(current, face);
            // OCCT L139.
            v = if (pf1.y - v_deg).abs() > (pf2.y - v_deg).abs() {
                pf1.y
            } else {
                pf2.y
            };
            break;
        }
    }
    // OCCT L143: P = GAS.Value(0., V).
    let p = surface_value(gas, 0.0, v);

    // OCCT L145-152.
    if (axis.axis.location - p).dot(axis.x_direction) < 0.0 {
        std::f64::consts::PI
    } else {
        0.0
    }
}

/// OCCT static EvalParameters (L155-378): the intersection parameters of
/// the bisectrice with the edge, ordered along the bisectrice.
///
/// Architecture difference #6: the OCCT `const handle(Geom2d_Curve)& Bis`
/// is mutated in the elongate branch through the shared handle
/// (`TBis->SetTrim(...)`, cxx L222-225); the rcad mirror takes `&mut`.
fn eval_parameters(edge: &Shape, face: &Shape, bis: &mut Curve2d, seq: &mut Vec<DVec3>) {
    let degener = brep_tool_degenerated(edge);
    // OCCT L176-179: Tol = 1.e-6 (BRepFill_Precision()); TolC = 0.
    let tol = TOL_1E6;
    let tol_c = 0.0;

    if !degener {
        // OCCT L183-185: C = BRep_Tool::Curve(Edge, L, f, l); CT = new
        // Geom_TrimmedCurve(C, f, l); CT->Transform(L.Transformation()).
        let ct = edge_trimmed_curve(edge).expect("BRep_Tool::Curve");
        // OCCT L187-188: C2d = GeomProjLib::Curve2d(CT, Plane) — the
        // projection of the 3d curves in the plane xOy.
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let mut c2d = geom_proj_lib_curve2d(&ct, &plane);

        // OCCT L189-190: Geom2dAdaptor_Curve AC(C2d); ABis(Bis).
        let ac_dom = c2d.default_domain();
        let (ac_f, ac_l) = (ac_dom[0], ac_dom[1]);
        let mut ac = Geom2dAdaptorCurve::new(c2d.clone(), ac_f, ac_l);
        let ab_dom = bis.default_domain();
        let (ab_f, ab_l) = (ab_dom[0], ab_dom[1]);
        let mut ab_is = Geom2dAdaptorCurve::new(bis.clone(), ab_f, ab_l);

        // OCCT L192: Intersector = Geom2dInt_GInter(ABis, AC, TolC, Tol).
        let mut intersector = TheIntPCurvePCurveOfGInter::new();
        let domain_bis = bounded_domain(bis, tol_c);
        let domain_ac = bounded_domain(&c2d, tol_c);
        intersector.perform(&ab_is, &domain_bis, &ac, &domain_ac, tol_c, tol);

        if !intersector.base.is_done() {
            // OCCT L194-197: throw StdFail_NotDone.
            panic!("StdFail_NotDone: BRepFill_TrimSurfaceTool::IntersectWith");
        }

        let mut nb_points = intersector.base.nb_points();

        if nb_points < 1 {
            // OCCT L203-210: try to elongate curves and enlarge tolerance
            // (not done right away in order not to get extra solutions).
            let c_type = Curve2dAdaptor::get_type(&ac);
            let bis_type = Curve2dAdaptor::get_type(&ab_is);
            let can_elongate_c = c_type != Curve2dType::BezierCurve
                && c_type != Curve2dType::BSplineCurve
                && c_type != Curve2dType::OffsetCurve
                && c_type != Curve2dType::OtherCurve;
            let can_elongate_bis = bis_type != Curve2dType::BezierCurve
                && bis_type != Curve2dType::BSplineCurve
                && bis_type != Curve2dType::OffsetCurve
                && bis_type != Curve2dType::OtherCurve;

            // OCCT L214-215: TBis = down_cast<Geom2d_TrimmedCurve>(Bis);
            // TC2d = down_cast<Geom2d_TrimmedCurve>(C2d).
            if can_elongate_c {
                // OCCT L219-220: TC2d->SetTrim(First-Tol, Last+Tol);
                // AC.Load(TC2d).
                elongate(&mut c2d, tol);
                let dom2 = c2d.default_domain();
                let (f2, l2) = (dom2[0], dom2[1]);
                ac = Geom2dAdaptorCurve::new(c2d.clone(), f2, l2);
            }
            if can_elongate_bis {
                // OCCT L222-225: TBis->SetTrim(First-Tol, Last+Tol);
                // ABis.Load(TBis) — the mutation reaches the caller's Bis
                // through the shared handle (architecture difference #6).
                elongate(bis, tol);
                let dom2 = bis.default_domain();
                let (f2, l2) = (dom2[0], dom2[1]);
                ab_is = Geom2dAdaptorCurve::new(bis.clone(), f2, l2);
            }
            // OCCT L227: Intersector = Geom2dInt_GInter(ABis, AC, TolC, Tol * 10).
            let mut intersector = TheIntPCurvePCurveOfGInter::new();
            let domain_bis = bounded_domain(bis, tol_c);
            let domain_ac = bounded_domain(&c2d, tol_c);
            intersector.perform(&ab_is, &domain_bis, &ac, &domain_ac, tol_c, tol * 10.0);

            if !intersector.base.is_done() {
                // OCCT L229-232: throw StdFail_NotDone.
                panic!("StdFail_NotDone: BRepFill_TrimSurfaceTool::IntersectWith");
            }

            nb_points = intersector.base.nb_points();
        }

        // OCCT L237-247: the intersection points
        // (U1 = ParamOnFirst, U2 = ParamOnSecond).
        if nb_points > 0 {
            for i in 1..=nb_points {
                let u1 = intersector.base.point(i).param_on_first();
                let u2 = intersector.base.point(i).param_on_second();
                let p = DVec3::new(u1, u2, 0.0);
                seq.push(p);
            }
        }

        // OCCT L249-269: the intersection segments.
        let nb_segments = intersector.base.nb_segments();
        if nb_segments > 0 {
            for i in 1..=nb_segments {
                let seg = intersector.base.segment(i);
                let mut u1 = seg.first_point().param_on_first();
                u1 += seg.last_point().param_on_first();
                u1 /= 2.0;
                let mut u2 = seg.first_point().param_on_second();
                u2 += seg.last_point().param_on_second();
                u2 /= 2.0;
                let p = DVec3::new(u1, u2, 0.0);
                seq.push(p);
            }
        }
        // OCCT L270-271: order the sequence by increasing parameter on the
        // bissectrice.
        bubble(seq);

        // OCCT L273-287: remove double points.  The loop bound NbPoints is
        // the point count captured before the segment append (the OCCT
        // source as written).
        let mut i: i32 = 1;
        while i < nb_points as i32 {
            let p1 = seq[(i - 1) as usize];
            let p2 = seq[i as usize];
            if p2.x - p1.x < tol {
                seq.remove((i - 1) as usize);
                i -= 1;
                nb_points -= 1;
            }
            i += 1;
        }
    } else {
        // OCCT L291-295: the edge is degenerated: the point and it is
        // found if it is on the bissectrice.
        let p3d = brep_tool_pnt(&top_exp_first_vertex(edge)).expect("FirstVertex");
        let p2d = DVec2::new(p3d.x, p3d.y);

        let mut u_bis = bis.first_parameter();
        let mut p_bis = bis.point_at(u_bis);

        // OCCT L305-308: Precision::IsPositiveInfinite guard (inside
        // gp_Pnt2d::Distance: Infinite * Infinite => DefaultNumericError).
        if is_positive_infinite(p_bis.x)
            || is_positive_infinite(p_bis.y)
            || p_bis.distance(p2d) > tol
        {
            // OCCT L309-313.
            u_bis = bis.last_parameter();
            if u_bis >= f64::MAX / 4.0 {
                return;
            }
            p_bis = bis.point_at(u_bis);
            if p_bis.distance(p2d) > tol {
                return;
            }
        }

        // OCCT L322-325: GS = BRep_Tool::Surface(Face);
        // GeomAdaptor_Surface GAS(GS); gp_Ax3 Axis; Phase = 0 — the
        // surface Position (architecture difference #4).
        let gs = brep_tool_surface(face).expect("BRep_Tool::Surface");

        let axis: Ax3;
        let mut phase = 0.0;

        match &gs {
            Surface3::Sphere(_) => {
                // OCCT L331-332: Axis = GAS.Sphere().Position().
                axis = quadric_position(&gs);
            }
            Surface3::Cone(_) => {
                // OCCT L334-341: Axis = GAS.Cone().Position();
                // Phase = EvalPhase(Edge, Face, GAS, Axis).
                axis = quadric_position(&gs);
                phase = eval_phase(edge, face, &gs, &axis);
            }
            Surface3::Torus(_) => {
                // OCCT L343-344: Axis = GAS.Torus().Position().
                axis = quadric_position(&gs);
            }
            Surface3::Cylinder(_) => {
                // OCCT L346-347: Axis = GAS.Cylinder().Position().
                axis = quadric_position(&gs);
            }
            Surface3::Revolution(_) => {
                // OCCT L349-359: the surface of revolution axis + EvalPhase.
                axis = quadric_position(&gs);
                phase = eval_phase(edge, face, &gs, &axis);
            }
            _ => {
                // OCCT L361-362: throw Standard_NotImplemented.
                panic!("Standard_NotImplemented: BRepFill_TrimSurfaceTool");
            }
        }

        // OCCT L365-366: D12d = Bis->DN(UBis, 1); D1 = (D12d.X(), D12d.Y(), 0.).
        let d12d = bis.derivative_at(u_bis);
        let d1 = DVec3::new(d12d.x, d12d.y, 0.0);

        // OCCT L368: U = Axis.XDirection().AngleWithRef(D1, X ^ Y).
        let y_dir = axis.x_direction.cross(axis.axis.direction);
        let mut u = angle_with_ref(&axis.x_direction, &d1, &axis.x_direction.cross(y_dir));
        // OCCT L369-373.
        u += phase;
        if u < 0.0 {
            u += 2.0 * std::f64::consts::PI;
        }

        // OCCT L375-376.
        let p = DVec3::new(bis.first_parameter(), u, 0.0);
        seq.push(p);
    }
}

/// OCCT Precision::IsPositiveInfinite (the +Infinite test — Precision.hxx).
fn is_positive_infinite(v: f64) -> bool {
    v >= f64::MAX / 4.0
}

/// OCCT `Geom2d_TrimmedCurve::SetTrim(First-Tol, Last+Tol)` (cxx L219 /
/// L224).  The OCCT down_cast requires a trimmed curve (the projection of
/// a trimmed curve and the bisectrice are trimmed in the OCCT flow); a
/// null down-cast in OCCT would crash — mirrored as a panic.
fn elongate(c: &mut Curve2d, tol: f64) {
    match c {
        Curve2d::Trimmed(t) => {
            let f = t.t_min;
            let l = t.t_max;
            t.t_min = f - tol;
            t.t_max = l + tol;
        }
        _ => panic!("BRepFill_TrimSurfaceTool: SetTrim on a non-trimmed curve (OCCT null down-cast)"),
    }
}

/// The bounded IntRes2d_Domain of a 2d curve (the trim_edge_tool.rs
/// construction).
fn bounded_domain(c: &Curve2d, tol: f64) -> Res2dDomain {
    let dom = c.default_domain();
    let (f, l) = (dom[0], dom[1]);
    Res2dDomain::bounded(c.point_at(f), f, tol, c.point_at(l), l, tol)
}

/// OCCT gp_Vec::AngleWithRef(V, Ref) (pure-math helper; the draft.rs
/// re-host).
fn angle_with_ref(v1: &DVec3, v2: &DVec3, reference: &DVec3) -> f64 {
    let n1 = v1.length();
    let n2 = v2.length();
    if n1 <= 1.0e-15 || n2 <= 1.0e-15 {
        return 0.0;
    }
    let cross = v1.cross(*v2);
    let sign = if reference.dot(cross) < 0.0 { -1.0 } else { 1.0 };
    sign * cross.length().atan2(v1.dot(*v2))
}

/// OCCT GeomAdaptor_Surface quadric Position() — the (location, main
/// direction, X direction) frame of the quadric (architecture difference
/// #4: read off the kernel surface encodings).
fn quadric_position(s: &Surface3) -> Ax3 {
    match s {
        Surface3::Sphere(sp) => Ax3::from_pnt_n_vx(sp.center, sp.axis, sp.ref_dir),
        Surface3::Cone(co) => Ax3::from_pnt_n_vx(co.apex, co.axis, co.ref_dir),
        Surface3::Torus(to) => Ax3::from_pnt_n_vx(to.center, to.axis, to.ref_dir),
        Surface3::Cylinder(cy) => Ax3::from_pnt_n_vx(cy.origin, cy.axis, cy.ref_dir),
        Surface3::Revolution(rv) => Ax3::from_pnt_n_vx(rv.axis_origin, rv.axis_dir, DVec3::X),
        _ => panic!("GeomAdaptor_Surface::Position: unsupported surface type"),
    }
}

// ---------------------------------------------------------------------------
// BRepFill_TrimSurfaceTool (hxx L33-77)
// ---------------------------------------------------------------------------

/// OCCT BRepFill_TrimSurfaceTool (hxx L33-77) — the tool for the
/// construction of the regularities between two neighbor faces; Edge1 and
/// Edge2 are the parallel edges corresponding to the minimum iso on F1 and
/// F2 respectively; Inv1 and Inv2 show if Edge1 and Edge2 are returned
/// parallel.
pub struct BRepFillTrimSurfaceTool {
    my_face1: Shape, // OCCT: myFace1
    my_face2: Shape, // OCCT: myFace2
    my_edge1: Shape, // OCCT: myEdge1
    my_edge2: Shape, // OCCT: myEdge2
    my_inv1: bool,   // OCCT: myInv1
    my_inv2: bool,   // OCCT: myInv2
    my_bis: Curve2d, // OCCT: myBis
}

impl BRepFillTrimSurfaceTool {
    /// OCCT BRepFill_TrimSurfaceTool::BRepFill_TrimSurfaceTool(Bis, Face1,
    /// Face2, Edge1, Edge2, Inv1, Inv2) (cxx L72-94).  The OCCT_DEBUG
    /// block (L87-93) is compiled out.
    pub fn new(
        bis: &Curve2d,
        face1: &Shape,
        face2: &Shape,
        edge1: &Shape,
        edge2: &Shape,
        inv1: bool,
        inv2: bool,
    ) -> Self {
        BRepFillTrimSurfaceTool {
            my_face1: face1.clone(),
            my_face2: face2.clone(),
            my_edge1: edge1.clone(),
            my_edge2: edge2.clone(),
            my_inv1: inv1,
            my_inv2: inv2,
            my_bis: bis.clone(),
        }
    }

    /// OCCT BRepFill_TrimSurfaceTool::IntersectWith(EdgeOnF1, EdgeOnF2,
    /// Points) (cxx L382-404).
    pub fn intersect_with(&self, edge_on_f1: &Shape, edge_on_f2: &Shape, points: &mut Vec<DVec3>) {
        points.clear();
        let mut points2: Vec<DVec3> = Vec::new();

        // OCCT L389-390: EvalParameters(EdgeOnF1, myFace1, myBis, Points);
        // EvalParameters(EdgeOnF2, myFace2, myBis, Points2).  The elongate
        // branch mutates myBis through the shared handle in OCCT — the
        // `&mut` clone mirror of architecture difference #6.
        let mut bis = self.my_bis.clone();
        eval_parameters(edge_on_f1, &self.my_face1, &mut bis, points);
        eval_parameters(edge_on_f2, &self.my_face2, &mut bis, &mut points2);

        // OCCT L392-393: StdFail_NotDone_Raise_if(Lengths differ).
        if points.len() != points2.len() {
            panic!(
                "StdFail_NotDone: BRepFill_TrimSurfaceTool::IntersectWith: \
                 incoherent intersection"
            );
        }

        // OCCT L395-403: PSeq.SetZ(Points2.Value(i).Y()).
        let nb_points = points.len();
        for i in 1..=nb_points {
            let mut p_seq = points[i - 1];
            p_seq.z = points2[i - 1].y;
            points[i - 1] = p_seq;
        }
    }

    /// OCCT BRepFill_TrimSurfaceTool::IsOnFace(Point) (cxx L408-427).
    pub fn is_on_face(&self, point: DVec2) -> bool {
        // OCCT L410-411: gp_Pnt P(Point.X(), Point.Y(), 0.);
        // gp_Lin Line(P, gp::DZ()).
        let p = DVec3::new(point.x, point.y, 0.0);
        let line = Line3::new(p, DVec3::Z);

        let mut inter = BRepIntCurveSurfaceInter::new();

        // OCCT L417: Inter.Init(myFace1, Line, 1e-6); More() — the rcad
        // Inter takes the explicit curve/surface form (architecture
        // difference #5).
        if face_line_hit(&self.my_face1, &line, &mut inter) {
            return true;
        }

        // OCCT L424: Inter.Init(myFace2, Line, 1e-6); More().
        face_line_hit(&self.my_face2, &line, &mut inter)
    }

    /// OCCT BRepFill_TrimSurfaceTool::ProjOn(Point, Edge) (cxx L431-458).
    /// The OCCT_DEBUG block (L446-453) is compiled out.
    pub fn proj_on(&self, point: DVec2, edge: &Shape) -> f64 {
        // OCCT L433-437: C1 = BRep_Tool::Curve(Edge, L, f, l);
        // CT = new Geom_TrimmedCurve(C1, f, l); CT->Transform(...).
        let ct = edge_trimmed_curve(edge).expect("BRep_Tool::Curve");

        // OCCT L440-442: the projection of the 3d curve in the plane xOy.
        let plane = Plane::new(DVec3::ZERO, DVec3::Z);
        let c2d = geom_proj_lib_curve2d(&ct, &plane);

        // OCCT L445: Projector = Geom2dAPI_ProjectPointOnCurve(Point, C2d).
        let projector = Geom2dAPIProjectPointOnCurve::new_point_curve(point, &c2d);

        // OCCT L456-457: U = Projector.LowerDistanceParameter().
        projector.lower_distance_parameter()
    }

    /// OCCT BRepFill_TrimSurfaceTool::Project(U1, U2, Curve, PCurve1,
    /// PCurve2, Cont) (cxx L462-486) — returns
    /// (Curve, PCurve1, PCurve2, Cont).
    pub fn project(
        &self,
        u1: f64,
        u2: f64,
    ) -> (Option<Curve3>, Option<Curve2d>, Option<Curve2d>, i32) {
        // OCCT L469: CT = new Geom2d_TrimmedCurve(myBis, U1, U2).
        let ct = Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(self.my_bis.clone()),
            t_min: u1,
            t_max: u2,
        });
        // OCCT L470: ML = BRepFill_MultiLine(myFace1, myFace2, myEdge1,
        // myEdge2, myInv1, myInv2, CT).
        let ml = BRepFillMultiLine::new(
            &self.my_face1,
            &self.my_face2,
            &self.my_edge1,
            &self.my_edge2,
            self.my_inv1,
            self.my_inv2,
            &ct,
        );

        // OCCT L472: Cont = ML.Continuity().
        let cont = ml.continuity();

        if ml.is_particular_case() {
            // OCCT L474-476: ML.Curves(Curve, PCurve1, PCurve2).
            let (curve, pcurve1, pcurve2) = ml.curves();
            (curve, pcurve1, pcurve2, cont)
        } else {
            // OCCT L478-485: AppSeew = BRepFill_ApproxSeewing(ML);
            // Curve = AppSeew.Curve(); PCurve1 = AppSeew.CurveOnF1();
            // PCurve2 = AppSeew.CurveOnF2().
            let app_seew = BRepFillApproxSeewing::new(&ml);
            (
                app_seew.curve(),
                app_seew.curve_on_f1(),
                app_seew.curve_on_f2(),
                cont,
            )
        }
    }
}

/// The OCCT `Inter.Init(Face, Line, Tol)` + `Inter.More()` pair (cxx
/// L417/L424) on the rcad Inter (architecture difference #5): the face
/// surface + UV bounds stand in for the patch walk.
fn face_line_hit(face: &Shape, line: &Line3, inter: &mut BRepIntCurveSurfaceInter) -> bool {
    let Some(surf) = brep_tool_surface(face) else {
        return false;
    };
    let (u_min, u_max, v_min, v_max) = face_uv_bounds(&surf);
    let curve = Curve3::Line(line.clone());
    inter.load(face, TOL_1E6);
    inter.init_curve(&curve, &surf, u_min, u_max, v_min, v_max);
    inter.more()
}

/// The face UV bounds (the BRepTools::UVBounds reduction — the natural
/// domain of the surface encoding).
fn face_uv_bounds(surf: &Surface3) -> (f64, f64, f64, f64) {
    let [u_min, u_max, v_min, v_max] = SurfaceDomain::default_domain_of(surf);
    (u_min, u_max, v_min, v_max)
}

/// OCCT BRepTools::UVBounds reduction — the Surface3 dispatch of the
/// natural domain (architecture difference #4 form).
trait SurfaceDomain {
    fn default_domain_of(&self) -> [f64; 4];
}

impl SurfaceDomain for Surface3 {
    fn default_domain_of(&self) -> [f64; 4] {
        match self {
            Surface3::Plane(p) => SurfaceEval::default_domain(p),
            Surface3::Cylinder(p) => SurfaceEval::default_domain(p),
            Surface3::Sphere(p) => SurfaceEval::default_domain(p),
            Surface3::Cone(p) => SurfaceEval::default_domain(p),
            Surface3::Torus(p) => SurfaceEval::default_domain(p),
            Surface3::Revolution(p) => SurfaceEval::default_domain(p),
            Surface3::Ellipsoid(p) => SurfaceEval::default_domain(p),
            Surface3::BSpline(p) => SurfaceEval::default_domain(p),
            _ => [0.0, 0.0, 0.0, 0.0],
        }
    }
}
