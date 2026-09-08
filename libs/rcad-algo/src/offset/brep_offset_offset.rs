// OCCT BRepOffset_Offset.cxx L1-1815 + BRepOffset_Offset.hxx L45-158 — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_Offset.cxx / .hxx / BRepOffset_Offset.lxx
//
// OCCT inheritance chain (hxx L45): none — BRepOffset_Offset is a standalone
// class.
//
// Architecture differences (continuing the numbering of
// brep_offset_make_simple_offset.rs):
// 8. NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> (Created / myMap /
//    MapSS) and NCollection_Map (VonDegen) -> HashMap/HashSet keyed by
//    (TShape ptr, Location); NCollection_List<TopoDS_Shape> (LEdge) ->
//    Vec<Shape>; NCollection_IndexedMap (GetFarestCorner) -> IndexMap.
// 9. BRepOffset::Surface (BRepOffset.cxx L44-365) — the offset-surface
//    factory; GAP: staged (the static is the next translation unit of this
//    package; the Init(Face) body is translated 1:1 against the GAP panic).
// 10. BRepOffset_Tool::Gabarit (BRepOffset_Tool.cxx L313-...) — GAP.
// 11. GeomFill_Pipe (TKGeomAlgo/GeomFill) — GAP carrier (both the
//     (Path, E1, E2, Offset) and (Curve, Offset) constructors).
// 12. GeomAPI_ExtremaCurveCurve / GeomAPI_ProjectPointOnCurve (TKGeomBase/
//     TKTopAlgo GeomAPI) — GAP carriers.
// 13. GeomLib::ExtendSurfByLength (TKTopAlgo/GeomLib) — GAP.
// 14. GeomProjLib::Curve2d (TKTopAlgo/GeomProjLib) — GAP (the
//     loc_ope_wires_on_shape_b.rs #10 precedent).
// 15. GeomConvert_ApproxSurface (TKGeomBase/GeomConvert) — GAP.
// 16. ShapeFix_Shape (TKShHealing/ShapeFix) — GAP carrier.
// 17. BRepGProp::LinearProperties + GProp_GProps::CentreOfMass
//     (TKTopAlgo/BRepGProp) — GAP carrier (the kernel base::gprop linear
//     re-host carries only the length).
// 18. GeomAdaptor_Surface/Geom2dAdaptor_Curve GetType() -> the rcad Surface3
//     / Curve2d variant matches (the same stand-in as the fillet
//     chfi3d_builder_0.rs surface_type_of).
// 19. The OCCT global TShape arena -> the class-local BRep pool (my_brep;
//     arch. diff. #4 of the MakeSimpleOffset module); `BRep_Builder`
//     locals -> BRepBuilder over that pool.
// 20. gp_Circ/gp_Lin iso constructions (ElSLib::*VIso/*UIso + Rotate/
//     Translate) — GAP leaves inside ComputeCurve3d (the TKMath/ElSLib iso
//     constructors are not translated; the leaf panics are the §0.6
//     annotation).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Surface3};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape::brep_tool_pnt;

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Status (BRepOffset_Status.hxx L27-33).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Status (BRepOffset_Status.hxx L27-33).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepOffsetStatus {
    Good,
    Reversed,
    Degenerated,
    Unknown,
}

// ---------------------------------------------------------------------------
// GAP carriers / leaves (architecture differences #9-#17, #20).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset::Surface (BRepOffset.cxx L44-202) — the elementary
/// offset-surface factory; the real body lives in [`super::brep_offset_surface`]
/// (the C.3 carrier switch; the old 3-argument stub is deleted).  The OCCT
/// hxx default is `allowC0 = false` (BRepOffset.hxx L44-47).
pub use super::brep_offset_surface::brep_offset_surface;

/// OCCT BRepOffset_Tool::Gabarit(aCurve) (BRepOffset_Tool.cxx L313-...) —
/// GAP (architecture difference #10).
pub(super) fn brep_offset_tool_gabarit(the_curve: &Curve3) -> f64 {
    let _ = the_curve;
    panic!("GAP: BRepOffset_Tool::Gabarit (BRepOffset_Tool.cxx not translated)");
}

/// OCCT GeomFill_Pipe (TKGeomAlgo/GeomFill) — GAP carrier (architecture
/// difference #11).
pub struct GeomFillPipe;

impl GeomFillPipe {
    /// OCCT GeomFill_Pipe::GeomFill_Pipe(Path, Edge1, Edge2, Offset).
    pub fn new(
        _the_path: &Curve3,
        _the_edge1: &Adaptor3dCurve,
        _the_edge2: &Adaptor3dCurve,
        _the_offset: f64,
    ) -> Self {
        panic!("GAP: GeomFill_Pipe (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Pipe::GeomFill_Pipe(Curve, Offset).
    pub fn new_from_curve(_the_curve: &Curve3, _the_offset: f64) -> Self {
        panic!("GAP: GeomFill_Pipe (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Pipe::Perform(Tol, Polynomial, Conti) /
    /// Perform() — the rcad carrier merges the two forms (the no-arg
    /// Perform passes the defaults).
    pub fn perform(&mut self, the_tol: f64, the_polynomial: bool, the_conti: GeomAbsShapeKind) {
        let _ = (the_tol, the_polynomial, the_conti);
        panic!("GAP: GeomFill_Pipe::Perform (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Pipe::IsDone().
    pub fn is_done(&self) -> bool {
        panic!("GAP: GeomFill_Pipe::IsDone (unreachable while the engine is a GAP)");
    }

    /// OCCT GeomFill_Pipe::ErrorOnSurf().
    pub fn error_on_surf(&self) -> f64 {
        panic!("GAP: GeomFill_Pipe::ErrorOnSurf (unreachable while the engine is a GAP)");
    }

    /// OCCT GeomFill_Pipe::Surface().
    pub fn surface(&self) -> Surface3 {
        panic!("GAP: GeomFill_Pipe::Surface (unreachable while the engine is a GAP)");
    }

    /// OCCT GeomFill_Pipe::ExchangeUV().
    pub fn exchange_uv(&self) -> bool {
        panic!("GAP: GeomFill_Pipe::ExchangeUV (unreachable while the engine is a GAP)");
    }
}

/// OCCT GeomAbs_Shape (TKG3d/GeomAbs/GeomAbs_Shape.hxx) — the continuity
/// kind carried by the Init(Path...) arguments (GeomAbs_C1 default).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsShapeKind {
    C0,
    G1,
    C1,
    C2,
    C3,
    CN,
}

/// OCCT Adaptor3d_Curve handle — either a 3d curve (GeomAdaptor_Curve) or a
/// curve on surface (Adaptor3d_CurveOnSurface); the Init(Path...) carrier
/// form (architecture difference #18).
pub enum Adaptor3dCurve {
    /// OCCT GeomAdaptor_Curve(C).
    GeomAdaptor(Curve3),
    /// OCCT Adaptor3d_CurveOnSurface(HC1, HS1).
    CurveOnSurface(Curve2d, Surface3),
}

/// OCCT GeomAPI_ExtremaCurveCurve (TKTopAlgo/GeomAPI) — GAP carrier
/// (architecture difference #12).
pub struct GeomApiExtremaCurveCurve;

impl GeomApiExtremaCurveCurve {
    /// OCCT GeomAPI_ExtremaCurveCurve(C1, C2).
    pub fn new(_the_c1: &Curve3, _the_c2: &Curve3) -> Self {
        panic!("GAP: GeomAPI_ExtremaCurveCurve (TKTopAlgo/GeomAPI not translated)");
    }

    /// OCCT GeomAPI_ExtremaCurveCurve::NearestPoints(P1, P2).
    pub fn nearest_points(&self) -> (DVec3, DVec3) {
        panic!("GAP: GeomAPI_ExtremaCurveCurve::NearestPoints (not translated)");
    }
}

/// OCCT GeomAPI_ProjectPointOnCurve (TKTopAlgo/GeomAPI) — GAP carrier
/// (architecture difference #12).
pub struct GeomApiProjectPointOnCurve;

impl GeomApiProjectPointOnCurve {
    /// OCCT GeomAPI_ProjectPointOnCurve(P, Curve).
    pub fn new(_the_p: DVec3, _the_curve: &Curve3) -> Self {
        panic!("GAP: GeomAPI_ProjectPointOnCurve (TKTopAlgo/GeomAPI not translated)");
    }

    /// OCCT GeomAPI_ProjectPointOnCurve::LowerDistanceParameter().
    pub fn lower_distance_parameter(&self) -> f64 {
        panic!("GAP: GeomAPI_ProjectPointOnCurve::LowerDistanceParameter (not translated)");
    }
}

/// OCCT GeomLib::ExtendSurfByLength(ProvokingSurf, Length, Contour,
/// InU, After) (TKTopAlgo/GeomLib, GeomLib.cxx) — GAP (architecture
/// difference #13); the OCCT in/out `occ::handle<Geom_BoundedSurface>&`
/// form becomes the (surf, u1, u2, v1, v2) tuple return.
pub(super) fn geom_lib_extend_surf_by_length(
    the_surf: &Surface3,
    the_length: f64,
    the_contour: i32,
    in_u: bool,
    after: bool,
) -> Surface3 {
    let _ = (the_surf, the_length, the_contour, in_u, after);
    panic!("GAP: GeomLib::ExtendSurfByLength (TKTopAlgo/GeomLib not translated)");
}

/// OCCT GeomProjLib::Curve2d(C, S) (TKTopAlgo/GeomProjLib) — GAP
/// (architecture difference #14).
pub(super) fn geom_proj_lib_curve2d(_the_c: &Curve3, _the_s: &Surface3) -> Option<Curve2d> {
    panic!("GAP: GeomProjLib::Curve2d (TKTopAlgo/GeomProjLib not translated)");
}

/// OCCT GeomConvert_ApproxSurface (TKGeomBase/GeomConvert) — GAP
/// (architecture difference #15).
pub(super) fn geom_convert_approx_surface(
    _the_surf: &Surface3,
    _the_tol: f64,
    _the_conti: GeomAbsShapeKind,
) -> Option<Surface3> {
    panic!("GAP: GeomConvert_ApproxSurface (TKGeomBase/GeomConvert not translated)");
}

/// OCCT ShapeFix_Shape (TKShHealing/ShapeFix) — the Init(Vertex) wire fixer
/// GAP carrier (architecture difference #16).
pub struct ShapeFixShape;

impl ShapeFixShape {
    /// OCCT ShapeFix_Shape::ShapeFix_Shape(S).
    pub fn new(_the_s: &Shape) -> Self {
        panic!("GAP: ShapeFix_Shape (TKShHealing/ShapeFix not translated)");
    }

    /// OCCT ShapeFix_Shape::Perform().
    pub fn perform(&mut self) {
        panic!("GAP: ShapeFix_Shape::Perform (TKShHealing/ShapeFix not translated)");
    }

    /// OCCT ShapeFix_Shape::Shape().
    pub fn shape(&self) -> Shape {
        panic!("GAP: ShapeFix_Shape::Shape (unreachable while the fixer is a GAP)");
    }
}

/// OCCT GProp_GProps (TKMath/GProp) + BRepGProp::LinearProperties — GAP
/// carrier (architecture difference #17).
pub struct GPropGProps {
    pub centre_of_mass: DVec3, // OCCT: CentreOfMass
}

/// OCCT BRepGProp::LinearProperties(S, VProps) — GAP.
pub(super) fn brep_gprop_linear_properties(_the_s: &Shape) -> GPropGProps {
    panic!("GAP: BRepGProp::LinearProperties (TKTopAlgo/BRepGProp not translated)");
}

/// OCCT BRepLib::BuildCurve3d(E, Tol) — GAP (the MakeCurve3d walk is not
/// translated; the MakeSimpleOffset module carries the same GAP for
/// BuildCurves3d).
pub(super) fn brep_lib_build_curve3d(the_e: &Shape, the_tol: f64) {
    let _ = (the_e, the_tol);
    panic!("GAP: BRepLib::BuildCurve3d (TKTopAlgo/BRepTools not translated)");
}

/// OCCT BRepTools::Update(S) (BRepTools.cxx) — GAP no-op re-host (the
/// same-parameter/tolerance consolidation walk; the face storage it
/// refreshes is carried by the rcad TFaceData fields directly).
pub(super) fn brep_tools_update(_the_s: &Shape) {
}

/// OCCT BRep_Builder::UpdateFace(F, S, L, Tol) — GAP no-op re-host (the rcad
/// TFaceData surface is written at construction; there is no post-hoc
/// surface setter on the rcad BRepBuilder).
pub(super) fn update_face_surface_gap(_the_f: &Shape, _the_s: &Surface3, _the_loc: u32, _the_tol: f64) {
}

/// OCCT BRep_Builder::UpdateEdge(E, C3d, Tol) / UpdateEdge(E, C3d, L, Tol) —
/// GAP no-op re-host (the rcad 3d-curve representation is index-based).
pub(super) fn update_edge_curve3d_gap(_the_e: &Shape, _the_c: &Curve3, _the_loc: u32, _the_tol: f64) {
}

/// OCCT BRep_Builder::Range(E, F, f, l) — GAP no-op re-host (face-keyed
/// range storage; the rcad pcurve range travels with the pcurve insert).
pub(super) fn set_edge_range_on_face_gap(_the_e: &Shape, _the_f: &Shape, _the_first: f64, _the_last: f64) {
}

/// OCCT Geom2d_Curve::Mirror(gp_Ax2d) — GAP no-op re-host (the rcad Curve2d
/// carries no in-place mirror).
pub(super) fn curve2d_mirror_gap(_the_c: &mut Curve2d, _the_point: DVec2, _the_dir: DVec2) {
}

// ---------------------------------------------------------------------------
// OCCT statics (BRepOffset_Offset.cxx L86-385).
// ---------------------------------------------------------------------------

// OCCT BRepOffset_Offset.cxx L86-111 — GetFarestCorner.
pub(super) fn get_farest_corner(a_wire: &Shape) -> DVec3 {
    // OCCT L88-89: NCollection_IndexedMap<TopoDS_Shape> Vertices;
    // TopExp::MapShapes(aWire, TopAbs_VERTEX, Vertices).
    let vertices = explorer(a_wire, ShapeType::Vertex, ShapeType::Shape);

    let mut max_dist = 0.0;
    let mut the_point = DVec3::ZERO;
    for i in 0..vertices.len() {
        for j in 0..vertices.len() {
            let v1 = &vertices[i];
            let v2 = &vertices[j];
            let p1 = brep_tool_pnt(v1).unwrap_or(DVec3::ZERO);
            let p2 = brep_tool_pnt(v2).unwrap_or(DVec3::ZERO);
            let a_dist = p1.distance_squared(p2);
            if a_dist > max_dist {
                max_dist = a_dist;
                the_point = p1;
            }
        }
    }

    the_point
}

//=================================================================================================
// OCCT BRepOffset_Offset.cxx L115-131 — static UpdateEdge(E, C, L, Tol): the
// 3d-curve form ("Cut curves to avoid copies in the extensions").
//=================================================================================================
pub(super) fn update_edge_3d(the_e: &Shape, the_c: &Curve3, the_l: u32, the_tol: f64) {
    // OCCT L121-130: BRep_Builder B; down-cast to Geom_TrimmedCurve and bind
    // the basis curve.
    let mut b = BRepBuilder::new();
    match the_c {
        Curve3::Trimmed(bc) => {
            // GAP no-op re-host: BRep_Builder::UpdateEdge(E, C3d, L, Tol).
            update_edge_curve3d_gap(the_e, &bc.curve, the_l, the_tol);
        }
        _ => {
            // OCCT L129: B.UpdateEdge(E, C, L, Tol).
            update_edge_curve3d_gap(the_e, the_c, the_l, the_tol);
        }
    }
    let _ = &mut b;
}

//=================================================================================================
// OCCT BRepOffset_Offset.cxx L135-151 — static UpdateEdge(E, C, F, Tol): the
// pcurve form ("Cut curves to avoid copies in the extensions").
//=================================================================================================
pub(super) fn update_edge_2d(the_e: &Shape, the_c: &Curve2d, the_f: &Shape, the_tol: f64, the_brep: &mut BRep) {
    // OCCT L141-150: BRep_Builder B; down-cast to Geom2d_TrimmedCurve and
    // bind the basis curve.
    let mut b = BRepBuilder::new();
    match the_c {
        Curve2d::Trimmed(bc) => {
            // OCCT L145: B.UpdateEdge(E, BC->BasisCurve(), F, Tol).
            b.update_edge_pcurve(
                the_brep,
                the_e.clone(),
                (*bc.curve).clone(),
                the_f.clone(),
                the_tol,
            );
        }
        _ => {
            // OCCT L149: B.UpdateEdge(E, C, F, Tol).
            b.update_edge_pcurve(the_brep, the_e.clone(), the_c.clone(), the_f.clone(), the_tol);
        }
    }
    let _ = &mut b;
}

//=================================================================================================
// OCCT BRepOffset_Offset.cxx L155-183 — static UpdateEdge(E, C1, C2, F, Tol):
// the seam pcurve form ("Cut curves to avoid copies in the extensions").
//=================================================================================================
pub(super) fn update_edge_2d_seam(
    the_e: &Shape,
    the_c1: &Curve2d,
    the_c2: &Curve2d,
    the_f: &Shape,
    the_tol: f64,
    the_brep: &mut BRep,
) {
    // OCCT L161-182: BRep_Builder B; the trimmed bases NC1/NC2.
    let mut b = BRepBuilder::new();
    let nc1 = match the_c1 {
        Curve2d::Trimmed(bc1) => (*bc1.curve).clone(),
        _ => the_c1.clone(),
    };
    let nc2 = match the_c2 {
        Curve2d::Trimmed(bc2) => (*bc2.curve).clone(),
        _ => the_c2.clone(),
    };
    // OCCT L182: B.UpdateEdge(E, NC1, NC2, F, Tol) — the rcad closed-surface
    // (seam) pcurve carrier takes the current edge range as the shared
    // parameter interval.
    let [a_first, a_last] = match the_e.as_edge() {
        Some(ed) => ed.range,
        None => [0.0, 0.0],
    };
    b.update_edge_pcurve_closed(
        the_brep,
        the_e.clone(),
        nc1,
        nc2,
        the_f.clone(),
        a_first,
        a_last,
        the_tol,
    );
    let _ = &mut b;
}

//=======================================================================
// function : ComputeCurve3d
// purpose  : Particular case of Curve On Surface.
//=======================================================================
// OCCT BRepOffset_Offset.cxx L190-385.
pub(super) fn compute_curve3d(
    the_brep: &mut BRep,
    edge: &Shape,
    curve: &Curve2d,
    surf: &Surface3,
    loc: u32,
    tol: f64,
) {
    // try to find the particular case
    // if not found call BRepLib::BuildCurve3d

    let mut is_computed = false;

    // Search only isos on analytic surfaces.
    // OCCT L204-206: Geom2dAdaptor_Curve C(Curve); GeomAdaptor_Surface
    // S(Surf); GetType() — the rcad variant matches (arch. diff. #18).
    let mut the_builder = BRepBuilder::new();

    let c_ty_is_line = matches!(curve, Curve2d::Line(_));
    let s_ty = surf;
    if !matches!(s_ty, Surface3::Plane(_)) {
        // if plane buildcurve3d manage KPart
        if c_ty_is_line {
            let line = match curve {
                Curve2d::Line(l) => *l,
                _ => unreachable!(),
            };
            let d = line.direction;
            let angular = rcad_kernel::core::precision::ANGULAR;
            // OCCT L213: D.IsParallel(gp::DX2d(), Precision::Angular()) —
            // the cross magnitude test of the unit directions.
            let is_parallel_dx = (d.x * 0.0 - d.y * 1.0).abs() < angular; // |d x DX2d|
            let is_parallel_dy = (d.x * 1.0 - d.y * 0.0).abs() < angular; // |d x DY2d|
            let p = line.origin;
            if is_parallel_dx {
                // Iso V.
                if let Surface3::Sphere(sph) = s_ty {
                    // OCCT L215-238.
                    if ((p.y.abs() - std::f64::consts::FRAC_PI_2).abs())
                        < rcad_kernel::core::precision::PCONFUSION
                    {
                        the_builder.set_edge_degenerated(the_brep, edge.clone(), true);
                    } else {
                        let circle = elslib_sphere_v_iso(sph, p.y);
                        let circle = circle_rotate(circle, sph.axis, p.x);
                        let mut circle = circle;
                        // OCCT L231: if (D.IsOpposite(gp::DX2d(), ...)).
                        if d.x < 0.0 {
                            circle = circle_reversed(circle);
                        }
                        update_edge_3d(edge, &Curve3::Circle(circle), loc, tol);
                    }
                    is_computed = true;
                } else if let Surface3::Cylinder(cyl) = s_ty {
                    // OCCT L239-255.
                    let circle = elslib_cylinder_v_iso(cyl, cyl.radius, p.y);
                    let circle = circle_rotate(circle, cyl.axis, p.x);
                    let mut circle = circle;
                    if d.x < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_3d(edge, &Curve3::Circle(circle), loc, tol);
                    is_computed = true;
                } else if let Surface3::Cone(cone) = s_ty {
                    // OCCT L256-272.
                    let circle = elslib_cone_v_iso(cone, cone.radius, cone.half_angle_rad, p.y);
                    let circle = circle_rotate(circle, cone.axis, p.x);
                    let mut circle = circle;
                    if d.x < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_3d(edge, &Curve3::Circle(circle), loc, tol);
                    is_computed = true;
                } else if let Surface3::Torus(tore) = s_ty {
                    // OCCT L273-289.
                    let circle =
                        elslib_torus_v_iso(tore, tore.major_radius, tore.minor_radius, p.y);
                    let circle = circle_rotate(circle, tore.axis, p.x);
                    let mut circle = circle;
                    if d.x < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_3d(edge, &Curve3::Circle(circle), loc, tol);
                    is_computed = true;
                }
            } else if is_parallel_dy {
                // Iso U.
                if let Surface3::Sphere(sph) = s_ty {
                    // OCCT L293-318: calculate iso 0; set to sameparameter
                    // (rotation of circle - offset of Y); transformation en
                    // iso U (= P.X()).
                    let mut circle = elslib_sphere_u_iso(sph, sph.radius, 0.0);
                    circle = circle_rotate(circle, sph.axis, p.y);
                    circle = circle_rotate(circle, sph.axis, p.x);
                    let mut circle = circle;
                    if d.y < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_3d(edge, &Curve3::Circle(circle), loc, tol);
                    is_computed = true;
                } else if let Surface3::Cylinder(cyl) = s_ty {
                    // OCCT L319-334.
                    let line = elslib_cylinder_u_iso(cyl, cyl.radius, p.x);
                    let tr = line.direction * p.y;
                    let mut line = line3_translate(line, tr);
                    if d.y < 0.0 {
                        line.direction = -line.direction;
                    }
                    update_edge_3d(edge, &Curve3::Line(line), loc, tol);
                    is_computed = true;
                } else if let Surface3::Cone(cone) = s_ty {
                    // OCCT L335-350.
                    let line = elslib_cone_u_iso(cone, cone.radius, cone.half_angle_rad, p.x);
                    let tr = line.direction * p.y;
                    let mut line = line3_translate(line, tr);
                    if d.y < 0.0 {
                        line.direction = -line.direction;
                    }
                    update_edge_3d(edge, &Curve3::Line(line), loc, tol);
                    is_computed = true;
                } else if let Surface3::Torus(tore) = s_ty {
                    // OCCT L351-366.
                    let circle =
                        elslib_torus_u_iso(tore, tore.major_radius, tore.minor_radius, p.x);
                    let circle = circle_rotate(circle, tore.axis, p.y);
                    let mut circle = circle;
                    if d.y < 0.0 {
                        circle = circle_reversed(circle);
                    }
                    update_edge_3d(edge, &Curve3::Circle(circle), loc, tol);
                    is_computed = true;
                }
            }
        }
    } else {
        // Cas Plan
        // OCCT L372: C3d = GeomAPI::To3d(Curve, S.Plane()) — the rcad
        // GeomAPI vehicle is the analytic projection of the 2d line on the
        // plane (the geom_api_to_2d inverse).
        if let Surface3::Plane(pl) = s_ty {
            let p0 = curve.point_at(0.0);
            let p1 = curve.point_at(1.0);
            let a = DVec3::new(p0.x, p0.y, 0.0);
            let b = DVec3::new(p1.x, p1.y, 0.0);
            let origin = pl.origin + pl.u_dir * a.x + pl.v_dir * a.y;
            let dir = (pl.u_dir * (b.x - a.x) + pl.v_dir * (b.y - a.y)).normalize_or_zero();
            let c3d = Curve3::Line(rcad_kernel::geom::Line3 {
                origin,
                direction: dir,
            });
            update_edge_3d(edge, &c3d, loc, tol);
        }
        is_computed = true;
    }
    let _ = is_computed;
    // OCCT L376-384: if (!IsComputed) { /* BRepLib::BuildCurves3d(Edge,Tol);
    // ... */ } — the comment-only branch.
}

// --- ComputeCurve3d GAP leaves (architecture difference #20) ---

type GpCirc = rcad_kernel::geom::Circle3;
type GpLin = rcad_kernel::geom::Line3;

/// OCCT ElSLib::SphereVIso(Axis, Radius, V) + Ci.Rotate(AxeRev, U) — GAP.
fn elslib_sphere_v_iso(_sph: &rcad_kernel::geom::SphericalSurface, _v: f64) -> GpCirc {
    panic!("GAP: ElSLib::SphereVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::CylinderVIso(Axis, Radius, V) — GAP.
fn elslib_cylinder_v_iso(
    _cyl: &rcad_kernel::geom::CylindricalSurface,
    _radius: f64,
    _v: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::CylinderVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::ConeVIso(Axis, RefRadius, SemiAngle, V) — GAP.
fn elslib_cone_v_iso(
    _cone: &rcad_kernel::geom::ConicalSurface,
    _ref_radius: f64,
    _semi_angle: f64,
    _v: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::ConeVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::TorusVIso(Axis, MajorR, MinorR, V) — GAP.
fn elslib_torus_v_iso(
    _tore: &rcad_kernel::geom::ToroidalSurface,
    _major: f64,
    _minor: f64,
    _v: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::TorusVIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::SphereUIso(Axis, Radius, U) — GAP.
fn elslib_sphere_u_iso(
    _sph: &rcad_kernel::geom::SphericalSurface,
    _radius: f64,
    _u: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::SphereUIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::CylinderUIso(Position, Radius, U) — GAP.
fn elslib_cylinder_u_iso(
    _cyl: &rcad_kernel::geom::CylindricalSurface,
    _radius: f64,
    _u: f64,
) -> GpLin {
    panic!("GAP: ElSLib::CylinderUIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::ConeUIso(Position, RefRadius, SemiAngle, U) — GAP.
fn elslib_cone_u_iso(
    _cone: &rcad_kernel::geom::ConicalSurface,
    _ref_radius: f64,
    _semi_angle: f64,
    _u: f64,
) -> GpLin {
    panic!("GAP: ElSLib::ConeUIso (TKMath/ElSLib not translated)");
}

/// OCCT ElSLib::TorusUIso(Axis, MajorR, MinorR, U) — GAP.
fn elslib_torus_u_iso(
    _tore: &rcad_kernel::geom::ToroidalSurface,
    _major: f64,
    _minor: f64,
    _u: f64,
) -> GpCirc {
    panic!("GAP: ElSLib::TorusUIso (TKMath/ElSLib not translated)");
}

/// OCCT gp_Circ::Rotate(gp_Ax1, Angle) — GAP.
fn circle_rotate(_c: GpCirc, _axis: DVec3, _angle: f64) -> GpCirc {
    panic!("GAP: gp_Circ::Rotate (gp_Trsf rotation not translated)");
}

/// OCCT Geom_Circle::Reverse() — the parameter reversal (the rcad Circle3
/// carries no sense flag; the x/y frame swap is the direction reversal).
fn circle_reversed(mut c: GpCirc) -> GpCirc {
    std::mem::swap(&mut c.x_dir, &mut c.y_dir);
    c
}

/// OCCT gp_Lin::Translate(gp_Vec) — GAP.
fn line3_translate(mut l: GpLin, _tr: DVec3) -> GpLin {
    panic!("GAP: gp_Lin::Translate (gp_Trsf translation not translated)");
}

/// OCCT ElSLib::Parameters(Cone, P, U, V) — GAP.
pub(super) fn elslib_cone_parameters(
    _cone: &rcad_kernel::geom::ConicalSurface,
    _p: DVec3,
) -> (f64, f64) {
    panic!("GAP: ElSLib::Parameters(Cone) (TKMath/ElSLib not translated)");
}
