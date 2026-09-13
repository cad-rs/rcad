// OCCT LocOpe_Generator.cxx L1234-1447 (static helpers) — 1:1 translation
// (part b of the LocOpe_Generator port; the class itself lives in
// loc_ope_generator.rs).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Generator.cxx
//
// Contents:
//   - static bool ToFuse(const TopoDS_Face&, const TopoDS_Face&)   L1234-1287
//   - static bool ToFuse(const TopoDS_Edge&, const TopoDS_Edge&)   L1291-1352
//   - static bool ToFuse(const TopoDS_Edge&, const TopoDS_Face&,
//                        const TopoDS_Vertex&, const map&)         L1356-1377
//   - static double NewParameter(Edg, Vtx, NewEdg, NewVtx)         L1381-1447
//
// Architecture differences (referenced from the affected functions):
// 1. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate).
// 2. BRep_Tool::Surface(F, loc) / Curve(E, loc, f, l): the feat pipeline
//    shapes carry identity locations (arch. diff. #1 of
//    loc_ope_find_edges.rs), so the OCCT "apply loc.Transformation()"
//    branches are no-ops here — the geometry travels untransformed.
// 3. Geom surface/curve DynamicType comparison maps to the rcad enum
//    discriminant comparison (std::mem::discriminant); Geom_Plane /
//    Geom_Line / Geom_Circle map to Surface3::Plane / Curve3::Line /
//    Curve3::Circle.
// 4. The gp statics (gp_Ax3::IsCoplanar, gp_Ax1::IsCoaxial,
//    gp_Dir::IsEqual / Angle / AngleWithRef, gp_XYZ::CrossCrossed) are
//    re-hosted below as pure-math helpers from gp_Ax3.hxx L603-620,
//    gp_Ax1.cxx L30-44, gp_Dir.hxx L191-195, gp_Dir.cxx L55-87,
//    gp_XYZ.hxx L551-564.
// 5. ElCLib::LineParameter / CircleParameter (ElCLib.cxx L1192-1224) and
//    the ElCLib-local normalizeAngle (ElCLib.cxx L56-72) are re-hosted
//    below.  gp::Resolution() is Precision::Computational
//    (rcad_kernel::precision::COMPUTATIONAL).
// 6. Epsilon(value) (Standard_Real.hxx L242-246) — one ULP — is re-hosted
//    below (standard_epsilon) via f64::next_after.
//
// first consumer: LocOpe_Generator::Perform (loc_ope_generator.rs).

use crate::feat::brep_feat_builder::explorer;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{ Circle3, Curve3, Plane, Surface3, TrimmedSurface };
use rcad_kernel::precision::{ ANGULAR, COMPUTATIONAL, CONFUSION };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, TShape };
use std::collections::HashSet;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
pub(crate) fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopAbs::Reverse (TopAbs.hxx) — FORWARD<->REVERSED, INTERNAL/
/// EXTERNAL unchanged.
pub(crate) fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        o => o,
    }
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
pub(crate) fn with_orientation(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT TopoDS_Shape::Reversed() — a copy with reversed orientation.
pub(crate) fn shape_reversed(s: &Shape) -> Shape {
    with_orientation(s, top_abs_reverse(s.orientation))
}

/// OCCT BRep_Tool::Pnt(vtx).
pub(crate) fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(shape) — the vertex/edge/face tolerance.
pub(crate) fn brep_tool_tolerance(s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
pub(crate) fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Range(edg, f, l).
pub(crate) fn brep_tool_range(edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Surface(fac) — the face surface (local coordinates,
/// architecture difference #2).
pub(crate) fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::Curve(edg, f, l) — the 3D curve and parameter range
/// (architecture difference #2).
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edg, face, f, l) — the pcurve of the edge
/// on the face with its range (same re-host as loc_ope_gluer.rs).
pub(crate) fn brep_tool_curve_on_surface(edg: &Shape, face: &Shape) -> Option<(rcad_kernel::geom::Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT BRep_Tool::Parameter(vtx, edg) — the vertex parameter stored on the
/// edge (OCCT asserts the representation exists; the missing entry maps to
/// 0.0 with the assert surface noted).
pub(crate) fn brep_tool_parameter(vtx: &Shape, edg: &Shape) -> f64 {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.vertex_params.get(&vtx.ptr_id()).copied().unwrap_or(0.0),
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::IsClosed(edg, fac) — true when the edge carries two
/// pcurves on the face (the seam of a periodic surface;
/// BRep_Tool.cxx IsClosed(E, F)).
pub(crate) fn brep_tool_is_closed_edge_on_face(edg: &Shape, face: &Shape) -> bool {
    let Some((_, pcurve1)) = edge_pcurves_on_face(edg, face) else {
        return false;
    };
    pcurve1.is_some()
}

/// The (first, second) pcurve of the edge on the face; the second is present
/// only for the BRep_CurveOnClosedSurface representation (the seam).
fn edge_pcurves_on_face(
    edg: &Shape,
    face: &Shape,
) -> Option<(Option<rcad_kernel::geom::Curve2d>, Option<rcad_kernel::geom::Curve2d>)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let key = shape_key(face);
            let first = ed.pcurves.get(&key).map(|(c, _, _)| c.clone());
            // The closed-surface representation holds both pcurves; it also
            // feeds the keyed map (with the first pcurve), so the second is
            // detected through the representations list.
            let mut second = None;
            for rep in &ed.representations {
                if let rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    face: fkey,
                    pcurve1,
                    ..
                } = rep
                {
                    if *fkey == key {
                        second = Some(pcurve1.clone());
                    }
                }
            }
            if first.is_none() && second.is_none() {
                None
            } else {
                Some((first, second))
            }
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::UVPoints(E, F, PFirst, PLast) — the pcurve values at the
/// edge range bounds (the reversed + closed-surface second-pcurve branch of
/// BRep_Tool.cxx L1005-1089 included).
pub(crate) fn brep_tool_uv_points(edg: &Shape, face: &Shape) -> (DVec2, DVec2) {
    use rcad_kernel::geom::Curve2dEval;
    let (f, l) = brep_tool_range(edg);
    if let Some((first, second)) = edge_pcurves_on_face(edg, face) {
        let eisreversed = edg.orientation == Orientation::Reversed;
        let pick = if eisreversed { second.or(first) } else { first.or(second) };
        if let Some(c) = pick {
            return (c.point_at(f), c.point_at(l));
        }
    }
    (DVec2::ZERO, DVec2::ZERO)
}

/// OCCT Geom_Surface::DynamicType() == Geom_RectangularTrimmedSurface ->
/// BasisSurface() — strip the trimmed wrapper.
pub(crate) fn basis_surface(s: &Surface3) -> Surface3 {
    match s {
        Surface3::Trimmed(t) => {
            let TrimmedSurface { basis, .. } = t;
            (**basis).clone()
        }
        other => other.clone(),
    }
}

/// OCCT Geom_Curve::DynamicType() == Geom_TrimmedCurve -> BasisCurve() —
/// strip the trimmed wrapper.
pub(crate) fn basis_curve(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Trimmed(tc) => (*tc.curve).clone(),
        other => other.clone(),
    }
}

/// OCCT gp_Dir::IsEqual(theOther, theAngularTolerance) (gp_Dir.hxx L191-195
/// IsParallel sibling): the angle within the tolerance (pure-math helper,
/// architecture difference #4).
fn dir_is_equal(d1: DVec3, d2: DVec3, the_angular_tolerance: f64) -> bool {
    let dot = d1.dot(d2).clamp(-1.0, 1.0);
    let cross_modulus = d1.cross(d2).length();
    let an_ang = cross_modulus.atan2(dot); // gp_Dir::Angle: [0, PI]
    an_ang <= the_angular_tolerance
}

/// OCCT gp_Ax3::IsCoplanar(theOther, theLinearTolerance,
/// theAngularTolerance) (gp_Ax3.hxx L603-620) over the rcad Plane carrier
/// (architecture difference #4).
fn ax3_is_coplanar(
    pos1: &Plane,
    pos2: &Plane,
    the_linear_tolerance: f64,
    the_angular_tolerance: f64,
) -> bool {
    let a_vec = pos2.origin - pos1.origin;
    let d1 = pos1.normal.dot(a_vec).abs();
    let d2 = pos2.normal.dot(a_vec).abs();
    d1 <= the_linear_tolerance
        && d2 <= the_linear_tolerance
        && dir_is_parallel(pos1.normal, pos2.normal, the_angular_tolerance)
}

/// OCCT gp_Dir::IsParallel(theOther, theAngularTolerance) (gp_Dir.hxx
/// L191-195).
fn dir_is_parallel(d1: DVec3, d2: DVec3, the_angular_tolerance: f64) -> bool {
    let dot = d1.dot(d2).clamp(-1.0, 1.0);
    let cross_modulus = d1.cross(d2).length();
    let an_ang = cross_modulus.atan2(dot); // gp_Dir::Angle: [0, PI]
    an_ang <= the_angular_tolerance || std::f64::consts::PI - an_ang <= the_angular_tolerance
}

/// OCCT gp_Ax1::IsCoaxial(Other, AngularTolerance, LinearTolerance)
/// (gp_Ax1.cxx L30-44).
fn ax1_is_coaxial(
    loc1: DVec3,
    dir1: DVec3,
    loc2: DVec3,
    dir2: DVec3,
    angular_tolerance: f64,
    linear_tolerance: f64,
) -> bool {
    let xyz1 = (loc1 - loc2).cross(dir2).length(); // D1
    let xyz2 = (loc2 - loc1).cross(dir1).length(); // D2
    dir_is_equal(dir1, dir2, angular_tolerance)
        && xyz1 <= linear_tolerance
        && xyz2 <= linear_tolerance
}

/// OCCT gp_Dir::AngleWithRef(theOther, theVRef) (gp_Dir.cxx L55-87) — the
/// signed angle of <me> -> theOther measured about theVRef, in [-PI, PI].
fn dir_angle_with_ref(d: DVec3, other: DVec3, vref: DVec3) -> f64 {
    let xyz = d.cross(other);
    let cosinus = d.dot(other).clamp(-1.0, 1.0);
    let sinus = xyz.length();
    let ang = if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else if cosinus < 0.0 {
        std::f64::consts::PI - sinus.asin()
    } else {
        sinus.asin()
    };
    if xyz.dot(vref) >= 0.0 {
        ang
    } else {
        -ang
    }
}

/// OCCT ElCLib-local normalizeAngle (ElCLib.cxx L56-72) — normalize to
/// [0, 2*PI] preserving the exact 2*PI seam value.
fn normalize_angle(the_angle: &mut f64) {
    use std::f64::consts::PI;
    let pipi = PI + PI;
    let negative_resolution = -COMPUTATIONAL;
    while *the_angle < negative_resolution {
        *the_angle += pipi;
    }
    while *the_angle > pipi * (1.0 + COMPUTATIONAL) {
        *the_angle -= pipi;
    }
    if *the_angle < 0.0 {
        *the_angle = 0.0;
    }
}

/// OCCT Epsilon(theValue) (Standard_Real.hxx L242-246) — one ULP of the
/// value (architecture difference #6); f64::next_after is spelled through
/// the bit representation (the std intrinsic is not available on this
/// toolchain).
pub(crate) fn standard_epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        let bits = the_value.to_bits();
        let next = if bits == f64::INFINITY.to_bits() {
            the_value
        } else {
            f64::from_bits(bits + 1)
        };
        next - the_value
    } else {
        the_value - f64::from_bits(the_value.to_bits() - 1)
    }
}

/// OCCT GeomProjLib::Curve2d(C, f, l, S, tol) re-host over
/// rcad_kernel::base::geom_proj_lib::curve2d (the surface natural domain is
/// supplied the way the OCCT 5-argument overload obtains it internally).
pub(crate) fn geomproj_lib_curve2d(
    c: &Curve3,
    f: f64,
    l: f64,
    s: &Surface3,
    _tol: f64,
) -> Option<rcad_kernel::geom::Curve2d> {
    use rcad_kernel::geom::SurfaceEval;
    let dom = s.default_domain();
    rcad_kernel::base::geom_proj_lib::curve2d(c, f, l, s, dom[0], dom[1], dom[2], dom[3])
}

/// OCCT BRepTools::UVBounds(F, Umin, Umax, Vmin, Vmax) (BRepTools.cxx
/// L64-80) re-host — the bounds of the pcurve endpoints over the face edges.
pub(crate) fn brep_tools_uv_bounds(fac: &Shape) -> [f64; 4] {
    use rcad_kernel::geom::Curve2dEval;
    let mut umin = f64::MAX;
    let mut umax = f64::MIN;
    let mut vmin = f64::MAX;
    let mut vmax = f64::MIN;
    for edg in explorer(fac, ShapeType::Edge, ShapeType::Shape) {
        let Some((c2d, f, l)) = brep_tool_curve_on_surface(&edg, fac) else {
            continue;
        };
        let p1 = c2d.point_at(f);
        let p2 = c2d.point_at(l);
        for p in [p1, p2] {
            umin = umin.min(p.x);
            umax = umax.max(p.x);
            vmin = vmin.min(p.y);
            vmax = vmax.max(p.y);
        }
    }
    [umin, umax, vmin, vmax]
}

/// OCCT BRepTools::IsReallyClosed(E, F) (BRepTools.cxx L1204-1220) re-host —
/// the edge is closed on the face and appears exactly twice in the
/// exploration.
pub(crate) fn brep_tools_is_really_closed(e: &Shape, f: &Shape) -> bool {
    if !brep_tool_is_closed_edge_on_face(e, f) {
        return false;
    }
    let mut nbocc = 0;
    for cur in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        if cur.is_same(e) {
            nbocc += 1;
        }
    }
    nbocc == 2
}

/// OCCT static ToFuse(const TopoDS_Face& F1, const TopoDS_Face& F2)
/// (cxx L1234-1287).
pub(crate) fn tofuse_face_face(f1: &Shape, f2: &Shape) -> bool {
    if f1.is_null() || f2.is_null() {
        return false;
    }

    // OCCT cxx L1244-1245.
    let tollin = CONFUSION;
    let tolang = ANGULAR;

    // OCCT cxx L1247-1268: S1 = Surface(F1, loc1); S2 = Surface(F2, loc2);
    // trimmed -> basis (architecture difference #2: identity locations).
    let s1 = brep_tool_surface(f1).map(|s| basis_surface(&s));
    let s2 = brep_tool_surface(f2).map(|s| basis_surface(&s));
    let (Some(s1), Some(s2)) = (s1, s2) else {
        // OCCT would dereference null handles — the rcad surfaces are
        // Option-carrying; the DynamicType comparison below fails for None.
        return false;
    };

    let typ_s1 = std::mem::discriminant(&s1);
    let typ_s2 = std::mem::discriminant(&s2);

    if typ_s1 != typ_s2 {
        return false;
    }

    let mut val_ret = false;
    if let Surface3::Plane(pl1_data) = &s1 {
        // OCCT cxx L1271-1283: pl1/pl2 transformed by the (identity)
        // locations; IsCoplanar.
        let pl1 = *pl1_data;
        let Surface3::Plane(pl2_data) = &s2 else {
            return false;
        };
        let pl2 = *pl2_data;
        if ax3_is_coplanar(&pl1, &pl2, tollin, tolang) {
            val_ret = true;
        }
    }

    val_ret
}

/// OCCT static ToFuse(const TopoDS_Edge& E1, const TopoDS_Edge& E2)
/// (cxx L1291-1352).
pub(crate) fn tofuse_edge_edge(e1: &Shape, e2: &Shape) -> bool {
    if e1.is_null() || e2.is_null() {
        return false;
    }

    // OCCT cxx L1302-1304.
    let tollin = CONFUSION;
    let tolang = ANGULAR;

    // OCCT cxx L1306-1321: C1/C2 with the (identity) location applied;
    // trimmed -> basis (architecture difference #2).
    let Some((c1_raw, _, _)) = brep_tool_curve(e1) else {
        return false;
    };
    let Some((c2_raw, _, _)) = brep_tool_curve(e2) else {
        return false;
    };
    let c1 = basis_curve(&c1_raw);
    let c2 = basis_curve(&c2_raw);

    let typ_c1 = std::mem::discriminant(&c1);
    let typ_c2 = std::mem::discriminant(&c2);

    if typ_c1 != typ_c2 {
        return false;
    }

    let mut val_ret = false;
    if let Curve3::Line(li1_data) = &c1 {
        // OCCT cxx L1340-1349: li1.Position().IsCoaxial(li2.Position(),
        // tolang, tollin).
        let Curve3::Line(li2_data) = &c2 else {
            return false;
        };
        if ax1_is_coaxial(li1_data.origin, li1_data.direction, li2_data.origin, li2_data.direction, tolang, tollin)
        {
            val_ret = true;
        }
    }

    val_ret
}

/// OCCT static ToFuse(const TopoDS_Edge& E, const TopoDS_Face& F,
/// const TopoDS_Vertex& V,
/// const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& toRemove)
/// (cxx L1356-1377).
pub(crate) fn tofuse_edge_face_vertex(
    e: &Shape,
    f: &Shape,
    v: &Shape,
    to_remove: &HashSet<(u64, u32)>,
) -> bool {
    // OCCT cxx L1361-1376.
    for eee in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        if !eee.is_same(e) {
            let (vf, vl) = top_exp_vertices(&eee);
            if (vf.is_same(v) || vl.is_same(v)) && !to_remove.contains(&shape_key(&eee)) {
                return false;
            }
        }
    }
    true
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) — the vertices of the edge
/// (TopExp.cxx: Vfirst = FirstVertex FORWARD / Vlast = LastVertex).
pub(crate) fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    match e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => (Shape::null(), Shape::null()),
    }
}

/// OCCT static NewParameter(Edg, Vtx, NewEdg, NewVtx) (cxx L1381-1447).
pub(crate) fn new_parameter(
    edg: &Shape,
    _vtx: &Shape,
    new_edg: &Shape,
    new_vtx: &Shape,
) -> f64 {
    // OCCT cxx L1392: gp_Pnt P = BRep_Tool::Pnt(NewVtx).
    let p = brep_tool_pnt(new_vtx);

    // OCCT cxx L1394-1405: C = Curve(Edg, loc, f, l); location applied
    // (identity, architecture difference #2); trimmed -> basis.
    let Some((c_raw, _, _)) = brep_tool_curve(edg) else {
        return 0.0;
    };
    let c = basis_curve(&c_raw);

    // OCCT cxx L1407-1410: Geom_Line -> ElCLib::Parameter(lin, P).
    if let Curve3::Line(lin) = &c {
        return (p - lin.origin).dot(lin.direction);
    }
    // OCCT cxx L1411-1444: Geom_Circle -> ElCLib::Parameter(circ, P) and
    // the 2*PI seam adjustment.
    if let Curve3::Circle(cir) = &c {
        let Circle3 {
            center,
            normal,
            x_dir,
            y_dir,
            radius: _,
        } = cir;
        let mut prm = elclib_circle_parameter(*center, *normal, *x_dir, *y_dir, p);
        // "Vtx vient d'une exploration de Edg orientee FORWARD".
        let orient = top_abs_reverse(_vtx.orientation);
        if orient == Orientation::Forward || orient == Orientation::Reversed {
            // OCCT cxx L1420: exp(NewEdg.Oriented(FORWARD), VERTEX).
            let mut prmmax: Option<f64> = None;
            for cur in explorer(
                &with_orientation(new_edg, Orientation::Forward),
                ShapeType::Vertex,
                ShapeType::Shape,
            ) {
                if cur.orientation == orient {
                    prmmax = Some(brep_tool_parameter(&cur, new_edg));
                    break;
                }
            }
            if let Some(prmmax) = prmmax {
                if (prmmax - prm).abs() <= standard_epsilon(2.0 * std::f64::consts::PI) {
                    if orient == Orientation::Reversed {
                        prm -= 2.0 * std::f64::consts::PI;
                    } else {
                        prm += 2.0 * std::f64::consts::PI;
                    }
                }
            }
        }
        return prm;
    }
    // OCCT cxx L1446: return 0.
    0.0
}

/// OCCT ElCLib::CircleParameter(Pos, P) (ElCLib.cxx L1199-1224) — pure-math
/// re-host (architecture difference #5).  gp::Resolution() maps to
/// Precision::Computational.
fn elclib_circle_parameter(center: DVec3, normal: DVec3, x_dir: DVec3, y_dir: DVec3, p: DVec3) -> f64 {
    let _ = y_dir;
    let a_vec = p - center;
    if a_vec.length_squared() < COMPUTATIONAL {
        // coinciding points -> infinite number of parameters
        return 0.0;
    }
    // Project vector on circle's plane: VProj = dir.CrossCrossed(aVec, dir)
    // (gp_XYZ.hxx L551-564).
    let cross = normal.cross(a_vec);
    let a_vproj = normal.cross(cross);
    if a_vproj.length_squared() < COMPUTATIONAL {
        return 0.0;
    }
    // Angle between X direction and projected vector.
    let mut teta = dir_angle_with_ref(x_dir, a_vproj, normal);
    normalize_angle(&mut teta);
    teta
}
