// OCCT LocOpe_DPrism.hxx L38-85 + LocOpe_DPrism.cxx L58-753 — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_DPrism.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_DPrism.cxx
//
// OCCT inheritance chain (LocOpe_DPrism.hxx L38): none — LocOpe_DPrism is a
// standalone class in the OCCT 8.0 sources at $OCCT_SRC (no
// LocOpe_GeneratedShape virtuals to translate; see loc_ope_prism.rs).
// The IntPerf() member declared in the OCCT header (hxx L69) has no
// definition in LocOpe_DPrism.cxx — nothing to translate.
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//    TopTools_ShapeMapHasher> (myMap) — HashMap keyed by (TShape ptr,
//    Location); never iterated.
// 2. BRepFill_Evolved (TKBool/BRepFill) — the evolved-prism engine's single
//    1:1 body lives in crate::brep_fill::brep_fill_evolved (its correct
//    location per the module map); the local GAP carrier of this file was
//    retired and Perform maps to the OCCT Face overload
//    (LocOpe_DPrism.cxx L110: myDPrism.Perform(mySpine, myProfile,
//    gp::XOY()) with TopoDS_Face mySpine, hxx L73). Everything downstream
//    (the whole IsDone() body of both constructors) is translated 1:1
//    against it.
// 3. BRepLib_MakeVertex / BRepLib_MakeEdge(V1, V2) / BRepLib_MakeWire —
//    re-hosted below over a local rcad topods::BRep pool (the same vehicle
//    as feat::loc_ope_build_shape::perform): MakeEdge(V1, V2) is the
//    Geom_Line(gp_Lin(P1, P2-P1)) edge with range [0, |P2-P1|]
//    (BRepLib_MakeEdge.cxx L185-198).
// 4. BRepTools::UVBounds(Spine, ...) — the TFaceData uv_domain /
//    surface.default_domain() stand-in (the same vehicle as
//    feat::loc_ope_cs_intersector::face_uv_domain); the OCCT pcurve-walk
//    bounds computation is a GAP.
// 5. BRepLib::UpdateTolerances(myRes) — GAP panic (TKTopAlgo/BRepLib not
//    translated); unreachable while the BRepFill_Evolved Perform is a GAP.
// 6. TopoDS_Face::EmptyCopied() — the rcad BRep::empty_copied pool re-host
//    (same TShape data, no children, orientation/location carried).
// 7. TopExp_Explorer ExpS + ReInit() — the rcad explorer returns an owned
//    Vec; each OCCT pass (More/Next or ReInit) is a fresh iteration of the
//    same Vec.
// 8. gp_Pln Directness (Curves/BarycCurve): the rcad Plane always carries a
//    direct frame (v_dir = normal x u_dir), so the OCCT
//    "if (!P.Direct()) Normale.Reverse()" has no rcad false case.
// 9. The OCCT_DEBUG trace block of BarycCurve (cxx L722-731, with the
//    commented-out Normale.Reverse()) is dead source text — not translated.
//
// first consumer: BRepFeat_MakeDPrism (3b).

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
// OCCT BRepFill_Evolved (TKBool/BRepFill/BRepFill_Evolved.hxx / .cxx) — the
// single 1:1 body lives in crate::brep_fill::brep_fill_evolved; the local
// GAP carrier of this file was retired.
use crate::brep_fill::brep_fill_evolved::BRepFillEvolved;
use crate::brep_fill::brep_fill_trim_edge_tool::GeomAbsJoinType;
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{ Curve3, Line3, Surface3, TrimmedCurve3, CurveEval, SurfaceEval };
use rcad_kernel::math::gp::{ Ax1, Ax3, GP_RESOLUTION };
use rcad_kernel::topo::topods::{ BRep, BRepBuilder };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, TShape };
use std::collections::HashMap;

/// OCCT #define NECHANT 7 // voir BRepFeat.cxx (LocOpe_DPrism.cxx L58).
const NECHANT: i32 = 7;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
fn with_orientation(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT TopoDS_Shape::Reversed() — a copy with reversed orientation.
fn shape_reversed(s: &Shape) -> Shape {
    let mut c = s.clone();
    c.orientation = match c.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        o => o,
    };
    c
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Surface(face) (architecture difference: identity
/// location on standalone feat shapes).
fn brep_tool_surface(face: &Shape) -> Surface3 {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone().expect("Standard_NoSuchObject"),
        _ => panic!("Standard_NoSuchObject"),
    }
}

/// OCCT BRep_Tool::Tolerance(face).
fn brep_tool_tolerance(face: &Shape) -> f64 {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRepLib_MakeVertex(P) (BRepLib_MakeVertex.cxx) — a vertex carrying
/// P at the Precision::Confusion() default tolerance, over the local pool
/// (architecture difference #3).
fn brep_lib_make_vertex(pool: &mut BRep, b: &mut BRepBuilder, p: DVec3) -> Shape {
    b.add_vertex(pool, p, rcad_kernel::precision::CONFUSION)
}

/// OCCT BRepLib_MakeEdge(V1, V2) (BRepLib_MakeEdge.cxx L185-198) — the
/// Geom_Line(gp_Lin(P1, gp_Vec(P1, P2))) edge with range [0, |P2-P1|] and
/// the vertices attached (architecture difference #3).
fn brep_lib_make_edge_v_v(pool: &mut BRep, b: &mut BRepBuilder, v1: &Shape, v2: &Shape) -> Shape {
    let p1 = brep_tool_pnt(v1);
    let p2 = brep_tool_pnt(v2);
    let l = p1.distance(p2);
    if l <= GP_RESOLUTION {
        // OCCT cxx L190-194: myError = BRepLib_LineThroughIdenticPoints.
        panic!("BRepLib_LineThroughIdenticPoints");
    }
    let gl = Curve3::Line(Line3::new(p1, p2 - p1));
    // OCCT cxx L197: Init(GL, V1, V2, 0, l).
    b.add_edge(pool, Some(gl), v1.clone(), v2.clone(), [0.0, l])
}

/// OCCT BRepTools::UVBounds(S, Umin, Umax, Vmin, Vmax) (architecture
/// difference #4 — the uv_domain / default_domain stand-in).
fn brep_tools_uv_bounds(face: &Shape) -> [f64; 4] {
    match face.data.as_ref() {
        TShape::Face(fd) => match (&fd.uv_domain, &fd.surface) {
            (Some(d), _) => *d,
            (None, Some(s)) => s.default_domain(),
            (None, None) => [0.0, 0.0, 0.0, 0.0],
        },
        _ => [0.0, 0.0, 0.0, 0.0],
    }
}

/// OCCT BRepLib::UpdateTolerances(S) (architecture difference #5 — GAP).
fn brep_lib_update_tolerances(_s: &mut Shape) {
    panic!("GAP: BRepLib::UpdateTolerances (TKTopAlgo/BRepLib not translated)");
}

/// OCCT LocOpe_DPrism (LocOpe_DPrism.hxx L38-85).
pub struct LocOpeDPrism {
    my_d_prism: BRepFillEvolved, // OCCT: myDPrism (BRepFill_Evolved)
    my_res: Shape,               // OCCT: myRes
    my_spine: Shape,             // OCCT: mySpine (TopoDS_Face)
    my_profile: Shape,           // OCCT: myProfile (TopoDS_Wire)
    my_profile1: Shape,          // OCCT: myProfile1 (TopoDS_Edge)
    my_profile2: Shape,          // OCCT: myProfile2 (TopoDS_Edge)
    my_profile3: Shape,          // OCCT: myProfile3 (TopoDS_Edge)
    my_height: f64,              // OCCT: myHeight
    my_first_shape: Shape,       // OCCT: myFirstShape
    my_last_shape: Shape,        // OCCT: myLastShape
    // OCCT: myCurvs (NCollection_Sequence) — the hxx field has no use in
    // LocOpe_DPrism.cxx; kept as the form mirror.
    #[allow(dead_code)]
    my_curvs: Vec<Curve3>,
    // OCCT: myMap (NCollection_DataMap) — arch. diff. #1
    my_map: HashMap<(u64, u32), Vec<Shape>>,
}

impl LocOpeDPrism {
    /// OCCT LocOpe_DPrism::LocOpe_DPrism(Spine, Height1, Height2, Angle)
    /// (cxx L62-361) — Rust has no overloading: the `_heights_angle`
    /// suffix.
    pub fn with_heights_angle(
        spine: &Shape,
        height1: f64,
        height2: f64,
        angle: f64,
    ) -> Self {
        // OCCT cxx L66: mySpine(Spine).
        let my_spine = spine.clone();

        // OCCT cxx L70-72.
        let my_height = height1 + height2;
        let y = height1 * angle.sin();
        let z = height1 * angle.cos();

        // OCCT cxx L74.
        let mut pool = BRep::new();
        let mut b = BRepBuilder::new();
        let vert2 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., y, z));

        // OCCT cxx L76-79.
        let y1 = -height2 * angle.sin();
        let z1 = -height2 * angle.cos();
        let vert1 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., y1, z1));

        // OCCT cxx L81.
        let my_profile2 = brep_lib_make_edge_v_v(&mut pool, &mut b, &vert1, &vert2);

        // OCCT cxx L83-90.
        let mut umax = 0.;
        let mut umin = 0.;
        let mut vmin = 0.;
        let mut vmax = 0.;

        let bounds = brep_tools_uv_bounds(spine);
        umin = bounds[0];
        umax = bounds[1];
        vmin = bounds[2];
        vmax = bounds[3];
        let mut deltay = (umax - umin).max(vmax - vmin) + y.abs();
        deltay *= 2.;

        // OCCT cxx L93-94.
        let vert3 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., y + deltay, z));
        let my_profile3 = brep_lib_make_edge_v_v(&mut pool, &mut b, &vert2, &vert3);

        // OCCT cxx L96-103.
        umax = 0.;
        umin = 0.;
        vmin = 0.;
        vmax = 0.;

        let bounds = brep_tools_uv_bounds(spine);
        umin = bounds[0];
        umax = bounds[1];
        vmin = bounds[2];
        vmax = bounds[3];
        let mut deltay1 = (umax - umin).max(vmax - vmin) + y1.abs();
        deltay1 *= 2.;

        // OCCT cxx L105-106.
        let vert4 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., y1 + deltay1, z1));
        let my_profile1 = brep_lib_make_edge_v_v(&mut pool, &mut b, &vert4, &vert1);

        // OCCT cxx L108.
        let my_profile = b.build_wire(
            &mut pool,
            vec![my_profile1.clone(), my_profile2.clone(), my_profile3.clone()],
        );

        // OCCT cxx L110 (arch. diff. #2).
        let mut my_d_prism = BRepFillEvolved::new();
        my_d_prism.perform_with_face_spine(
            &my_spine,
            &my_profile,
            &Ax3::new(),
            GeomAbsJoinType::Arc,
            false,
        );

        let mut my_map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        let mut my_res = Shape::null();
        let mut my_first_shape = Shape::null();
        let mut my_last_shape = Shape::null();

        if my_d_prism.is_done() {
            // OCCT cxx L114-118.
            let mut bs = LocOpeBuildShape::new();
            let c = b.make_compound(&mut pool, Vec::new());
            let mut lfaces: Vec<Shape> = Vec::new();
            let mut lcomplete: Vec<Shape> = Vec::new();

            // OCCT cxx L121-123: ExpS(mySpine, TopAbs_EDGE); View.
            let exp_s: Vec<Shape> = explorer(&my_spine, ShapeType::Edge, ShapeType::Shape);
            let mut view: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
            // OCCT cxx L124-135.
            for es in &exp_s {
                let lffs = my_d_prism.generated_shapes(es, &my_profile1);
                for value in &lffs {
                    // OCCT cxx L130-133.
                    if view.insert(shape_key(value)) {
                        b.add_to_compound(&mut pool, c.clone(), value.clone());
                    }
                }
            }

            // OCCT cxx L137-141: theMapEF over C.
            let mut the_map_ef: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            map_shapes_and_ancestors(
                &c,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_map_ef,
            );
            view.clear();

            // OCCT cxx L144-177.
            for i in 1..=the_map_ef.len() {
                let (_, entry) = the_map_ef.get_index(i - 1).expect("theMapEF entry");
                // OCCT cxx L146.
                if entry.1.len() == 1 {
                    // OCCT cxx L148-149.
                    let edg = entry.0.clone();
                    let fac = entry.1.first().expect("First()").clone();
                    // OCCT cxx L150.
                    if view.insert(shape_key(&fac)) {
                        // OCCT cxx L152-153 (arch. diff. #6).
                        let new_face = pool.empty_copied(&fac);
                        // OCCT cxx L155-174.
                        let fac_forward = with_orientation(&fac, Orientation::Forward);
                        'wires: for wire in
                            explorer(&fac_forward, ShapeType::Wire, ShapeType::Shape)
                        {
                            // OCCT cxx L159: exp2 over the wire edges.
                            let mut hit = false;
                            for e2 in explorer(&wire, ShapeType::Edge, ShapeType::Shape) {
                                // OCCT cxx L162.
                                if e2.is_same(&edg) {
                                    // OCCT cxx L164-167: B.Add(newFace,
                                    // exp.Current()) — the first wire of the
                                    // empty-copied face (the rcad outer_wire
                                    // slot).
                                    let fd = pool.face_mut(new_face.clone());
                                    fd.outer_wire = wire.clone();
                                    lfaces.push(new_face.clone());
                                    lcomplete.push(new_face.clone());
                                    hit = true;
                                    break;
                                }
                            }
                            // OCCT cxx L170-173: if (exp2.More()) break.
                            if hit {
                                break 'wires;
                            }
                        }
                    }
                }
            }

            // OCCT cxx L179-180.
            bs.perform(&lfaces);
            my_first_shape = bs.shape().cloned().unwrap_or_else(Shape::null);

            // OCCT cxx L182: B.MakeCompound(D).
            let d = b.make_compound(&mut pool, Vec::new());

            // OCCT cxx L184-185: ExpS.ReInit(); View.Clear().
            view.clear();

            // OCCT cxx L187-198.
            for es in &exp_s {
                let lfls = my_d_prism.generated_shapes(es, &my_profile3);
                for value in &lfls {
                    // OCCT cxx L192-195.
                    if view.insert(shape_key(value)) {
                        b.add_to_compound(&mut pool, d.clone(), value.clone());
                    }
                }
            }

            // OCCT cxx L200-203.
            lfaces.clear();
            the_map_ef.clear();
            map_shapes_and_ancestors(
                &d,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_map_ef,
            );
            view.clear();

            // OCCT cxx L205-238.
            for i in 1..=the_map_ef.len() {
                let (_, entry) = the_map_ef.get_index(i - 1).expect("theMapEF entry");
                // OCCT cxx L207.
                if entry.1.len() == 1 {
                    // OCCT cxx L209-210.
                    let edg = entry.0.clone();
                    let fac = entry.1.first().expect("First()").clone();
                    // OCCT cxx L211.
                    if view.insert(shape_key(&fac)) {
                        // OCCT cxx L213-214 (arch. diff. #6).
                        let new_face = pool.empty_copied(&fac);
                        // OCCT cxx L216-235.
                        let fac_forward = with_orientation(&fac, Orientation::Forward);
                        'wires: for wire in
                            explorer(&fac_forward, ShapeType::Wire, ShapeType::Shape)
                        {
                            // OCCT cxx L220: exp2 over the wire edges.
                            let mut hit = false;
                            for e2 in explorer(&wire, ShapeType::Edge, ShapeType::Shape) {
                                // OCCT cxx L223.
                                if e2.is_same(&edg) {
                                    // OCCT cxx L225-228.
                                    let fd = pool.face_mut(new_face.clone());
                                    fd.outer_wire = wire.clone();
                                    lfaces.push(new_face.clone());
                                    lcomplete.push(new_face.clone());
                                    hit = true;
                                    break;
                                }
                            }
                            // OCCT cxx L231-234.
                            if hit {
                                break 'wires;
                            }
                        }
                    }
                }
            }
            // OCCT cxx L239-240.
            bs.perform(&lfaces);
            my_last_shape = bs.shape().cloned().unwrap_or_else(Shape::null);

            // OCCT cxx L242.
            view.clear();

            // OCCT cxx L244-355.
            for es in &exp_s {
                let lffs = my_d_prism.generated_shapes(es, &my_profile2);

                // OCCT cxx L249-255: the first EDGE of lffs (the OCCT list
                // iterator stops there).
                let mut removed_edge: Option<Shape> = None;
                for value in &lffs {
                    if value.shape_type() == ShapeType::Edge {
                        removed_edge = Some(value.clone());
                        break;
                    }
                }
                // OCCT cxx L256: if (it.More()).
                if let Some(removed_edge) = removed_edge {
                    // OCCT cxx L258-262.
                    let mut new_face = Shape::null();
                    let new_wire = b.make_wire(&mut pool);
                    let mut orref = Orientation::Forward;
                    // OCCT cxx L264-311.
                    'faces: for value in &lffs {
                        if value.shape_type() == ShapeType::Face {
                            // OCCT cxx L268-269: theWire = the FORWARD-
                            // oriented face's first wire.
                            let value_forward = with_orientation(value, Orientation::Forward);
                            let the_wire = explorer(
                                &value_forward,
                                ShapeType::Wire,
                                ShapeType::Shape,
                            )
                            .into_iter()
                            .next()
                            .expect("exp.Current()");
                            if new_face.is_null() {
                                // OCCT cxx L272: S = BRep_Tool::Surface(...).
                                let mut s = brep_tool_surface(value);
                                // OCCT cxx L273-276: RectangularTrimmedSurface
                                // -> BasisSurface.
                                if let Surface3::Trimmed(ts) = &s {
                                    s = (*ts.basis).clone();
                                }
                                // OCCT cxx L277-280.
                                if !matches!(s, Surface3::Plane(_)) {
                                    break 'faces;
                                }

                                // OCCT cxx L282 (arch. diff.: a wire-less
                                // face; the wire attaches at the B.Add site).
                                let made = pool.add_tface_tol(
                                    Some(s),
                                    Shape::null(),
                                    Vec::new(),
                                    None,
                                    None,
                                    Vec::new(),
                                    false,
                                    brep_tool_tolerance(value),
                                );
                                // OCCT cxx L283.
                                new_face = with_orientation(&made, Orientation::Forward);
                                // OCCT cxx L284.
                                orref = the_wire.orientation;
                                // OCCT cxx L285-291.
                                let wire_forward =
                                    with_orientation(&the_wire, Orientation::Forward);
                                for e in explorer(&wire_forward, ShapeType::Edge, ShapeType::Shape)
                                {
                                    if !e.is_same(&removed_edge) {
                                        b.add_to_wire(&mut pool, new_wire.clone(), e);
                                    }
                                }
                            } else {
                                // OCCT cxx L295-308.
                                let wire_forward =
                                    with_orientation(&the_wire, Orientation::Forward);
                                for e in explorer(&wire_forward, ShapeType::Edge, ShapeType::Shape)
                                {
                                    if !e.is_same(&removed_edge) {
                                        if the_wire.orientation != orref {
                                            // Les 2 faces planes ont des
                                            // normales opposees
                                            b.add_to_wire(&mut pool, new_wire.clone(), e);
                                        } else {
                                            b.add_to_wire(
                                                &mut pool,
                                                new_wire.clone(),
                                                shape_reversed(&e),
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // OCCT cxx L312-319.
                    if !new_face.is_null() {
                        let fd = pool.face_mut(new_face.clone());
                        fd.outer_wire = with_orientation(&new_wire, orref);
                        lcomplete.push(new_face.clone());
                        my_map.insert(shape_key(es), Vec::new());
                        my_map
                            .get_mut(&shape_key(es))
                            .expect("myMap(ES)")
                            .push(new_face);
                    } else {
                        // OCCT cxx L322-329.
                        for value in &lffs {
                            if view.insert(shape_key(value))
                                && value.shape_type() == ShapeType::Face
                            {
                                lcomplete.push(value.clone());
                            }
                        }
                    }
                } else {
                    // OCCT cxx L333-339.
                    for value in &lffs {
                        if view.insert(shape_key(value))
                            && value.shape_type() == ShapeType::Face
                        {
                            lcomplete.push(value.clone());
                        }
                    }
                }

                // OCCT cxx L342-354.
                for es2 in explorer(es, ShapeType::Vertex, ShapeType::Shape) {
                    let ls2 = my_d_prism.generated_shapes(&es2, &my_profile2);
                    for value in &ls2 {
                        if view.insert(shape_key(value))
                            && value.shape_type() == ShapeType::Face
                        {
                            lcomplete.push(value.clone());
                        }
                    }
                }
            }

            // OCCT cxx L357-359.
            bs.perform(&lcomplete);
            my_res = bs.shape().cloned().unwrap_or_else(Shape::null);
            brep_lib_update_tolerances(&mut my_res);
        }

        LocOpeDPrism {
            my_d_prism,
            my_res,
            my_spine,
            my_profile,
            my_profile1,
            my_profile2,
            my_profile3,
            my_height,
            my_first_shape,
            my_last_shape,
            my_curvs: Vec::new(),
            my_map,
        }
    }

    /// OCCT LocOpe_DPrism::LocOpe_DPrism(Spine, Height, Angle)
    /// (cxx L365-552).
    pub fn with_height_angle(spine: &Shape, height: f64, angle: f64) -> Self {
        // OCCT cxx L366-369.
        let my_spine = spine.clone();
        let my_height = height;
        let y = height * angle.sin();
        let z = height * angle.cos();

        // OCCT cxx L373-375.
        let mut pool = BRep::new();
        let mut b = BRepBuilder::new();
        let vert1 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., 0., 0.));
        let vert2 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., y, z));
        let my_profile2 = brep_lib_make_edge_v_v(&mut pool, &mut b, &vert1, &vert2);

        // OCCT cxx L377-380.
        let bounds = brep_tools_uv_bounds(spine);
        let umin = bounds[0];
        let umax = bounds[1];
        let vmin = bounds[2];
        let vmax = bounds[3];
        let mut deltay = (umax - umin).max(vmax - vmin) + y.abs();
        deltay *= 2.;

        // OCCT cxx L382-383.
        let vert3 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., y + deltay, z));
        let my_profile3 = brep_lib_make_edge_v_v(&mut pool, &mut b, &vert2, &vert3);

        // OCCT cxx L385-386.
        let vert4 = brep_lib_make_vertex(&mut pool, &mut b, DVec3::new(0., deltay, 0.));
        let my_profile1 = brep_lib_make_edge_v_v(&mut pool, &mut b, &vert4, &vert1);

        // OCCT cxx L388.
        let my_profile = b.build_wire(
            &mut pool,
            vec![my_profile1.clone(), my_profile2.clone(), my_profile3.clone()],
        );

        // OCCT cxx L389 (arch. diff. #2).
        let mut my_d_prism = BRepFillEvolved::new();
        my_d_prism.perform_with_face_spine(
            &my_spine,
            &my_profile,
            &Ax3::new(),
            GeomAbsJoinType::Arc,
            false,
        );

        // myMap is never written by this constructor body (the OCCT member
        // stays empty for the (Spine, Height, Angle) form).
        let my_map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        let mut my_res = Shape::null();
        let mut my_first_shape = Shape::null();
        let mut my_last_shape = Shape::null();

        if my_d_prism.is_done() {
            // OCCT cxx L393-397.
            let mut bs = LocOpeBuildShape::new();
            let c = b.make_compound(&mut pool, Vec::new());
            let mut lfaces: Vec<Shape> = Vec::new();
            let mut lcomplete: Vec<Shape> = Vec::new();

            // OCCT cxx L401-402.
            let exp_s: Vec<Shape> = explorer(&my_spine, ShapeType::Edge, ShapeType::Shape);
            let mut view: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
            // OCCT cxx L403-414.
            for es in &exp_s {
                let lffs = my_d_prism.generated_shapes(es, &my_profile1);
                for value in &lffs {
                    if view.insert(shape_key(value)) {
                        b.add_to_compound(&mut pool, c.clone(), value.clone());
                    }
                }
            }

            // OCCT cxx L416-420: theMapEF over C.
            let mut the_map_ef: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            map_shapes_and_ancestors(
                &c,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_map_ef,
            );
            view.clear();

            // OCCT cxx L423-456.
            for i in 1..=the_map_ef.len() {
                let (_, entry) = the_map_ef.get_index(i - 1).expect("theMapEF entry");
                // OCCT cxx L425.
                if entry.1.len() == 1 {
                    // OCCT cxx L427-428.
                    let edg = entry.0.clone();
                    let fac = entry.1.first().expect("First()").clone();
                    // OCCT cxx L429.
                    if view.insert(shape_key(&fac)) {
                        // OCCT cxx L431-432 (arch. diff. #6).
                        let new_face = pool.empty_copied(&fac);
                        // OCCT cxx L434-453.
                        let fac_forward = with_orientation(&fac, Orientation::Forward);
                        'wires: for wire in
                            explorer(&fac_forward, ShapeType::Wire, ShapeType::Shape)
                        {
                            // OCCT cxx L438.
                            let mut hit = false;
                            for e2 in explorer(&wire, ShapeType::Edge, ShapeType::Shape) {
                                // OCCT cxx L441.
                                if e2.is_same(&edg) {
                                    // OCCT cxx L443-446.
                                    let fd = pool.face_mut(new_face.clone());
                                    fd.outer_wire = wire.clone();
                                    lfaces.push(new_face.clone());
                                    lcomplete.push(new_face.clone());
                                    hit = true;
                                    break;
                                }
                            }
                            // OCCT cxx L449-452.
                            if hit {
                                break 'wires;
                            }
                        }
                    }
                }
            }

            // OCCT cxx L458-459.
            bs.perform(&lfaces);
            my_first_shape = bs.shape().cloned().unwrap_or_else(Shape::null);

            // OCCT cxx L461: B.MakeCompound(D).
            let d = b.make_compound(&mut pool, Vec::new());

            // OCCT cxx L463-464.
            view.clear();

            // OCCT cxx L466-477.
            for es in &exp_s {
                let lfls = my_d_prism.generated_shapes(es, &my_profile3);
                for value in &lfls {
                    if view.insert(shape_key(value)) {
                        b.add_to_compound(&mut pool, d.clone(), value.clone());
                    }
                }
            }

            // OCCT cxx L479-482.
            lfaces.clear();
            the_map_ef.clear();
            map_shapes_and_ancestors(
                &d,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_map_ef,
            );
            view.clear();

            // OCCT cxx L484-517.
            for i in 1..=the_map_ef.len() {
                let (_, entry) = the_map_ef.get_index(i - 1).expect("theMapEF entry");
                // OCCT cxx L486.
                if entry.1.len() == 1 {
                    // OCCT cxx L488-489.
                    let edg = entry.0.clone();
                    let fac = entry.1.first().expect("First()").clone();
                    // OCCT cxx L490.
                    if view.insert(shape_key(&fac)) {
                        // OCCT cxx L492-493 (arch. diff. #6).
                        let new_face = pool.empty_copied(&fac);
                        // OCCT cxx L495-514.
                        let fac_forward = with_orientation(&fac, Orientation::Forward);
                        'wires: for wire in
                            explorer(&fac_forward, ShapeType::Wire, ShapeType::Shape)
                        {
                            // OCCT cxx L499.
                            let mut hit = false;
                            for e2 in explorer(&wire, ShapeType::Edge, ShapeType::Shape) {
                                // OCCT cxx L502.
                                if e2.is_same(&edg) {
                                    // OCCT cxx L504-507.
                                    let fd = pool.face_mut(new_face.clone());
                                    fd.outer_wire = wire.clone();
                                    lfaces.push(new_face.clone());
                                    lcomplete.push(new_face.clone());
                                    hit = true;
                                    break;
                                }
                            }
                            // OCCT cxx L510-513.
                            if hit {
                                break 'wires;
                            }
                        }
                    }
                }
            }
            // OCCT cxx L518-519.
            bs.perform(&lfaces);
            my_last_shape = bs.shape().cloned().unwrap_or_else(Shape::null);

            // OCCT cxx L521-546.
            view.clear();
            for es in &exp_s {
                let ls = my_d_prism.generated_shapes(es, &my_profile2);
                for value in &ls {
                    // OCCT cxx L528-531.
                    if view.insert(shape_key(value)) {
                        lcomplete.push(value.clone());
                    }
                }
                // OCCT cxx L533-545.
                for es2 in explorer(es, ShapeType::Vertex, ShapeType::Shape) {
                    let ls2 = my_d_prism.generated_shapes(&es2, &my_profile2);
                    for value in &ls2 {
                        if view.insert(shape_key(value))
                            && value.shape_type() == ShapeType::Face
                        {
                            lcomplete.push(value.clone());
                        }
                    }
                }
            }

            // OCCT cxx L548-550.
            bs.perform(&lcomplete);
            my_res = bs.shape().cloned().unwrap_or_else(Shape::null);
            brep_lib_update_tolerances(&mut my_res);
        }

        LocOpeDPrism {
            my_d_prism,
            my_res,
            my_spine,
            my_profile,
            my_profile1,
            my_profile2,
            my_profile3,
            my_height,
            my_first_shape,
            my_last_shape,
            my_curvs: Vec::new(),
            my_map,
        }
    }

    /// OCCT LocOpe_DPrism::IsDone() (cxx L556-559).
    pub fn is_done(&self) -> bool {
        self.my_d_prism.is_done()
    }

    /// OCCT LocOpe_DPrism::Shape() (cxx L563-570).
    pub fn shape(&self) -> Shape {
        if !self.my_d_prism.is_done() {
            // OCCT cxx L565-568: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.clone()
    }

    /// OCCT LocOpe_DPrism::Spine() (cxx L574-577).
    pub fn spine(&self) -> Shape {
        self.my_spine.clone()
    }

    /// OCCT LocOpe_DPrism::Profile() (cxx L581-584).
    pub fn profile(&self) -> Shape {
        self.my_profile.clone()
    }

    /// OCCT LocOpe_DPrism::FirstShape() (cxx L588-591).
    pub fn first_shape(&self) -> Shape {
        self.my_first_shape.clone()
    }

    /// OCCT LocOpe_DPrism::LastShape() (cxx L595-598).
    pub fn last_shape(&self) -> Shape {
        self.my_last_shape.clone()
    }

    /// OCCT LocOpe_DPrism::Shapes(S) (cxx L602-616) — the OCCT returns a
    /// reference (bound to a temporary in the GeneratedShapes branch); rcad
    /// returns the owned list.
    pub fn shapes(&self, s: &Shape) -> Vec<Shape> {
        if !self.my_d_prism.is_done() {
            // OCCT cxx L604-607: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        // OCCT cxx L608-611.
        if let Some(list) = self.my_map.get(&shape_key(s)) {
            return list.clone();
        }
        // OCCT cxx L614.
        self.my_d_prism.generated_shapes(s, &self.my_profile2)
    }

    /// OCCT LocOpe_DPrism::Curves(SCurves) (cxx L620-690).
    pub fn curves(&self, s_curves: &mut Vec<Curve3>) {
        // Retrieves dy and dz with myProfile2
        //
        // OCCT cxx L623-628.
        let (v1_opt, v2_opt) = crate::feat::loc_ope::top_exp_vertices_edge(&self.my_profile);
        let v1 = v1_opt.expect("V1");
        let v2 = v2_opt.expect("V2");
        let p1 = brep_tool_pnt(&v1);
        let p2 = brep_tool_pnt(&v2);
        let dy = p2.y - p1.y;
        let dz = p2.z - p1.z;
        // OCCT cxx L629.
        s_curves.clear();
        // OCCT cxx L630-634: S = BRep_Tool::Surface(mySpine);
        // RectangularTrimmedSurface -> BasisSurface.
        let mut s = brep_tool_surface(&self.my_spine);
        if let Surface3::Trimmed(ts) = &s {
            s = (*ts.basis).clone();
        }

        // OCCT cxx L636-640.
        let pp = match &s {
            Surface3::Plane(pp) => *pp,
            _ => {
                // OCCT cxx L638-639: throw Standard_ConstructionError().
                panic!("Standard_ConstructionError");
            }
        };

        // OCCT cxx L642-647: Normale = P.Axis().Direction(); the rcad Plane
        // is always direct (arch. diff. #8), so the
        // "if (!P.Direct()) Normale.Reverse()" has no rcad false case.
        let normale = pp.normal;

        // OCCT cxx L649-654.
        let mut the_map: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();
        let spine_forward = with_orientation(&self.my_spine, Orientation::Forward);

        // OCCT cxx L656-689.
        for cur_edg in explorer(&spine_forward, ShapeType::Edge, ShapeType::Shape) {
            let edg = cur_edg.clone();
            // OCCT cxx L659-662.
            if !the_map.insert(shape_key(&edg)) {
                continue;
            }
            // OCCT cxx L663.
            if !brep_tool_degenerated(&edg) {
                // OCCT cxx L665-666: BRep_Tool::Curve + Transformed (the
                // loc_ope_find_edges arch. diff. #1 vehicle).
                let Some((c, f, l)) = crate::feat::loc_ope_find_edges::brep_tool_curve(&edg)
                else {
                    continue;
                };
                let c: Curve3 = c;
                // OCCT cxx L667-668.
                let u1 = -2. * self.my_height.abs();
                let u2 = 2. * self.my_height.abs();

                // OCCT cxx L670-687.
                for i in 0..=NECHANT {
                    let prm = ((NECHANT - i) as f64 * f + i as f64 * l) / (NECHANT as f64);
                    let pt = c.point_at(prm);
                    let mut d1 = c.derivative_at(prm);
                    if cur_edg.orientation == Orientation::Reversed {
                        d1 = -d1;
                    }
                    d1 = d1.normalize_or_zero();
                    let locy = normale.cross(d1).normalize_or_zero();
                    let ldir = locy * dy + normale * dz;
                    let lin = rcad_kernel::math::gp::Lin::from_pnt_dir(pt, ldir);
                    let lin = Curve3::Line(Line3::new(lin.pos, lin.dir));
                    let trlin = Curve3::Trimmed(TrimmedCurve3::new(lin, u1, u2));
                    s_curves.push(trlin);
                }
            }
        }
    }

    /// OCCT LocOpe_DPrism::BarycCurve() (cxx L694-753).
    pub fn baryc_curve(&self) -> Curve3 {
        // OCCT cxx L696-700.
        let (v1_opt, v2_opt) = crate::feat::loc_ope::top_exp_vertices_edge(&self.my_profile);
        let v1 = v1_opt.expect("V1");
        let v2 = v2_opt.expect("V2");
        let p1 = brep_tool_pnt(&v1);
        let p2 = brep_tool_pnt(&v2);
        let dz = p2.z - p1.z;

        // OCCT cxx L702-706.
        let mut s = brep_tool_surface(&self.my_spine);
        if let Surface3::Trimmed(ts) = &s {
            s = (*ts.basis).clone();
        }

        // OCCT cxx L708-719.
        let pp = match &s {
            Surface3::Plane(pp) => *pp,
            _ => {
                // OCCT cxx L710-711: throw Standard_ConstructionError().
                panic!("Standard_ConstructionError");
            }
        };
        // OCCT cxx L714-719 (arch. diff. #8).
        let normale = pp.normal;
        // OCCT cxx L720-731: the OCCT_DEBUG trace block over the REVERSED
        // spine (the Normale.Reverse() is commented out in the OCCT source)
        // is dead source text — arch. diff. #9.

        // OCCT cxx L732.
        let vec = normale * dz;

        // OCCT cxx L734-743.
        let mut bar = DVec3::ZERO;
        let mut spt: Vec<DVec3> = Vec::new();
        if !self.my_first_shape.is_null() {
            crate::feat::loc_ope::sample_edges(&self.my_first_shape, &mut spt);
        } else {
            crate::feat::loc_ope::sample_edges(&self.my_spine, &mut spt);
        }
        for pvt in &spt {
            // OCCT cxx L746-747: bar.ChangeCoord() += pvt.XYZ().
            bar += *pvt;
        }
        // OCCT cxx L749.
        bar /= spt.len() as f64;
        // OCCT cxx L750-752: newAx(bar, Vec); Geom_Line(newAx).
        let new_ax = Ax1::new(bar, vec);
        let the_lin = Curve3::Line(Line3::new(new_ax.location, new_ax.direction));
        the_lin
    }
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
