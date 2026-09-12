// OCCT BRepFeat_RibSlot.cxx L747-2668 — 1:1 translation (part b).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_RibSlot.cxx
//
// This file continues the `impl BRepFeatRibSlot` of brep_feat_rib_slot.rs
// (the 2000-line rule; architecture difference #1 of the parent module
// header). Content: ExtremeFaces (L747-1319), PtOnEdgeVertex (L1327-1411),
// SlidingProfile (L1418-1786), NoSlidingProfile (L1793-2668).
//
// The architecture-difference numbering of the parent module header applies
// (the vehicle re-hosts — IntCurvesFaceIntersector / Geom2dAPIInterCurveCurve
// / GeomAPIProjectPointOnCurve / FClass2d / BRepLib_MakeEdge family /
// geom_curve_reversed / gp_Dir comparisons — live in the parent file).

use crate::feat::brep_feat_form_2::{brep_algo_is_valid, CutVehicle};
use crate::feat::brep_feat_rib_slot::{
    brep_top_adaptor_fclass2d_perform, brep_tool_curve, brep_tool_degenerated,
    brep_tool_is_closed, brep_tool_pnt, brep_tool_surface, brep_tool_tolerance, data_map_bind,
    data_map_change_find, data_map_is_bound, geom_api_to_2d, geom_curve_reversed, geom_line_parts,
    make_edge_c_p_p, make_edge_c_v_v, make_edge_p_p, make_edge_v_v, make_vertex, map_add,
    shape_is_same, shape_key, top_exp_first_vertex, top_exp_last_vertex, BRepFeatRibSlot,
    Geom2dAPIInterCurveCurve, GeomAPIProjectPointOnCurve,
};
use crate::feat::loc_ope_cs_intersector::{IntCurvesFaceIntersector, LocOpeCSIntersector};
use crate::feat::loc_ope_find_edges::elclib_parameter_lin;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval, Line3, Plane, Surface3, TrimmedCurve3};
use rcad_kernel::topo::topods::{tshape_flags, BRep, BRepBuilder};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use rcad_kernel::precision::CONFUSION;

/// OCCT gp_Dir::IsEqual(Other, AngularTolerance) — the angle between the two
/// directions is within the tolerance (gp_Dir.cxx: A.Angle(B) <= Tolerance,
/// i.e. the dot product is >= cos(Tolerance)).
fn gp_dir_is_equal(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    a.dot(b) >= ang_tol.cos()
}

/// OCCT gp_Dir::IsOpposite(Other, AngularTolerance) — the angle is
/// PI - Tolerance or more (the dot product <= -cos(Tolerance)).
fn gp_dir_is_opposite(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    a.dot(b) <= -ang_tol.cos()
}

/// The edges of a wire in wire order — the BRepTools_WireExplorer re-host
/// (the rcad wire carries the traversal order; the same reduction as
/// topalgo::thru_sections wire_edges).
fn wire_explorer_edges(the_wire: &Shape) -> Vec<Shape> {
    match the_wire.data.as_ref() {
        rcad_kernel::topo::topods::TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT BRepLib_MakeFace myPln->Pln() + wire (BRepLib_MakeFace.cxx L262-272:
/// `Init(Pl, false, Confusion); Add(W); if (Inside && BRep_Tool::IsClosed(W))
/// CheckInside();`). The caller always passes Inside = true
/// (BRepFeat_RibSlot.cxx L1733 `BRepLib_MakeFace f(myPln->Pln(), WW, true)`).
fn make_face_plane_wire(
    the_b: &mut BRepBuilder,
    the_pool: &mut BRep,
    the_pln: &Plane,
    the_wire: &Shape,
) -> Shape {
    let fac = the_b.make_face(the_pool, Some(Surface3::Plane(*the_pln)), the_wire.clone());
    // OCCT L267-270.
    if brep_tool_is_closed(the_wire) {
        check_inside(the_pool, &fac);
    }
    fac
}

/// OCCT BRepLib_MakeFace::CheckInside (BRepLib_MakeFace.cxx L905-924):
/// "Reverses the current face if not a bounded area" — when the infinite
/// point classifies IN, `myShape` is replaced by an `EmptyCopied` face whose
/// children are added REVERSED. `BRep_TFace::EmptyCopy` (BRep_TFace.cxx
/// L36-43) copies the Surface/Location/Tolerance and no children, so the net
/// effect is that every wire of the face is reversed and the surface kept.
pub(crate) fn check_inside(the_pool: &mut BRep, the_fac: &Shape) {
    let Some(surf) = brep_tool_surface(the_fac) else {
        return;
    };
    let locations = [glam::DAffine3::IDENTITY];
    let src = crate::topalgo::shape_source::FaceShapeSource::new(the_fac, surf, &locations);
    // BRepTopAdaptor_FClass2d FClass(F, 0.).
    let a_cl = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(&src, 0, 0.0);
    if a_cl.perform_infinite_point(&src) != rcad_kernel::topods::State::In {
        return;
    }
    // B.Add(S, it.Value().Reversed()) for every child (the wires).
    let fd = the_pool.face_mut(the_fac.clone());
    fd.outer_wire.orientation = reverse_orientation(fd.outer_wire.orientation);
    for w in fd.inner_wires.iter_mut() {
        w.orientation = reverse_orientation(w.orientation);
    }
}

/// OCCT TopAbs::Reverse(O) — FORWARD <-> REVERSED, INTERNAL <-> EXTERNAL
/// (TopAbs.cxx L36-52).
fn reverse_orientation(o: rcad_kernel::topods::Orientation) -> rcad_kernel::topods::Orientation {
    use rcad_kernel::topods::Orientation;
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// OCCT BRepLib_MakeFace(const TopoDS_Wire& W, const bool OnlyPlane = false)
/// (BRepLib_MakeFace.cxx L189-262) — the face carrying the surface FOUND
/// through the wire, not a surface-less face:
///
/// ```text
///   BRepLib_FindSurface FS(W, -1, OnlyPlane, true);
///   if (!FS.Found()) { myError = BRepLib_NotPlanar; return; }
///   double tol = std::max(1.2 * FS.ToleranceReached(), FS.Tolerance());
///   B.MakeFace(F, FS.Surface(), FS.Location(), tol);
///   B.Add(F, W);
/// ```
///
/// With `OnlyPlane == false` the wire is added unchanged (the degenerate-edge
/// filtering of cxx L209-260 belongs to the `OnlyPlane == true` arm).
/// BRepFeat_RibSlot relies on this: NoSlidingProfile cxx L2643 and
/// SlidingProfile cxx L1764 both call `BRepLib_MakeFace ff(ww)` on the
/// boolean-cut wire and assert `BRepAlgo::IsValid` on the result.
fn make_face_wire(the_b: &mut BRepBuilder, the_pool: &mut BRep, the_wire: &Shape) -> Shape {
    // The wire reaches here from the CutVehicle result pool (OCCT carries the
    // TShape graph by pointer, rcad's `Shape::index` is pool-local), so adopt
    // its subgraph into `the_pool` before any pool-based lookup happens.
    let the_wire = &crate::feat::brep_feat_form_2::adopt_subgraph_into(the_pool, the_wire);
    // BRepLib_FindSurface FS(W, -1, OnlyPlane = false, OnlyClosed = true).
    let mut fs =
        crate::topalgo::brep_lib_find_surface::BRepLibFindSurface::new_closed(
            the_pool, the_wire, -1.0, false, true,
        );
    if !fs.found() {
        // OCCT: myError = BRepLib_NotPlanar; the shape stays null.
        return Shape::null();
    }
    // OCCT L204: tol = max(1.2 * ToleranceReached(), Tolerance()).
    let tol = (1.2 * fs.tolerance_reached()).max(fs.tolerance());
    let surface = fs.surface();
    let location = fs.location();
    // OCCT L206: B.MakeFace(F, FS.Surface(), FS.Location(), tol); L262: B.Add(F, W).
    let fac = the_b.make_face(the_pool, surface, the_wire.clone());
    let fd = the_pool.face_mut(fac.clone());
    fd.surface_location = location;
    fd.tolerance = tol;
    fac
}

/// OCCT w.Closed(BRep_Tool::IsClosed(w)) — the flag is SET to the value
/// (set or cleared).
fn wire_set_closed(the_pool: &mut BRep, the_wire: &Shape, the_closed: bool) {
    let wd = the_pool.wire_mut(the_wire.clone());
    if the_closed {
        wd.flags |= tshape_flags::CLOSED;
    } else {
        wd.flags &= !tshape_flags::CLOSED;
    }
}

impl BRepFeatRibSlot {
    // -----------------------------------------------------------------
    // ExtremeFaces (cxx L747-1319) — the base faces of the rib.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::ExtremeFaces (cxx L747-1319). The
    // OCCT locals (Data and the explorers) are initialized and reassigned
    // along the OCCT control flow — the never-read initial assignments are
    // the OCCT shape (kept with an allow).
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub(crate) fn extreme_faces(
        &self,
        revol_rib: bool,
        bnd: f64,
        pln: &Plane,
        first_edge: &mut Shape,
        last_edge: &mut Shape,
        first_face: &mut Shape,
        last_face: &mut Shape,
        first_vertex: &mut Shape,
        last_vertex: &mut Shape,
        on_first_face: &mut bool,
        on_last_face: &mut bool,
        pt_on_first_edge: &mut bool,
        pt_on_last_edge: &mut bool,
        on_first_edge: &mut Shape,
        on_last_edge: &mut Shape,
    ) -> bool {
        // OCCT L769.
        let mut data = true;
        // OCCT L770-777.
        *first_face = Shape::null();
        *last_face = Shape::null();
        *first_edge = Shape::null();
        *last_edge = Shape::null();
        *pt_on_first_edge = false;
        *pt_on_last_edge = false;
        *on_first_edge = Shape::null();
        *on_last_edge = Shape::null();

        // OCCT L779-781: BRepIntCurveSurface_Inter inter (the per-use
        // IntCurvesFaceIntersector re-host); BRep_Builder B; TopExp_Explorer
        // ex1 — the explorer re-host materializes per iteration.
        let mut pool = BRep::new();

        // OCCT L783.
        let mut first_ok = false;
        let mut last_ok = false;

        // OCCT L785-791.
        let wire_edges = explorer_edges_of(&self.my_wire);
        let number_of_edges = wire_edges.len() as i32;

        // --- the wire includes only one edge (OCCT L794-1075) ---
        if number_of_edges == 1 {
            // OCCT L800-811: exp.ReInit(); E = the unique edge;
            // cc = BRep_Tool::Curve(E, f, l); p1/p2 = the limit points;
            // FirstPar = f; LastPar = l.
            let e = wire_edges[0].clone();
            let (cc, f, l) = match brep_tool_curve(&e) {
                Some(hit) => hit,
                None => panic!("null 3D curve (OCCT null-handle deref)"),
            };
            let p1 = brep_tool_pnt(&top_exp_first_vertex(&e, true));
            let p2 = brep_tool_pnt(&top_exp_last_vertex(&e, true));

            let first_par = f;
            let last_par = l;

            // OCCT L813-817: find if the 2 points limiting the unique edge of
            // the wire are on an edge or a vertex of the base shape.
            let mut pt_on_first_vertex = false;
            let mut pt_on_last_vertex = false;
            let mut on_first_vertex = Shape::null();
            let mut on_last_vertex = Shape::null();
            self.pt_on_edge_vertex(
                revol_rib,
                &self.my_sbase,
                p1,
                first_vertex,
                last_vertex,
                pt_on_first_edge,
                on_first_edge,
                &mut pt_on_first_vertex,
                &mut on_first_vertex,
            );
            self.pt_on_edge_vertex(
                revol_rib,
                &self.my_sbase,
                p2,
                first_vertex,
                last_vertex,
                pt_on_last_edge,
                on_last_edge,
                &mut pt_on_last_vertex,
                &mut on_last_vertex,
            );

            // OCCT L837.
            let mut map = crate::feat::brep_feat_builder::OcctShapeMap::new();

            // OCCT L839-916.
            if *pt_on_first_edge {
                if !pt_on_first_vertex {
                    // OCCT L843-875: find FirstFace : face of the base shape
                    // containing OnFirstEdge meeting ChoiceOfFaces.
                    let mut faces: Vec<Shape> = Vec::new();
                    faces.clear();
                    map.clear();
                    for ex4 in explorer_edges_of_faces(&self.my_sbase) {
                        let fx = ex4;
                        if !map_add(&mut map, &fx) {
                            continue;
                        }
                        for ex5 in explorer_edges_of(&fx) {
                            let ee = ex5;
                            if shape_is_same(&ee, on_first_edge) {
                                faces.push(fx.clone());
                            }
                        }
                    }
                    if !faces.is_empty() {
                        let fff = Self::choice_of_faces(
                            &mut faces,
                            &cc,
                            first_par + bnd / 50.0,
                            bnd / 50.0,
                            pln,
                        );
                        if !fff.is_null() {
                            *first_face = fff;
                        }
                    }
                } else {
                    // OCCT L876-911: pt_on_first_vertex — find FirstFace :
                    // face of the base shape containing OnFirstVertex
                    // meeting ChoiceOfFaces.
                    let mut faces: Vec<Shape> = Vec::new();
                    faces.clear();
                    map.clear();
                    for ex4 in explorer_edges_of_faces(&self.my_sbase) {
                        let fx = ex4;
                        if !map_add(&mut map, &fx) {
                            continue;
                        }
                        for ex5 in explorer_vertices_of(&fx) {
                            let vv = ex5;
                            if shape_is_same(&vv, &on_first_vertex) {
                                faces.push(fx.clone());
                                break;
                            }
                        }
                    }
                    if !faces.is_empty() {
                        let fff = Self::choice_of_faces(
                            &mut faces,
                            &cc,
                            first_par + bnd / 50.0,
                            bnd / 50.0,
                            pln,
                        );
                        if !fff.is_null() {
                            *first_face = fff;
                        }
                    }
                }
                // OCCT L912-915.
                *first_edge = e.clone();
                *first_vertex = make_vertex(&mut pool, p1);
                *on_first_face = true;
            }

            // OCCT L918-996.
            if *pt_on_last_edge {
                if !pt_on_last_vertex {
                    // OCCT L922-955: find LastFace : face of the base shape
                    // containing OnLastEdge meeting ChoiceOfFaces.
                    let mut faces: Vec<Shape> = Vec::new();
                    faces.clear();
                    map.clear();
                    for ex4 in explorer_edges_of_faces(&self.my_sbase) {
                        let fx = ex4;
                        if !map_add(&mut map, &fx) {
                            continue;
                        }
                        for ex5 in explorer_edges_of(&fx) {
                            let ee = ex5;
                            if shape_is_same(&ee, on_last_edge) {
                                faces.push(fx.clone());
                                break;
                            }
                        }
                    }
                    if !faces.is_empty() {
                        let fff = Self::choice_of_faces(
                            &mut faces,
                            &cc,
                            last_par - bnd / 50.0,
                            bnd / 50.0,
                            pln,
                        );
                        if !fff.is_null() {
                            *last_face = fff;
                        }
                    }
                } else {
                    // OCCT L956-991: pt_on_last_edge && pt_on_last_vertex —
                    // find LastFace : face of the base shape containing
                    // OnLastVertex meeting ChoiceOfFaces.
                    let mut faces: Vec<Shape> = Vec::new();
                    faces.clear();
                    map.clear();
                    for ex4 in explorer_edges_of_faces(&self.my_sbase) {
                        let fx = ex4;
                        if !map_add(&mut map, &fx) {
                            continue;
                        }
                        for ex5 in explorer_vertices_of(&fx) {
                            let vv = ex5;
                            if shape_is_same(&vv, &on_last_vertex) {
                                faces.push(fx.clone());
                                break;
                            }
                        }
                    }
                    if !faces.is_empty() {
                        let fff = Self::choice_of_faces(
                            &mut faces,
                            &cc,
                            last_par - bnd / 50.0,
                            bnd / 50.0,
                            pln,
                        );
                        if !fff.is_null() {
                            *last_face = fff;
                        }
                    }
                }
                // OCCT L992-995.
                *last_edge = e.clone();
                *last_vertex = make_vertex(&mut pool, p2);
                *on_last_face = true;
            }

            // OCCT L998-1001.
            if !first_face.is_null() && !last_face.is_null() {
                return true;
            }

            // --- FirstFace or LastFace was not found (OCCT L1003-1074) ---
            let mut asi = LocOpeCSIntersector::with_shape(&self.my_sbase);
            let mut scur: Vec<Option<Curve3>> = Vec::new();
            scur.clear();
            scur.push(Some(cc.clone()));
            asi.perform_cur(&scur);
            // OCCT L802: firstpoint, lastpoint (the OCCT locals are assigned
            // in the IsDone branch below; the zero init carries the
            // uninitialized read note).
            let mut firstpoint = DVec3::ZERO;
            let mut lastpoint = DVec3::ZERO;
            if asi.is_done() && asi.nb_points(1) >= 2 {
                // OCCT L1016-1032.
                let mut lastpar = asi.point(1, asi.nb_points(1)).parameter();
                let mut lastindex = asi.nb_points(1);
                if lastpar > l {
                    for jj in (1..=(asi.nb_points(1) - 1)).rev() {
                        let par = asi.point(1, jj).parameter();
                        if par <= l {
                            lastpar = par;
                            lastindex = jj;
                            break;
                        }
                    }
                }
                let firstindex = lastindex - 1;
                let firstpar = asi.point(1, firstindex).parameter();

                // OCCT L1034-1041.
                if first_face.is_null() {
                    // the intersector points always carry the face.
                    *first_face = asi
                        .point(1, firstindex)
                        .face()
                        .cloned()
                        .unwrap_or_else(Shape::null);
                    firstpoint = cc.point_at(firstpar);
                    *first_vertex = make_vertex(&mut pool, firstpoint);
                    *first_edge = e.clone();
                }

                // OCCT L1043-1050.
                if last_face.is_null() {
                    *last_face = asi
                        .point(1, lastindex)
                        .face()
                        .cloned()
                        .unwrap_or_else(Shape::null);
                    lastpoint = cc.point_at(lastpar);
                    *last_vertex = make_vertex(&mut pool, lastpoint);
                    *last_edge = e.clone();
                }
            } else {
                // OCCT L1052-1060: less than 2 intersection points.
                data = false;
                return data;
            }

            // OCCT L1062-1070.
            if !*on_first_face {
                *on_first_face = p1.distance(firstpoint) <= CONFUSION;
            }

            if !*on_last_face {
                *on_last_face = p2.distance(lastpoint) <= CONFUSION;
            }

            // OCCT L1072-1074.
            data = !((*first_face).is_null() || (*last_face).is_null());

            return data;
        }
        // --- the wire consists of several edges (OCCT L1077-1318) ---
        else {
            // OCCT L1083: BRepTools_WireExplorer ex(myWire).
            for e in wire_explorer_edges(&self.my_wire) {
                // OCCT L1086-1092.
                let (cur, mut f, mut l) = match brep_tool_curve(&e) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                f = f - bnd / 10000.0;
                l = l + bnd / 10000.0;
                let curve = Curve3::Trimmed(TrimmedCurve3::new(cur.clone(), f, l));
                // OCCT L1097: P2 = BRep_Tool::Pnt(TopExp::LastVertex(E, true)).
                let p2 = brep_tool_pnt(&top_exp_last_vertex(&e, true));
                // OCCT L1098-1106.
                let mut the_vertex = Shape::null();
                let mut the_edge = Shape::null();
                let mut the_face = Shape::null();
                let mut pt_on_edge = false;
                let mut pt_on_vertex = false;
                let mut on_edge = Shape::null();
                let mut on_vertex = Shape::null();
                let mut intpar = 0.0;
                // OCCT L1107: for (; ex1.More(); ex1.Next()) — the faces of
                // mySbase.
                for a_cur_face in explorer_edges_of_faces(&self.my_sbase) {
                    // OCCT L1110-1115: GeomAdaptor_Curve aGAC(curve);
                    // inter.Init(aCurFace, aGAC, Tolerance(aCurFace));
                    // if (!inter.More()) continue.
                    let mut the_int = IntCurvesFaceIntersector::new(
                        &a_cur_face,
                        brep_tool_tolerance(&a_cur_face),
                    );
                    let dom = curve.default_domain();
                    the_int.perform_curve(&curve, dom[0], dom[1]);
                    if the_int.nb_pnt() == 0 {
                        continue;
                    }
                    // OCCT L1116: for (; inter.More(); inter.Next()).
                    'pt: for pt_j in 1..=the_int.nb_pnt() {
                        // OCCT L1118.
                        let the_point = the_int.pnt(pt_j);
                        // OCCT L1119-1126.
                        if !first_vertex.is_null() {
                            let point = brep_tool_pnt(first_vertex);
                            if point.distance(the_point) <= brep_tool_tolerance(&a_cur_face) {
                                continue 'pt;
                            }
                        }
                        // OCCT L1127-1130.
                        intpar = Self::int_par(&curve, the_point);
                        the_edge = e.clone();
                        the_face = a_cur_face.clone();
                        the_vertex = make_vertex(&mut pool, the_point);
                        // OCCT L1131-1137.
                        if !first_ok && the_point.distance(p2) <= CONFUSION {
                            continue 'pt;
                        }
                        // OCCT L1140-1148: find thepoint on an edge or a
                        // vertex of face f.
                        self.pt_on_edge_vertex(
                            revol_rib,
                            &a_cur_face,
                            the_point,
                            first_vertex,
                            last_vertex,
                            &mut pt_on_edge,
                            &mut on_edge,
                            &mut pt_on_vertex,
                            &mut on_vertex,
                        );

                        // OCCT L1152-1220: the First block.
                        if first_edge.is_null()
                            && !the_edge.is_null()
                            && !the_face.is_null()
                            && !the_vertex.is_null()
                        {
                            *first_edge = the_edge.clone();
                            *first_face = the_face.clone();
                            *first_vertex = the_vertex.clone();
                            *pt_on_first_edge = pt_on_edge;
                            *on_first_edge = on_edge.clone();
                            the_edge = Shape::null();
                            the_face = Shape::null();
                            the_vertex = Shape::null();
                            if pt_on_edge && !pt_on_vertex {
                                // OCCT L1164-1188.
                                let mut faces: Vec<Shape> = Vec::new();
                                faces.clear();
                                faces.push(first_face.clone());
                                for ex2 in explorer_edges_of_faces(&self.my_sbase) {
                                    let fx = ex2;
                                    for ex3 in explorer_edges_of(&fx) {
                                        let e3 = ex3;
                                        if shape_is_same(&e3, &on_edge)
                                            && !shape_is_same(&fx, first_face)
                                        {
                                            faces.push(fx.clone());
                                        }
                                    }
                                }
                                let fff = Self::choice_of_faces(
                                    &mut faces,
                                    &curve,
                                    intpar + bnd / 10.0,
                                    bnd / 10.0,
                                    pln,
                                );
                                if !fff.is_null() {
                                    *first_face = fff;
                                }
                            } else if pt_on_edge && pt_on_vertex {
                                // OCCT L1189-1215.
                                let mut faces: Vec<Shape> = Vec::new();
                                faces.clear();
                                faces.push(first_face.clone());
                                for ex2 in explorer_edges_of_faces(&self.my_sbase) {
                                    let fx = ex2;
                                    for ex3 in explorer_vertices_of(&fx) {
                                        let v = ex3;
                                        if shape_is_same(&v, &on_vertex)
                                            && !shape_is_same(&fx, first_face)
                                        {
                                            faces.push(fx.clone());
                                        }
                                    }
                                }
                                let fff = Self::choice_of_faces(
                                    &mut faces,
                                    &curve,
                                    intpar + bnd / 10.0,
                                    bnd / 10.0,
                                    pln,
                                );
                                if !fff.is_null() {
                                    *first_face = fff;
                                }
                            }
                            // OCCT L1216-1219.
                            if !first_edge.is_null()
                                && !first_face.is_null()
                                && !first_vertex.is_null()
                            {
                                first_ok = true;
                            }
                        }
                        // OCCT L1221-1288: the Last block.
                        if last_edge.is_null()
                            && !the_edge.is_null()
                            && !the_face.is_null()
                            && !the_vertex.is_null()
                            && !first_edge.is_null()
                        {
                            *last_edge = the_edge.clone();
                            *last_face = the_face.clone();
                            *last_vertex = the_vertex.clone();
                            *pt_on_last_edge = pt_on_edge;
                            *on_last_edge = on_edge.clone();
                            if pt_on_edge && !pt_on_vertex {
                                // OCCT L1229-1255.
                                let mut faces: Vec<Shape> = Vec::new();
                                faces.clear();
                                faces.push(last_face.clone());
                                for ex2 in explorer_edges_of_faces(&self.my_sbase) {
                                    let fx = ex2;
                                    for ex3 in explorer_edges_of(&fx) {
                                        let e3 = ex3;
                                        if shape_is_same(&e3, &on_edge)
                                            && !shape_is_same(&fx, last_face)
                                        {
                                            faces.push(fx.clone());
                                        }
                                    }
                                }
                                let fff = Self::choice_of_faces(
                                    &mut faces,
                                    &curve,
                                    intpar - bnd / 10.0,
                                    bnd / 10.0,
                                    pln,
                                );
                                if !fff.is_null() {
                                    *last_face = fff;
                                }
                            } else if pt_on_edge && pt_on_vertex {
                                // OCCT L1256-1282.
                                let mut faces: Vec<Shape> = Vec::new();
                                faces.clear();
                                faces.push(last_face.clone());
                                for ex2 in explorer_edges_of_faces(&self.my_sbase) {
                                    let fx = ex2;
                                    for ex3 in explorer_vertices_of(&fx) {
                                        let v = ex3;
                                        if shape_is_same(&v, &on_vertex)
                                            && !shape_is_same(&fx, last_face)
                                        {
                                            faces.push(fx.clone());
                                        }
                                    }
                                }
                                let fff = Self::choice_of_faces(
                                    &mut faces,
                                    &curve,
                                    intpar - bnd / 10.0,
                                    bnd / 10.0,
                                    pln,
                                );
                                if !fff.is_null() {
                                    *last_face = fff;
                                }
                            }
                            // OCCT L1283-1287.
                            if !last_edge.is_null()
                                && !last_face.is_null()
                                && !last_vertex.is_null()
                            {
                                last_ok = true;
                            }
                            break 'pt;
                        }
                    }
                }
            }

            // OCCT L1293-1318.
            if first_ok && last_ok {
                let pp1 = brep_tool_pnt(&top_exp_first_vertex(first_edge, true));
                let pp2 = brep_tool_pnt(&top_exp_last_vertex(last_edge, true));
                let p1 = brep_tool_pnt(first_vertex);
                let p2 = brep_tool_pnt(last_vertex);
                if p1.distance(pp1) <= brep_tool_tolerance(first_face) {
                    *on_first_face = true;
                }
                if p2.distance(pp2) <= brep_tool_tolerance(last_face) {
                    *on_last_face = true;
                }
                return true;
            } else {
                return false;
            }
        }
    }

    // -----------------------------------------------------------------
    // PtOnEdgeVertex (cxx L1327-1411).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::PtOnEdgeVertex (cxx L1327-1411) — find if the
    /// 2 limit points of the unique edge of a wire are on an edge or a
    /// vertex of the base shape. The FirstVertex/LastVertex parameters are
    /// unnamed (unused) in OCCT; the PtOnEdge/OnEdge/PtOnVertex/OnVertex
    /// out-parameters are NOT reset by the OCCT body (the reset block is
    /// commented out at cxx L1344-1347).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pt_on_edge_vertex(
        &self,
        revol_rib: bool,
        shape: &Shape,
        point: DVec3,
        _first_vertex: &Shape,
        _last_vertex: &Shape,
        pt_on_edge: &mut bool,
        on_edge: &mut Shape,
        pt_on_vertex: &mut bool,
        on_vertex: &mut Shape,
    ) {
        // OCCT L1343: TestOK.
        let mut test_ok;
        // OCCT L1349-1352: EXP over the edges; Map.
        let mut map = crate::feat::brep_feat_builder::OcctShapeMap::new();
        for exp in explorer_edges_of(shape) {
            let e = exp;
            // OCCT L1355-1358.
            if !map_add(&mut map, &e) {
                continue;
            }
            // OCCT L1359-1365.
            if !revol_rib && brep_tool_degenerated(&e) {
                continue;
            }
            // OCCT L1366-1371: fff, lll; ccc = BRep_Tool::Curve(e, fff, lll)
            // [trimmed when not RevolRib] — the null-handle continuation is
            // not defined in OCCT (GeomAPI_ProjectPointOnCurve raises
            // Standard_NoSuchObject on the null curve).
            let Some((ccc, fff, lll)) = brep_tool_curve(&e) else {
                panic!("Standard_NoSuchObject (null 3D curve)");
            };
            let ccc = if !revol_rib {
                Curve3::Trimmed(TrimmedCurve3::new(ccc, fff, lll))
            } else {
                ccc
            };
            // OCCT L1372-1373.
            let proj = GeomAPIProjectPointOnCurve::new(point, &ccc);
            test_ok = false;
            // OCCT L1374-1387.
            if !revol_rib {
                if proj.nb_points() == 1 {
                    test_ok = true;
                }
            } else {
                if proj.nb_points() >= 1 {
                    test_ok = true;
                }
            }
            // OCCT L1388-1409.
            if test_ok && proj.distance(1) <= brep_tool_tolerance(&e) {
                *pt_on_edge = true;
                *on_edge = e.clone();
                let ev1 = top_exp_first_vertex(&e, true);
                let ev2 = top_exp_last_vertex(&e, true);
                let ep1 = brep_tool_pnt(&ev1);
                let ep2 = brep_tool_pnt(&ev2);
                if point.distance(ep1) <= brep_tool_tolerance(&ev1) {
                    *pt_on_vertex = true;
                    *on_vertex = ev1;
                    break;
                } else if point.distance(ep2) <= brep_tool_tolerance(&ev1) {
                    // the OCCT source tests Tolerance(ev1) in the second arm
                    // (cxx L1402) — kept literally.
                    *pt_on_vertex = true;
                    *on_vertex = ev2;
                    break;
                }
                break;
            }
        }
    }

    // -----------------------------------------------------------------
    // SlidingProfile (cxx L1418-1786) — the profile face in case of sliding.
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::SlidingProfile (cxx L1418-1786). The
    /// FirstVertex/LastVertex parameters are unnamed (unused) in OCCT.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn sliding_profile(
        &self,
        prof: &mut Shape,
        revol_rib: bool,
        my_tol: f64,
        concavite: &mut i32,
        my_pln: &Plane,
        bnd_face: &Shape,
        check_pnt: DVec3,
        first_face: &Shape,
        last_face: &Shape,
        _first_vertex: &Shape,
        _last_vertex: &Shape,
        first_edge: &Shape,
        last_edge: &Shape,
    ) -> bool {
        // OCCT L1438.
        let mut profile_ok = true;
        // --case of sliding : construction of the wire of the profile
        // --> 1 part bounding box + 1 part wire
        //   attention to the compatibility of orientations

        // OCCT L1443-1444: FN, LN; BRepLib_MakeWire WW.
        let mut pool = BRep::new();
        let mut ww = BRepBuilder::new();
        let wire = ww.make_wire(&mut pool);

        // OCCT L1446-1456: FN/LN from the face normals; the groove (cut)
        // flips to stay in the material.
        let mut fn_dir = self.normal(first_face, self.my_first_pnt);
        let mut ln_dir = self.normal(last_face, self.my_last_pnt);

        if !self.my_fuse {
            fn_dir = -fn_dir;
            ln_dir = -ln_dir;
        }

        // OCCT L1458-1462: ln2 = Geom_Line(myFirstPnt, FN);
        // ln1 = Geom_Line(myLastPnt, LN); Pt.
        let ln2 = Curve3::Line(Line3::new(self.my_first_pnt, fn_dir));
        let ln1 = Curve3::Line(Line3::new(self.my_last_pnt, ln_dir));
        let mut pt = DVec3::ZERO;

        // OCCT L1464-1467: the 2d carriers and the intersection.
        let ln2d1 = match geom_api_to_2d(&ln1, my_pln) {
            Some(c) => c,
            None => panic!("GeomAPI::To2d null handle"),
        };
        let ln2d2 = match geom_api_to_2d(&ln2, my_pln) {
            Some(c) => c,
            None => panic!("GeomAPI::To2d null handle"),
        };

        let inter = Geom2dAPIInterCurveCurve::new(&ln2d1, &ln2d2, CONFUSION);

        // OCCT L1469-1493.
        let mut test_ok = true;
        if revol_rib {
            // OCCT L1472-1473: d1/d2 = ln->Position().Direction().
            let (_o1, d1) = geom_line_parts(&ln1);
            let (_o2, d2) = geom_line_parts(&ln2);
            if gp_dir_is_opposite(d1, d2, my_tol) {
                // OCCT L1477-1478: ElCLib::Parameter(ln->Lin(), myPnt).
                let par1 = elclib_parameter_lin(&curve3_lin(&ln1), self.my_first_pnt);
                let par2 = elclib_parameter_lin(&curve3_lin(&ln2), self.my_last_pnt);
                if par1 >= my_tol || par2 >= my_tol {
                    *concavite = 2; // parallel and concave
                    // OCCT L1482-1483.
                    let e1 = make_edge_p_p(&mut pool, self.my_last_pnt, self.my_first_pnt);
                    ww.add_to_wire(&mut pool, wire.clone(), e1);
                }
            }
            if gp_dir_is_equal(d1, d2, my_tol) {
                if *concavite == 3 {
                    test_ok = false;
                }
            }
        }

        // OCCT L1495-1507.
        if test_ok && inter.nb_points() > 0 {
            let p2d = inter.point(1);
            // OCCT L1500: myPln->D0(P.X(), P.Y(), Pt).
            pt = my_pln.origin + my_pln.u_dir * p2d.x + my_pln.v_dir * p2d.y;
            let par = Self::int_par(&ln1, pt);
            if par > 0.0 {
                *concavite = 1; // concave
            }
        }

        // --- construction of the profile face (OCCT L1509-1670) ---
        if *concavite == 1 {
            // OCCT L1510-1518: concave — extend the first and last edges of
            // the wire to the intersection point.
            let e1 = make_edge_p_p(&mut pool, self.my_last_pnt, pt);
            ww.add_to_wire(&mut pool, wire.clone(), e1);
            let e2 = make_edge_p_p(&mut pool, pt, self.my_first_pnt);
            ww.add_to_wire(&mut pool, wire.clone(), e2);
        } else if *concavite == 3 {
            // OCCT L1519-1669: BndEdge : edges of intersection with the
            // bounding box.
            let mut bnd_edge1 = Shape::null();
            let mut bnd_edge2 = Shape::null();
            let mut bnd_pnt1 = DVec3::ZERO;
            let mut bnd_pnt2 = DVec3::ZERO;
            // OCCT L1524: LastPnt (the OCCT local is assigned in the ring
            // walk below; the zero init carries the uninitialized read
            // note).
            let mut last_pnt = DVec3::ZERO;
            // OCCT L1525-1529: the first wire of BndFace.
            let bnd_wire = first_wire_of(bnd_face);
            let bnd_wire_edges = wire_explorer_edges(&bnd_wire);
            // OCCT BRepTools_WireExplorer explo(BndWire) (cxx L1567-1569): the
            // enumeration is the CONNECTIVITY walk, not the stored edge list.
            // The BndFace wire comes from the boolean Common, whose storage
            // order differs from the traversal order — verified against OCCT
            // on featrf_a1: storage BW[0..3] = z-/x+/z+/x-, traversal from
            // BW[0] = z- -> x- -> z+ -> x+. Walking the stored list descends
            // BW[1] (the edge carrying BndEdge1) instead of BW[0], so the
            // profile wire gains the two extra boundary edges and comes out
            // self-intersecting.
            let bnd_edges = match brep_tool_surface(bnd_face) {
                Some(surf) => {
                    let locations = [glam::DAffine3::IDENTITY];
                    let src =
                        crate::topalgo::shape_source::FaceShapeSource::new(bnd_face, surf, &locations);
                    crate::topalgo::brep_top_adaptor::fclass2d::wire_explorer_order(
                        &src,
                        0,
                        &bnd_wire_edges,
                    )
                }
                None => bnd_wire_edges,
            };
            // OCCT L1530-1569.
            for e in &bnd_edges {
                let (c, first, last) = match brep_tool_curve(e) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                let c2d = match geom_api_to_2d(&c, my_pln) {
                    Some(c) => c,
                    None => panic!("GeomAPI::To2d null handle"),
                };
                let intcln1 = Geom2dAPIInterCurveCurve::new(&ln2d1, &c2d, CONFUSION);
                if intcln1.nb_points() > 0 {
                    let p2d = intcln1.point(1);
                    let p = my_pln.origin + my_pln.u_dir * p2d.x + my_pln.v_dir * p2d.y;
                    let parl = Self::int_par(&ln1, p);
                    let parc = Self::int_par(&c, p);
                    if parc >= first && parc <= last && parl >= 0.0 {
                        bnd_edge1 = e.clone();
                        bnd_pnt1 = p;
                    }
                }

                let intcln2 = Geom2dAPIInterCurveCurve::new(&ln2d2, &c2d, CONFUSION);
                if intcln2.nb_points() > 0 {
                    let p2d = intcln2.point(1);
                    let p = my_pln.origin + my_pln.u_dir * p2d.x + my_pln.v_dir * p2d.y;
                    let parl = Self::int_par(&ln2, p);
                    let parc = Self::int_par(&c, p);
                    if parc >= first && parc <= last && parl >= 0.0 {
                        bnd_edge2 = e.clone();
                        bnd_pnt2 = p;
                    }
                }
                // OCCT L1565-1568.
                if !bnd_edge1.is_null() && !bnd_edge2.is_null() {
                    break;
                }
            }

            // OCCT L1571-1579.
            if bnd_edge1.is_null() || bnd_edge2.is_null() {

                profile_ok = false;
                return profile_ok;
            }

            // OCCT L1581-1582.
            let e1 = make_edge_p_p(&mut pool, self.my_last_pnt, bnd_pnt1);
            ww.add_to_wire(&mut pool, wire.clone(), e1);

            // OCCT L1584-1669.
            if shape_is_same(&bnd_edge1, &bnd_edge2) {
                // OCCT L1586-1591: same edge -> simply determined path.
                let e2 = make_edge_p_p(&mut pool, bnd_pnt1, bnd_pnt2);
                ww.add_to_wire(&mut pool, wire.clone(), e2);
                let e3 = make_edge_p_p(&mut pool, bnd_pnt2, self.my_first_pnt);
                ww.add_to_wire(&mut pool, wire.clone(), e3);
            } else {
                // OCCT L1594-1614: walk to BndEdge1.
                let mut idx = 0usize;
                while idx < bnd_edges.len() {
                    let e = &bnd_edges[idx];
                    if shape_is_same(e, &bnd_edge1) {
                        // OCCT L1600-1610.
                        let pp = brep_tool_pnt(&top_exp_last_vertex(e, true));
                        if pp.distance(bnd_pnt1) >= brep_tool_tolerance(e) {
                            last_pnt = pp;
                        }
                        let e2 = make_edge_p_p(&mut pool, bnd_pnt1, last_pnt);
                        ww.add_to_wire(&mut pool, wire.clone(), e2);
                        break;
                    }
                    idx += 1;
                }
                // OCCT L1616-1627: the wrapped advance (Next + the null
                // Current re-Init is the modulo wrap).
                if idx + 1 < bnd_edges.len() {
                    idx += 1;
                } else {
                    idx = 0;
                }

                // OCCT L1629-1665: check if this is BndEdge2 -> if yes the
                // path is closed, else add edges.
                let mut fin = false;
                while !fin {
                    let e = &bnd_edges[idx];
                    if !shape_is_same(e, &bnd_edge2) {
                        // OCCT L1638-1642.
                        let pp = brep_tool_pnt(&top_exp_last_vertex(e, true));
                        let ee = make_edge_p_p(&mut pool, last_pnt, pp);
                        ww.add_to_wire(&mut pool, wire.clone(), ee);
                        last_pnt = pp;
                    } else {
                        // OCCT L1646-1651: the path is closed -> since met
                        // BndEdge2, end of borders on BndFace.
                        fin = true;
                        let ee = make_edge_p_p(&mut pool, last_pnt, bnd_pnt2);
                        ww.add_to_wire(&mut pool, wire.clone(), ee);
                        last_pnt = bnd_pnt2;
                    }
                    // OCCT L1653-1664: the wrapped advance.
                    if idx + 1 < bnd_edges.len() {
                        idx += 1;
                    } else {
                        idx = 0;
                    }
                }

                // OCCT L1667-1668.
                let e3 = make_edge_p_p(&mut pool, bnd_pnt2, self.my_first_pnt);
                ww.add_to_wire(&mut pool, wire.clone(), e3);
            }
        }

        // --- construction of the profile (OCCT L1672-1731) ---
        // Explore the wire provided by the user (BRepTools_WireExplorer :
        // correct order - without repetition <> TopExp : non ordered).
        // OCCT L1678-1679.
        let (first_curve, _ff, _ll) = match brep_tool_curve(first_edge) {
            Some(hit) => hit,
            None => panic!("null 3D curve (OCCT null-handle deref)"),
        };

        if !shape_is_same(first_edge, last_edge) {
            // OCCT L1683-1686.
            let fl_vert = top_exp_last_vertex(first_edge, true);
            let fl_pnt = brep_tool_pnt(&fl_vert);
            let ef = make_edge_c_p_p(&mut pool, &first_curve, self.my_first_pnt, fl_pnt);
            ww.add_to_wire(&mut pool, wire.clone(), ef);
            // OCCT L1687-1694: skip the wire edges up to FirstEdge.
            let user_edges = wire_explorer_edges(&self.my_wire);
            let mut ex_i = 0usize;
            while ex_i < user_edges.len() {
                let e = &user_edges[ex_i];
                if shape_is_same(e, first_edge) {
                    break;
                }
                ex_i += 1;
            }
            // OCCT L1695: EX.Next().
            ex_i += 1;
            // OCCT L1696-1707: add the edges until LastEdge.
            while ex_i < user_edges.len() {
                let e = &user_edges[ex_i];
                if !shape_is_same(e, last_edge) {
                    ww.add_to_wire(&mut pool, wire.clone(), e.clone());
                } else {
                    break;
                }
                ex_i += 1;
            }
            // OCCT L1708-1712.
            let (last_curve, _ff2, _ll2) = match brep_tool_curve(last_edge) {
                Some(hit) => hit,
                None => panic!("null 3D curve (OCCT null-handle deref)"),
            };
            let lf_vert = top_exp_first_vertex(last_edge, true);
            let lf_pnt = brep_tool_pnt(&lf_vert);
            let el = make_edge_c_p_p(&mut pool, &last_curve, lf_pnt, self.my_last_pnt);
            ww.add_to_wire(&mut pool, wire.clone(), el);
        } else {
            // OCCT L1714-1731: only one edge : particular processing.
            let fpar = Self::int_par(&first_curve, self.my_first_pnt);
            let lpar = Self::int_par(&first_curve, self.my_last_pnt);
            let c = if fpar > lpar {
                geom_curve_reversed(&first_curve)
            } else {
                first_curve.clone()
            };

            let ef = make_edge_c_p_p(&mut pool, &c, self.my_first_pnt, self.my_last_pnt);
            ww.add_to_wire(&mut pool, wire.clone(), ef);
        }

        // OCCT L1733-1734: BRepLib_MakeFace f(myPln->Pln(), WW, true).
        let mut f = BRepBuilder::new();
        let fac = make_face_plane_wire(&mut f, &mut pool, my_pln, &wire);

        // OCCT L1736-1744.
        if !brep_algo_is_valid(&fac) {

            profile_ok = false;
            return profile_ok;
        }

        // OCCT L1746-1774.
        if *concavite != 3 {
            // OCCT L1747-1749: concave — face is OK.
            *prof = fac;
        } else {
            // OCCT L1751-1773: not concave — CheckPnt : point slightly
            // inside the material side; BndFace : face/cut of the bounding
            // box in the plane of the profile.
            let u_v = plane_parameters(my_pln, check_pnt);
            let checkpnt2d = glam::DVec2::new(u_v.0, u_v.1);
            if brep_top_adaptor_fclass2d_perform(&fac, checkpnt2d) == rcad_kernel::topods::State::Out
            {
                // OCCT L1762-1767: if face is not the correct part of
                // BndFace take the complementary — BRepAlgoAPI_Cut.
                let c = CutVehicle::new(bnd_face, &fac);
                let c_shape = match c.shape() {
                    Some(s) => s.clone(),
                    None => Shape::null(),
                };
                let w = match explorer_wires_of(&c_shape).into_iter().next() {
                    Some(w) => w,
                    None => panic!("TopExp_Explorer exhausted (OCCT null deref)"),
                };
                let mut ffx = BRepBuilder::new();
                *prof = make_face_wire(&mut ffx, &mut pool, &w);
            } else {
                // OCCT L1769-1772: face is the correct part of BndFace.
                *prof = fac;
            }
        }

        // OCCT L1776-1785.
        if !brep_algo_is_valid(prof) {

            profile_ok = false;
            return profile_ok;
        }
        profile_ok
    }

    // -----------------------------------------------------------------
    // NoSlidingProfile (cxx L1793-2668).
    // -----------------------------------------------------------------

    /// OCCT BRepFeat_RibSlot::NoSlidingProfile (cxx L1793-2668). The
    /// FirstFace/LastFace/FirstVertex/LastVertex parameters are unnamed
    /// (unused) in OCCT. The never-read assignments are the OCCT shape (the
    /// theFV/theLastEdge Nullify sequence and the Pnt/D1 double assignments
    /// of cxx L1894-1897 / L1925-1929 — kept with an allow).
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub(crate) fn no_sliding_profile(
        &mut self,
        prof: &mut Shape,
        revol_rib: bool,
        my_tol: f64,
        concavite: &mut i32,
        my_pln: &Plane,
        bnd: f64,
        bnd_face: &Shape,
        check_pnt: DVec3,
        _first_face: &Shape,
        _last_face: &Shape,
        _first_vertex: &Shape,
        _last_vertex: &Shape,
        first_edge: &Shape,
        last_edge: &Shape,
        on_first_face: bool,
        on_last_face: bool,
    ) -> bool {
        // OCCT L1816.
        let mut profile_ok = true;

        // OCCT L1818-1830.
        let mut the_fv = Shape::null(); // theFV.Nullify()
        let mut the_firstpoint = DVec3::ZERO;
        let mut the_last_edge = Shape::null(); // theLastEdge.Nullify()
        let mut firstpoint;
        let mut lastpoint;
        let mut firstvect;
        let mut lastvect;
        let mut pool = BRep::new();
        let mut bb = BRepBuilder::new();
        let w = bb.make_wire(&mut pool);
        let mut false_first_edge = Shape::null();
        let mut false_last_edge = Shape::null();
        let mut false_only_one = Shape::null();

        // OCCT L1832-1839.
        let (first_curve, f1, _l1) = match brep_tool_curve(first_edge) {
            Some(hit) => hit,
            None => panic!("null 3D curve (OCCT null-handle deref)"),
        };
        let (last_curve, _f2, l2) = match brep_tool_curve(last_edge) {
            Some(hit) => hit,
            None => panic!("null 3D curve (OCCT null-handle deref)"),
        };

        // OCCT L1835-1839.
        firstpoint = first_curve.point_at(f1);
        firstvect = first_curve.derivative_at(f1);
        let mut lastln = Curve3::Line(Line3::new(firstpoint, -firstvect));
        lastpoint = last_curve.point_at(l2);
        lastvect = last_curve.derivative_at(l2);
        let mut firstln = Curve3::Line(Line3::new(lastpoint, lastvect));

        // OCCT L1841: Pt (assigned in the TestOK branch below; the zero init
        // carries the uninitialized read note).
        let mut pt = DVec3::ZERO;

        // OCCT L1843-1846.
        let mut ln2d1 = match geom_api_to_2d(&firstln, my_pln) {
            Some(c) => c,
            None => panic!("GeomAPI::To2d null handle"),
        };
        let mut ln2d2 = match geom_api_to_2d(&lastln, my_pln) {
            Some(c) => c,
            None => panic!("GeomAPI::To2d null handle"),
        };

        let inter = Geom2dAPIInterCurveCurve::new(&ln2d1, &ln2d2, CONFUSION);

        // OCCT L1848-1870.
        let mut test_ok = true;
        if revol_rib {
            let (_o1, d1) = geom_line_parts(&firstln);
            let (_o2, d2) = geom_line_parts(&lastln);
            if gp_dir_is_opposite(d1, d2, my_tol) {
                let par1 = elclib_parameter_lin(&curve3_lin(&firstln), self.my_first_pnt);
                let par2 = elclib_parameter_lin(&curve3_lin(&lastln), self.my_last_pnt);
                if par1 >= my_tol || par2 >= my_tol {
                    *concavite = 2; // parallel and concave
                }
            }
            if gp_dir_is_equal(d1, d2, my_tol) {
                if *concavite == 3 {
                    test_ok = false;
                }
            }
        }

        // OCCT L1872-1884.
        if test_ok && inter.nb_points() > 0 {
            let p2d = inter.point(1);
            pt = my_pln.origin + my_pln.u_dir * p2d.x + my_pln.v_dir * p2d.y;
            let par = Self::int_par(&firstln, pt);
            if par > 0.0 {
                *concavite = 1; // concave
            }
        }

        // --- construction of the face profile (OCCT L1886-2195) ---
        if *concavite == 3 {
            // OCCT L1889-1904.
            if on_first_face {
                false_first_edge = first_edge.clone();
                self.edge_extention(&mut false_first_edge, bnd, true);
                let vv1 = top_exp_first_vertex(&false_first_edge, true);
                firstpoint = brep_tool_pnt(&vv1);
                let (cc, f, _l) = match brep_tool_curve(&false_first_edge) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                firstpoint = cc.point_at(f);
                firstvect = cc.derivative_at(f);
                lastln = Curve3::Line(Line3::new(firstpoint, -firstvect));
                if shape_is_same(first_edge, last_edge) {
                    false_only_one = false_first_edge.clone();
                }
                ln2d2 = match geom_api_to_2d(&lastln, my_pln) {
                    Some(c) => c,
                    None => panic!("GeomAPI::To2d null handle"),
                };
            }
            // OCCT L1905-1932.
            if on_last_face {
                if !shape_is_same(first_edge, last_edge) {
                    false_last_edge = last_edge.clone();
                } else {
                    if false_only_one.is_null() {
                        false_only_one = last_edge.clone();
                    }
                    false_last_edge = false_only_one.clone();
                }
                self.edge_extention(&mut false_last_edge, bnd, false);
                if shape_is_same(first_edge, last_edge) {
                    false_only_one = false_last_edge.clone();
                }
                let vv2 = top_exp_last_vertex(&false_last_edge, true);
                lastpoint = brep_tool_pnt(&vv2);
                let (cc, _f, l) = match brep_tool_curve(&false_last_edge) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                lastpoint = cc.point_at(l);
                lastvect = cc.derivative_at(l);
                lastpoint = brep_tool_pnt(&vv2);
                firstln = Curve3::Line(Line3::new(lastpoint, lastvect));
                ln2d1 = match geom_api_to_2d(&firstln, my_pln) {
                    Some(c) => c,
                    None => panic!("GeomAPI::To2d null handle"),
                };
            }

            // OCCT L1934-1990: BndEdge1/BndEdge2 discovery.
            let mut bnd_edge1 = Shape::null();
            let mut bnd_edge2 = Shape::null();
            let mut bnd_pnt1 = DVec3::ZERO;
            let mut bnd_pnt2 = DVec3::ZERO;
            // OCCT L1935: LastPnt (uninitialized-read note as above).
            let mut last_pnt = DVec3::ZERO;
            let bnd_wire = first_wire_of(bnd_face);
            let bnd_edges = wire_explorer_edges(&bnd_wire);
            for e in &bnd_edges {
                let (c, first, last) = match brep_tool_curve(e) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                let c2d = match geom_api_to_2d(&c, my_pln) {
                    Some(c) => c,
                    None => panic!("GeomAPI::To2d null handle"),
                };
                let intcln1 = Geom2dAPIInterCurveCurve::new(&ln2d1, &c2d, CONFUSION);
                if intcln1.nb_points() > 0 {
                    let p2d = intcln1.point(1);
                    let p = my_pln.origin + my_pln.u_dir * p2d.x + my_pln.v_dir * p2d.y;
                    let parl = Self::int_par(&firstln, p);
                    let parc = Self::int_par(&c, p);
                    if parc >= first && parc <= last && parl >= 0.0 {
                        bnd_edge1 = e.clone();
                        bnd_pnt1 = p;
                    }
                }

                let intcln2 = Geom2dAPIInterCurveCurve::new(&ln2d2, &c2d, CONFUSION);
                if intcln2.nb_points() > 0 {
                    let p2d = intcln2.point(1);
                    let p = my_pln.origin + my_pln.u_dir * p2d.x + my_pln.v_dir * p2d.y;
                    let parl = Self::int_par(&lastln, p);
                    let parc = Self::int_par(&c, p);
                    if parc >= first && parc <= last && parl >= 0.0 {
                        bnd_edge2 = e.clone();
                        bnd_pnt2 = p;
                    }
                }
                // OCCT L1976-1979.
                if !bnd_edge1.is_null() && !bnd_edge2.is_null() {
                    break;
                }
            }

            // OCCT L1982-1990.
            if bnd_edge1.is_null() || bnd_edge2.is_null() {

                profile_ok = false;
                return profile_ok;
            }

            // OCCT L1992-2011: ee1.
            let ee1;
            if the_last_edge.is_null() {
                ee1 = make_edge_p_p(&mut pool, lastpoint, bnd_pnt1);
            } else {
                let v1 = top_exp_last_vertex(&the_last_edge, true);
                let v2 = make_vertex(&mut pool, bnd_pnt1);
                ee1 = make_edge_v_v(&mut pool, &v1, &v2);
            }
            bb.add_to_wire(&mut pool, w.clone(), ee1.clone());
            the_last_edge = ee1.clone();
            if the_fv.is_null() {
                the_fv = top_exp_first_vertex(&ee1, true);
                the_firstpoint = brep_tool_pnt(&the_fv);
            }

            // OCCT L2013-2194.
            if shape_is_same(&bnd_edge1, &bnd_edge2) {
                // OCCT L2015-2053: ee2 + ee3.
                let ee2;
                if the_last_edge.is_null() {
                    ee2 = make_edge_p_p(&mut pool, bnd_pnt1, bnd_pnt2);
                } else {
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let v2 = make_vertex(&mut pool, bnd_pnt2);
                    ee2 = make_edge_v_v(&mut pool, &v1, &v2);
                }
                bb.add_to_wire(&mut pool, w.clone(), ee2.clone());
                the_last_edge = ee2.clone();
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&ee2, true);
                    the_firstpoint = brep_tool_pnt(&the_fv);
                }
                let ee3;
                if the_last_edge.is_null() {
                    ee3 = make_edge_p_p(&mut pool, bnd_pnt2, firstpoint);
                } else {
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let v2 = make_vertex(&mut pool, firstpoint);
                    ee3 = make_edge_v_v(&mut pool, &v1, &v2);
                }
                bb.add_to_wire(&mut pool, w.clone(), ee3.clone());
                the_last_edge = ee3.clone();
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&ee3, true);
                    the_firstpoint = brep_tool_pnt(&the_fv);
                }
            } else {
                // OCCT L2057-2091: walk to BndEdge1.
                let mut idx = 0usize;
                while idx < bnd_edges.len() {
                    let e = &bnd_edges[idx];
                    if shape_is_same(e, &bnd_edge1) {
                        // OCCT L2063-2068.
                        let pp = brep_tool_pnt(&top_exp_last_vertex(e, true));
                        if pp.distance(bnd_pnt1) > brep_tool_tolerance(e) {
                            last_pnt = pp;
                        }
                        // OCCT L2069-2083.
                        let eee;
                        if the_last_edge.is_null() {
                            eee = make_edge_p_p(&mut pool, bnd_pnt1, last_pnt);
                        } else {
                            let v1 = top_exp_last_vertex(&the_last_edge, true);
                            let v2 = make_vertex(&mut pool, last_pnt);
                            eee = make_edge_v_v(&mut pool, &v1, &v2);
                        }
                        bb.add_to_wire(&mut pool, w.clone(), eee.clone());
                        the_last_edge = eee.clone();
                        if the_fv.is_null() {
                            the_fv = top_exp_first_vertex(&eee, true);
                            the_firstpoint = brep_tool_pnt(&the_fv);
                        }
                        break;
                    }
                    idx += 1;
                }

                // OCCT L2093-2104: the wrapped advance.
                if idx + 1 < bnd_edges.len() {
                    idx += 1;
                } else {
                    idx = 0;
                }

                // OCCT L2105-2172: the Fin loop.
                let mut fin = false;
                while !fin {
                    let e = &bnd_edges[idx];
                    if !shape_is_same(e, &bnd_edge2) {
                        // OCCT L2111-2133.
                        let pp = brep_tool_pnt(&top_exp_last_vertex(e, true));
                        let eee1;
                        if the_last_edge.is_null() {
                            eee1 = make_edge_p_p(&mut pool, last_pnt, pp);
                        } else {
                            let v1 = top_exp_last_vertex(&the_last_edge, true);
                            let v2 = make_vertex(&mut pool, pp);
                            eee1 = make_edge_v_v(&mut pool, &v1, &v2);
                        }
                        bb.add_to_wire(&mut pool, w.clone(), eee1.clone());
                        the_last_edge = eee1.clone();
                        if the_fv.is_null() {
                            the_fv = top_exp_first_vertex(&eee1, true);
                            the_firstpoint = brep_tool_pnt(&the_fv);
                        }
                        last_pnt = pp;
                    } else {
                        // OCCT L2135-2158.
                        fin = true;
                        let eee2;
                        if the_last_edge.is_null() {
                            eee2 = make_edge_p_p(&mut pool, last_pnt, bnd_pnt2);
                        } else {
                            let v1 = top_exp_last_vertex(&the_last_edge, true);
                            let v2 = make_vertex(&mut pool, bnd_pnt2);
                            eee2 = make_edge_v_v(&mut pool, &v1, &v2);
                        }
                        bb.add_to_wire(&mut pool, w.clone(), eee2.clone());
                        the_last_edge = eee2.clone();
                        if the_fv.is_null() {
                            the_fv = top_exp_first_vertex(&eee2, true);
                            the_firstpoint = brep_tool_pnt(&the_fv);
                        }
                        last_pnt = bnd_pnt2;
                    }
                    // OCCT L2160-2171: the wrapped advance.
                    if idx + 1 < bnd_edges.len() {
                        idx += 1;
                    } else {
                        idx = 0;
                    }
                }

                // OCCT L2174-2193: ee3.
                let eee3;
                if the_last_edge.is_null() {
                    eee3 = make_edge_p_p(&mut pool, bnd_pnt2, firstpoint);
                } else {
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let v2 = make_vertex(&mut pool, firstpoint);
                    eee3 = make_edge_v_v(&mut pool, &v1, &v2);
                }
                bb.add_to_wire(&mut pool, w.clone(), eee3.clone());
                the_last_edge = eee3.clone();
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&eee3, true);
                    the_firstpoint = brep_tool_pnt(&the_fv);
                }
            }
        }

        // OCCT L2197-2219: Concavite == 1 (the first closing edge).
        if *concavite == 1 {
            let eee4;
            if the_last_edge.is_null() {
                eee4 = make_edge_p_p(&mut pool, pt, firstpoint);
            } else {
                let v1 = top_exp_last_vertex(&the_last_edge, true);
                let v2 = make_vertex(&mut pool, firstpoint);
                eee4 = make_edge_v_v(&mut pool, &v1, &v2);
            }
            bb.add_to_wire(&mut pool, w.clone(), eee4.clone());
            if the_fv.is_null() {
                the_fv = top_exp_first_vertex(&eee4, true);
                the_firstpoint = brep_tool_pnt(&the_fv);
            }
            the_last_edge = eee4.clone();
        }

        // OCCT L2221-2580: the user-wire edges.
        if shape_is_same(first_edge, last_edge) {
            // OCCT L2223-2227.
            if !data_map_is_bound(&self.my_lfmap, first_edge) {
                data_map_bind(&mut self.my_lfmap, first_edge);
            }
            // OCCT L2228-2318.
            if on_first_face || on_last_face {
                // OCCT L2230-2271.
                let (cc, _f, _l) = match brep_tool_curve(&false_only_one) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                let the_edge;
                if !the_last_edge.is_null() {
                    // OCCT L2235-2247.
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let pp = brep_tool_pnt(&top_exp_last_vertex(&false_only_one, true));
                    let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                        the_fv.clone()
                    } else {
                        top_exp_last_vertex(&false_only_one, true)
                    };
                    the_edge = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                } else {
                    // OCCT L2251-2263.
                    let v1 = top_exp_first_vertex(&false_only_one, true);
                    let pp = brep_tool_pnt(&top_exp_last_vertex(&false_only_one, true));
                    let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                        the_fv.clone()
                    } else {
                        top_exp_last_vertex(&false_only_one, true)
                    };
                    the_edge = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                }
                // OCCT L2265-2271.
                data_map_change_find(&mut self.my_lfmap, shape_key(first_edge))
                    .push(the_edge.clone());
                bb.add_to_wire(&mut pool, w.clone(), the_edge.clone());
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&the_edge, true);
                }
                the_last_edge = the_edge.clone();
            } else {
                // OCCT L2273-2318.
                let (cc, _f, _l) = match brep_tool_curve(first_edge) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                let the_edge;
                if !the_last_edge.is_null() {
                    // OCCT L2280-2293.
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    // Attention case Wire Reversed -> LastVertex without true
                    let pp = brep_tool_pnt(&top_exp_last_vertex(first_edge, false));
                    let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                        the_fv.clone()
                    } else {
                        top_exp_last_vertex(first_edge, false)
                    };
                    the_edge = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                } else {
                    // OCCT L2297-2309.
                    let v1 = top_exp_first_vertex(first_edge, true);
                    let pp = brep_tool_pnt(&top_exp_last_vertex(first_edge, true));
                    let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                        the_fv.clone()
                    } else {
                        top_exp_last_vertex(first_edge, true)
                    };
                    the_edge = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                }
                // OCCT L2311-2317.
                data_map_change_find(&mut self.my_lfmap, shape_key(first_edge))
                    .push(the_edge.clone());
                bb.add_to_wire(&mut pool, w.clone(), the_edge.clone());
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&the_edge, true);
                }
                the_last_edge = the_edge.clone();
            }
        } else {
            // OCCT L2320-2580.
            // OCCT L2322-2326.
            if !data_map_is_bound(&self.my_lfmap, first_edge) {
                data_map_bind(&mut self.my_lfmap, first_edge);
            }
            // OCCT L2327-2350.
            if !on_first_face {
                let (cc, _f, _l) = match brep_tool_curve(first_edge) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                let the_edge;
                if !the_last_edge.is_null() {
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let v2 = top_exp_last_vertex(first_edge, true);
                    the_edge = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                } else {
                    the_edge = first_edge.clone();
                }
                data_map_change_find(&mut self.my_lfmap, shape_key(first_edge))
                    .push(the_edge.clone());
                bb.add_to_wire(&mut pool, w.clone(), the_edge.clone());
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&the_edge, true);
                }
                the_last_edge = the_edge.clone();
            } else {
                // OCCT L2351-2374.
                let (cc, _f, _l) = match brep_tool_curve(&false_first_edge) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                let the_edge;
                if !the_last_edge.is_null() {
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let v2 = top_exp_last_vertex(&false_first_edge, true);
                    the_edge = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                } else {
                    the_edge = false_first_edge.clone();
                }
                data_map_change_find(&mut self.my_lfmap, shape_key(first_edge))
                    .push(the_edge.clone());
                bb.add_to_wire(&mut pool, w.clone(), the_edge.clone());
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&the_edge, true);
                }
                the_last_edge = the_edge.clone();
            }

            // OCCT L2376-2424: the middle-wire walk.
            let user_edges = wire_explorer_edges(&self.my_wire);
            let mut ex_i = 0usize;
            while ex_i < user_edges.len() {
                let e = &user_edges[ex_i];
                if shape_is_same(e, first_edge) {
                    break;
                }
                ex_i += 1;
            }

            // OCCT L2386: ex.Next().
            ex_i += 1;

            // OCCT L2388-2424.
            while ex_i < user_edges.len() {
                let e = &user_edges[ex_i];
                if !shape_is_same(e, last_edge) {
                    // OCCT L2393-2397.
                    if !data_map_is_bound(&self.my_lfmap, e) {
                        data_map_bind(&mut self.my_lfmap, e);
                    }
                    let (cc, _f, _l) = match brep_tool_curve(e) {
                        Some(hit) => hit,
                        None => panic!("null 3D curve (OCCT null-handle deref)"),
                    };
                    let eee;
                    if !the_last_edge.is_null() {
                        let v1 = top_exp_last_vertex(&the_last_edge, true);
                        let v2 = top_exp_last_vertex(e, true);
                        eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                    } else {
                        eee = e.clone();
                    }
                    data_map_change_find(&mut self.my_lfmap, shape_key(e))
                        .push(eee.clone());
                    bb.add_to_wire(&mut pool, w.clone(), eee.clone());
                    if the_fv.is_null() {
                        the_fv = top_exp_first_vertex(&eee, true);
                    }
                    the_last_edge = eee.clone();
                } else {
                    break;
                }
                ex_i += 1;
            }

            // OCCT L2426-2529.
            if !on_last_face {
                if !shape_is_same(first_edge, last_edge) {
                    // OCCT L2430-2478.
                    let edg = user_edges
                        .get(ex_i)
                        .cloned()
                        .expect("WireExplorer exhausted (OCCT null deref)");
                    if !data_map_is_bound(&self.my_lfmap, &edg) {
                        data_map_bind(&mut self.my_lfmap, &edg);
                    }
                    let (cc, _f, _l) = match brep_tool_curve(&edg) {
                        Some(hit) => hit,
                        None => panic!("null 3D curve (OCCT null-handle deref)"),
                    };
                    let eee;
                    if !the_last_edge.is_null() {
                        let v1 = top_exp_last_vertex(&the_last_edge, true);
                        let pp = brep_tool_pnt(&top_exp_last_vertex(&edg, true));
                        let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                            the_fv.clone()
                        } else {
                            top_exp_last_vertex(&edg, true)
                        };
                        eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                    } else {
                        let v1 = top_exp_first_vertex(&edg, true);
                        let pp = brep_tool_pnt(&top_exp_last_vertex(&edg, true));
                        let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                            the_fv.clone()
                        } else {
                            top_exp_last_vertex(&edg, true)
                        };
                        eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                    }
                    data_map_change_find(&mut self.my_lfmap, shape_key(&edg))
                        .push(eee.clone());
                    bb.add_to_wire(&mut pool, w.clone(), eee.clone());
                    if the_fv.is_null() {
                        the_fv = top_exp_first_vertex(&eee, true);
                    }
                    the_last_edge = eee.clone();
                } else {
                    // OCCT L2481-2528.
                    let eee;
                    if !data_map_is_bound(&self.my_lfmap, last_edge) {
                        data_map_bind(&mut self.my_lfmap, last_edge);
                    }
                    let (cc, _f, _l) = match brep_tool_curve(&false_only_one) {
                        Some(hit) => hit,
                        None => panic!("null 3D curve (OCCT null-handle deref)"),
                    };
                    if !the_last_edge.is_null() {
                        let v1 = top_exp_last_vertex(&the_last_edge, true);
                        let pp = brep_tool_pnt(&top_exp_last_vertex(&false_only_one, true));
                        let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                            the_fv.clone()
                        } else {
                            top_exp_last_vertex(&false_only_one, true)
                        };
                        eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                    } else {
                        let v1 = top_exp_first_vertex(&false_only_one, true);
                        let pp = brep_tool_pnt(&top_exp_last_vertex(&false_only_one, true));
                        let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                            the_fv.clone()
                        } else {
                            top_exp_last_vertex(&false_only_one, true)
                        };
                        eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                    }
                    data_map_change_find(&mut self.my_lfmap, shape_key(last_edge))
                        .push(eee.clone());
                    bb.add_to_wire(&mut pool, w.clone(), eee.clone());
                    if the_fv.is_null() {
                        the_fv = top_exp_first_vertex(&eee, true);
                    }
                    the_last_edge = eee.clone();
                }
            } else {
                // OCCT L2530-2579.
                let eee;
                if !data_map_is_bound(&self.my_lfmap, last_edge) {
                    data_map_bind(&mut self.my_lfmap, last_edge);
                }
                let (cc, _f, _l) = match brep_tool_curve(&false_last_edge) {
                    Some(hit) => hit,
                    None => panic!("null 3D curve (OCCT null-handle deref)"),
                };
                if !the_last_edge.is_null() {
                    let v1 = top_exp_last_vertex(&the_last_edge, true);
                    let pp = brep_tool_pnt(&top_exp_last_vertex(&false_last_edge, true));
                    let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                        the_fv.clone()
                    } else {
                        top_exp_last_vertex(&false_last_edge, true)
                    };
                    eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                } else {
                    let v1 = top_exp_first_vertex(&false_last_edge, true);
                    let pp = brep_tool_pnt(&top_exp_last_vertex(&false_last_edge, true));
                    let v2 = if !the_fv.is_null() && the_firstpoint.distance(pp) <= my_tol {
                        the_fv.clone()
                    } else {
                        top_exp_last_vertex(&false_last_edge, true)
                    };
                    eee = make_edge_c_v_v(&mut pool, &cc, &v1, &v2);
                }
                data_map_change_find(&mut self.my_lfmap, shape_key(last_edge))
                    .push(eee.clone());
                bb.add_to_wire(&mut pool, w.clone(), eee.clone());
                if the_fv.is_null() {
                    the_fv = top_exp_first_vertex(&eee, true);
                }
                the_last_edge = eee.clone();
            }
        }

        // OCCT L2582-2609: Concavite == 1 (the closing edge).
        if *concavite == 1 {
            let eef;
            if the_last_edge.is_null() {
                eef = make_edge_p_p(&mut pool, lastpoint, pt);
            } else {
                // OCCT L2592-2601.
                let v1 = top_exp_last_vertex(&the_last_edge, true);
                let vv = make_vertex(&mut pool, pt);
                let mut v2 = vv.clone();
                if !the_fv.is_null() && pt.distance(the_firstpoint) <= my_tol {
                    v2 = the_fv.clone();
                }
                eef = make_edge_v_v(&mut pool, &v1, &v2);
            }
            bb.add_to_wire(&mut pool, w.clone(), eef.clone());
            if the_fv.is_null() {
                the_fv = top_exp_first_vertex(&eef, true);
            }
            the_last_edge = eef.clone();
        }

        // OCCT L2611-2616: Concavite == 2.
        if *concavite == 2 {
            let ee = make_edge_p_p(&mut pool, lastpoint, firstpoint);
            bb.add_to_wire(&mut pool, w.clone(), ee);
        }

        // OCCT L2618-2620.
        let closed = brep_tool_is_closed(&w);
        wire_set_closed(&mut pool, &w, closed);
        let mut fa = BRepBuilder::new();
        let fac = make_face_plane_wire(&mut fa, &mut pool, my_pln, &w);

        // OCCT L2622-2630.
        if !brep_algo_is_valid(&fac) {

            profile_ok = false;
            return profile_ok;
        }

        // OCCT L2632-2656.
        if *concavite == 3 {
            let u_v = plane_parameters(my_pln, check_pnt);
            let checkpnt2d = glam::DVec2::new(u_v.0, u_v.1);
            if brep_top_adaptor_fclass2d_perform(&fac, checkpnt2d) == rcad_kernel::topods::State::Out
            {
                // OCCT L2641-2646: BRepAlgoAPI_Cut c(BndFace, fac);
                // UpdateDescendants(c, c.Shape(), false).
                let c = CutVehicle::new(bnd_face, &fac);
                let c_shape = match c.shape() {
                    Some(s) => s.clone(),
                    None => Shape::null(),
                };
                self.update_descendants_bop(&c, &c_shape, false);
                let ww_shape = match explorer_wires_of(&c_shape).into_iter().next() {
                    Some(w) => w,
                    None => panic!("TopExp_Explorer exhausted (OCCT null deref)"),
                };
                let mut ff = BRepBuilder::new();
                *prof = make_face_wire(&mut ff, &mut pool, &ww_shape);
            } else {
                *prof = fac;
            }
        } else {
            *prof = fac;
        }

        // OCCT L2658-2667.
        if !brep_algo_is_valid(prof) {

            profile_ok = false;
            return profile_ok;
        }
        profile_ok
    }
}

// ---------------------------------------------------------------------------
// Local exploration / geometry helpers (the TopExp_Explorer re-host readers).
// ---------------------------------------------------------------------------

/// OCCT TopExp_Explorer(S, TopAbs_FACE) — the faces of S (the crate explorer
/// re-host, brep_feat_builder::explorer).
fn explorer_edges_of_faces(the_s: &Shape) -> Vec<Shape> {
    crate::feat::brep_feat_builder::explorer(the_s, ShapeType::Face, ShapeType::Shape)
}

/// OCCT TopExp_Explorer(S, TopAbs_EDGE).
fn explorer_edges_of(the_s: &Shape) -> Vec<Shape> {
    crate::feat::brep_feat_builder::explorer(the_s, ShapeType::Edge, ShapeType::Shape)
}

/// OCCT TopExp_Explorer(S, TopAbs_VERTEX).
fn explorer_vertices_of(the_s: &Shape) -> Vec<Shape> {
    crate::feat::brep_feat_builder::explorer(the_s, ShapeType::Vertex, ShapeType::Shape)
}

/// OCCT TopExp_Explorer(S, TopAbs_WIRE).
fn explorer_wires_of(the_s: &Shape) -> Vec<Shape> {
    crate::feat::brep_feat_builder::explorer(the_s, ShapeType::Wire, ShapeType::Shape)
}

/// OCCT TopExp_Explorer(BndFace, TopAbs_WIRE).Current() — the first wire
/// (the exhausted-explorer continuation is not defined in OCCT: TopoDS::Wire
/// on the null current).
fn first_wire_of(the_face: &Shape) -> Shape {
    explorer_wires_of(the_face)
        .into_iter()
        .next()
        .expect("BndFace first wire (OCCT null-handle deref)")
}

/// OCCT Geom_Line carrier extraction.
fn curve3_lin(the_c: &Curve3) -> Line3 {
    match the_c {
        Curve3::Line(l) => *l,
        _ => panic!("Geom_Line carrier expected"),
    }
}

/// OCCT ElSLib::Parameters(myPln->Pln(), P, u, v) — the plane parameters
/// (the ElSLib.cxx plane branch).
fn plane_parameters(the_pln: &Plane, the_p: DVec3) -> (f64, f64) {
    rcad_kernel::math::el::elslib_plane_parameters(
        the_p, the_pln.origin, the_pln.u_dir, the_pln.v_dir,
    )
}
