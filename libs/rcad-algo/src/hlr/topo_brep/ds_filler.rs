// OCCT HLRTopoBRep_DSFiller (TKHLR/HLRTopoBRep/HLRTopoBRep_DSFiller.hxx
// L1-79 + .cxx L1-758) — fills a HLRTopoBRep_Data with the outlines of a
// shape (HLRTopoBRep_DSFiller::Insert / InsertFace / MakeVertex /
// InsertVertex / ProcessEdges).
//
// Translation notes:
// - OCCT `static` class methods map to free functions in this module;
//   Insert is the public entry, the rest are the OCCT private members
//   (pub(crate) for the sibling-module consumers).
// - The OCCT global TShape graph travels as an explicit `&mut BRep` — every
//   `BRep_Builder` mutation (MakeEdge / Add / Range / UpdateEdge /
//   UpdateVertex / Remove) edits that BRep in place (the edge_mut_inplace
//   identity-preserving contract), so the DS Shape handles observe the same
//   TShape mutations as the OCCT handles do.
// - The MST map is the OCCT
//   `NCollection_DataMap<TopoDS_Shape, BRepTopAdaptor_Tool,
//   TopTools_ShapeMapHasher>` encoded as `Vec<(Shape, BRepTopAdaptorTool)>`
//   keyed by `Shape::ptr_id()` (the data.rs map precedent).
// - The OCCT `TopoDS_Shape` handle copies in the DS lists share the TShape
//   with the BRep; the rcad `Shape` clones do the same, so the IntL merge
//   loop mutations performed on local copies are visible through the lists.

use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3, Curve2d, Curve3};
use rcad_kernel::p_confusion;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, Shape, TShape};

use crate::geomalgo::approx_int::ApproxParamType;
use crate::geomalgo::brep_approx::ApproxLine;
use crate::geomalgo::brep_approx_approx::BRepApproxApprox;
use crate::hlr::contap::contour::Contour;
use crate::hlr::contap::i_type::IType;
use crate::hlr::contap::point::Point as ContapPoint;
use crate::topalgo::brep_adaptor::curve2d::BRepCurve2d;
use crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool;

use super::data::Data;

/// OCCT NCollection_DataMap::IsBound over the MST map (the
/// TopTools_ShapeMapHasher identity — the data.rs ptr_id precedent).
fn mst_is_bound(mst: &[(Shape, BRepTopAdaptorTool)], s: &Shape) -> bool {
    let id = s.ptr_id();
    mst.iter().any(|(k, _)| k.ptr_id() == id)
}

/// The MST map position of the shape (OCCT ChangeFind / Find).
fn mst_pos(mst: &[(Shape, BRepTopAdaptorTool)], s: &Shape) -> usize {
    let id = s.ptr_id();
    mst.iter()
        .position(|(k, _)| k.ptr_id() == id)
        .expect("Standard_NoSuchObject: MST key not bound")
}

/// OCCT MST.Bind(S1, BRT) — insert the tool under the shape key.
fn mst_bind(mst: &mut Vec<(Shape, BRepTopAdaptorTool)>, s: &Shape, tool: BRepTopAdaptorTool) {
    mst.push((s.clone(), tool));
}

/// OCCT TopExp_Explorer(S, TopAbs_FACE) — the recursive face enumeration.
/// The rcad walk follows the TShape child lists; the cumulative
/// Location/Orientation composition of the OCCT explorer is the identity
/// pass-through here (the same convention as the wire-edge enumeration of
/// topol_tool_brep::initialize_surface).
fn top_exp_faces(s: &Shape, out: &mut Vec<Shape>) {
    match &*s.data {
        TShape::Compound(children) => {
            for c in children {
                top_exp_faces(c, out);
            }
        }
        TShape::CompSolid(cs) => {
            for c in cs {
                top_exp_faces(c, out);
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                top_exp_faces(sh, out);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                top_exp_faces(f, out);
            }
        }
        TShape::Face(_) => out.push(s.clone()),
        // Wire / Edge / Vertex: the explorer finds nothing below them.
        _ => {}
    }
}

/// The flat knot vector of the rcad kernel BSpline representation from the
/// OCCT compressed (knots, multiplicities) arrays (each knot repeated its
/// multiplicity) — the bop/int_tools/face_make_curve.rs expand_knots
/// precedent.
fn expand_knots(knots: &[f64], mults: &[usize]) -> Vec<f64> {
    let mut out = Vec::new();
    for (i, &k) in knots.iter().enumerate() {
        let m = if i < mults.len() { mults[i] } else { 1 };
        for _ in 0..m {
            out.push(k);
        }
    }
    out
}

//=======================================================================
// function : Insert
// purpose  : explore the faces and insert them
// OCCT HLRTopoBRep_DSFiller.cxx L62-113
//=======================================================================

/// OCCT HLRTopoBRep_DSFiller::Insert(S, FO, DS, MST, nbIso).
pub fn insert(
    brep: &mut BRep,
    s: &Shape,
    fo: &mut Contour,
    ds: &mut Data,
    mst: &mut Vec<(Shape, BRepTopAdaptorTool)>,
    nb_iso: usize,
) {
    // NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> ShapeMap — the
    // visited-face set (the ptr_id encoding of the data.rs map precedent).
    let mut shape_map: Vec<u64> = Vec::new();
    // TopExp_Explorer ex(S, TopAbs_FACE).
    let mut faces = Vec::new();
    top_exp_faces(s, &mut faces);
    ds.clear();
    let with_pcurve = true; // instead of nbIso != 0;
    let mut f: usize = 0;

    // The BRep handle for the MST tools: the OCCT tools read the global
    // TShape graph; the rcad tool holds a refcounted clone of the BRep (the
    // brep_adaptor::curve2d ownership note) — all pre-existing face data is
    // shared through the TShape Arcs.
    let tool_brep = Arc::new(brep.clone());

    for ex_current in faces {
        // if (ShapeMap.Add(ex.Current())).
        let id = ex_current.ptr_id();
        if !shape_map.contains(&id) {
            shape_map.push(id);
            f += 1;
            let mut s1 = ex_current.clone();
            s1.orientation = Orientation::Forward;
            // occ::handle<BRepTopAdaptor_TopolTool> Domain;
            // occ::handle<Adaptor3d_Surface> Surface; — resolved from the
            // bound MST tool below.
            let surface;
            let domain_pos: usize;
            if mst_is_bound(mst, &s1) {
                // BRepTopAdaptor_Tool& BRT = MST.ChangeFind(S1);
                let pos = mst_pos(mst, &s1);
                let brt = &mut mst[pos].1;
                surface = brt
                    .get_surface()
                    .expect("BRepTopAdaptor_Tool::GetSurface")
                    .clone();
                domain_pos = pos;
            } else {
                // BRepTopAdaptor_Tool BRT(S1, Precision::PConfusion());
                let brt =
                    BRepTopAdaptorTool::new_face(tool_brep.clone(), &s1, p_confusion());
                surface = brt
                    .get_surface()
                    .expect("BRepTopAdaptor_Tool::GetSurface")
                    .clone();
                // MST.Bind(S1, BRT).
                mst_bind(mst, &s1, brt);
                domain_pos = mst_pos(mst, &s1);
            }
            // FO.Perform(Surface, Domain) — Domain = BRT.GetTopolTool().
            fo.perform(&surface, mst[domain_pos].1.get_topol_tool());
            if fo.is_done() {
                if !fo.is_empty() {
                    insert_face(brep, f, &s1, fo, ds, with_pcurve);
                }
            }
            if nb_iso != 0 {
                super::face_iso_liner::perform(f, &s1, ds, nb_iso);
            }
        }
        // ex.Next().
    }
    process_edges(brep, ds);
}

//=======================================================================
// function : InsertFace
// purpose  : private, insert the outlines of a face
// OCCT HLRTopoBRep_DSFiller.cxx L120-592
//=======================================================================

/// OCCT HLRTopoBRep_DSFiller::InsertFace(FI, F, FO, DS, withPCurve).
pub(crate) fn insert_face(
    brep: &mut BRep,
    _fi: usize,
    f: &Shape,
    fo: &mut Contour,
    ds: &mut Data,
    with_pcurve: bool,
) {
    let mut b = BRepBuilder::new();

    // Insert the intersections of FO in DS

    // const double tol = BRep_Tool::Tolerance(F);
    let tol = brep.tolerance(f);
    // NCollection_List<TopoDS_Shape>& IntL = DS.AddIntL(F);
    // NCollection_List<TopoDS_Shape>& OutL = DS.AddOutL(F);
    // (the rcad lists are re-fetched per access — AddIntL/AddOutL bind the
    // FaceData once and the fetches address the same list).
    ds.add_int_l(f);
    ds.add_out_l(f);

    // TopoDS_Vertex VF, VL; — the loop-carried vertices (null until set).
    let mut vf = Shape::null();
    let mut vl = Shape::null();
    // The commented-out VM / closed-edge explorer block (cxx L133-145) is
    // not translated (OCCT keeps it commented out).

    let nb_lines = fo.nb_lines();
    for cur_line in 1..=nb_lines {
        let line = fo.line(cur_line);
        let nb_points = line.nb_vertex();
        if line.type_contour() == IType::Restriction {
            // OutLine on restriction
            // TopoDS_Edge E = (*(BRepAdaptor_Curve2d*)(Line.Arc().get())).Edge().
            let e = line
                .arc()
                .as_any()
                .downcast_ref::<BRepCurve2d>()
                .expect("Standard_TypeMismatch: BRepAdaptor_Curve2d")
                .edge()
                .clone();
            ds.add_out_l(f).push(e.clone());
            // TopExp::Vertices(E, VF, VL) — the first and last stored vertex.
            vf = brep.first_vertex(&e);
            vl = brep.last_vertex(&e);
            // insert the Internal points.

            for cur_point in 1..=nb_points {
                let p = line.vertex(cur_point).clone();
                if p.is_internal() {
                    // P.Value().IsEqual(BRep_Tool::Pnt(VF),
                    // BRep_Tool::Tolerance(VF)).
                    if p.value().distance(brep.vertex_position(&vf))
                        <= brep.vertex_tolerance(&vf)
                    {
                        if p.value().distance(brep.vertex_position(&vl))
                            <= brep.vertex_tolerance(&vl)
                        {
                            insert_vertex(brep, &p, tol, &e, ds);
                        }
                    }
                }
            }
        } else {
            for cur_point in 1..=nb_points {
                // const Contap_Point PF = Line.Vertex(CurPoint);
                let pf = line.vertex(cur_point).clone();
                if pf.is_internal() && cur_point != 1 {
                    vf = vl.clone();
                } else {
                    vf = make_vertex(brep, &pf, tol, ds);
                }
                let par_f = pf.parameter_on_line();

                if cur_point < nb_points {
                    // const Contap_Point PL = Line.Vertex(CurPoint + 1);
                    let pl = line.vertex(cur_point + 1).clone();
                    vl = make_vertex(brep, &pl, tol, ds);
                    let par_l = pl.parameter_on_line();

                    if (par_l - par_f) > p_confusion() {
                        // occ::handle<Geom_Curve> C;
                        // occ::handle<Geom2d_Curve> C2d;
                        let mut c: Option<Curve3> = None;
                        let mut c2d: Option<Curve2d> = None;
                        let mut first = par_f;
                        let mut last = par_l;
                        let mut insuffisant_number_of_points = false;

                        match line.type_contour() {
                            IType::Lin => {
                                // C = new Geom_Line(Line.Line()).
                                let c_lin = Curve3::Line(line.line());
                                c = Some(c_lin);
                                if with_pcurve {
                                    // occ::handle<Geom_Surface> S =
                                    //   BRep_Tool::Surface(F);
                                    let s = brep
                                        .face_surface(f)
                                        .expect("BRep_Tool::Surface")
                                        .clone();
                                    let _tol = 1e-7; // double Tol = 1e-7; (the rcad curve2d_simple carries TOL_DEFAULT)
                                    // C2d = GeomProjLib::Curve2d(C, first,
                                    // last, S, Tol) — the rcad kernel
                                    // carries the 1e-7 tolerance as its
                                    // TOL_DEFAULT and the surface natural
                                    // domain (curve2d_simple).
                                    c2d = rcad_kernel::base::geom_proj_lib::curve2d_simple(
                                        c.as_ref().unwrap(),
                                        first,
                                        last,
                                        &s,
                                    );
                                }
                            }
                            IType::Circle => {
                                // C = new Geom_Circle(Line.Circle()).
                                let c_circ = Curve3::Circle(line.circle());
                                c = Some(c_circ);
                                if with_pcurve {
                                    // TopLoc_Location Loc;
                                    // occ::handle<Geom_Surface> S =
                                    //   BRep_Tool::Surface(F, Loc);
                                    // if (!Loc.IsIdentity())
                                    //   S = S->Transformed(
                                    //     Loc.Transformation()); — the rcad
                                    // world-surface accessor applies the
                                    // location.
                                    let s = brep
                                        .face_surface_world(f)
                                        .expect("BRep_Tool::Surface")
                                        .clone();
                                    let _tol = 1e-7; // double Tol = 1e-7; (the rcad curve2d_simple carries TOL_DEFAULT)
                                    // C2d = GeomProjLib::Curve2d(C, first,
                                    // last, S, Tol).
                                    c2d = rcad_kernel::base::geom_proj_lib::curve2d_simple(
                                        c.as_ref().unwrap(),
                                        first,
                                        last,
                                        &s,
                                    );
                                }
                            }
                            IType::Walking => {
                                // copy the points
                                // int ipF = int(parF); int ipL = int(parL).
                                let ip_f = par_f as i32 as usize;
                                let ip_l = par_l as i32 as usize;

                                if ip_l < ip_f + 1 {
                                    // (ipL - ipF < 1)
                                    insuffisant_number_of_points = true;
                                    // std::cout<<"\n !! Pb ds
                                    // HLRTopoBRep_DSFiller.cxx (Contour App
                                    // Nbp <3)"<<std::endl; — commented out.
                                }
                                // The commented-out tangent-interpolation
                                // block `else if(ipL-ipF < 6)` (cxx L249-320)
                                // is not translated (OCCT keeps it commented
                                // out).
                                else if ip_l < ip_f + 5 {
                                    // (ipL - ipF < 5)
                                    // const int nbp = ipL - ipF + 1;
                                    let nbp = ip_l - ip_f + 1;
                                    // NCollection_Array1<double> knots(1, nbp);
                                    // NCollection_Array1<int> mults(1, nbp);
                                    // NCollection_Array1<gp_Pnt> Points(1, nbp);
                                    let mut knots: Vec<f64> = Vec::with_capacity(nbp);
                                    let mut mults: Vec<usize> = Vec::with_capacity(nbp);
                                    let mut points: Vec<DVec3> = Vec::with_capacity(nbp);

                                    for i in 1..=nbp {
                                        knots.push(i as f64); // knots.SetValue(i, (double)i)
                                        mults.push(1); // mults.SetValue(i, 1)
                                        // Points.SetValue(i,
                                        //   Line.Point(i + ipF - 1).Value()).
                                        points.push(line.point(i + ip_f - 1).value());
                                    }
                                    // mults(1) = mults(nbp) = 2.
                                    mults[0] = 2;
                                    mults[nbp - 1] = 2;
                                    // C = new Geom_BSplineCurve(Points, knots,
                                    // mults, 1) — the rcad kernel BSpline
                                    // carries the flat knot vector.
                                    c = Some(Curve3::BSpline(BSplineCurve3 {
                                        degree: 1,
                                        knots: expand_knots(&knots, &mults),
                                        control_points: points,
                                        weights: vec![],
                                        is_periodic: false,
                                    }));

                                    if with_pcurve {
                                        // NCollection_Array1<gp_Pnt2d>
                                        //   Points2d(1, nbp);
                                        let mut points2d: Vec<DVec2> =
                                            Vec::with_capacity(nbp);
                                        for i in 1..=nbp {
                                            // double u, v;
                                            // Line.Point(i + ipF - 1)
                                            //   .ParametersOnS2(u, v);
                                            let (u, v) = line
                                                .point(i + ip_f - 1)
                                                .parameters_on_surface(false);
                                            // Points2d.SetValue(i,
                                            //   gp_Pnt2d(u, v)).
                                            points2d.push(DVec2::new(u, v));
                                        }
                                        // C2d = new Geom2d_BSplineCurve(
                                        //   Points2d, knots, mults, 1).
                                        c2d = Some(Curve2d::BSpline(BSplineCurve2 {
                                            degree: 1,
                                            knots: {
                                                let kn = knots.clone();
                                                let mu = mults.clone();
                                                expand_knots(&kn, &mu)
                                            },
                                            control_points: points2d,
                                            weights: vec![],
                                        }));
                                    }
                                    first = 1.0;
                                    last = nbp as f64;
                                } else {
                                    // const int nbp = ipL - ipF + 1;
                                    let nbp = ip_l - ip_f + 1;
                                    let mut knots: Vec<f64> = Vec::with_capacity(nbp);
                                    let mut mults: Vec<usize> = Vec::with_capacity(nbp);
                                    let mut points: Vec<DVec3> = Vec::with_capacity(nbp);

                                    // double Maxx, Maxy, Maxz, Maxu, Maxv;
                                    // double Minx, Miny, Minz, Minu, Minv;
                                    let mut max_x: f64;
                                    let mut max_y: f64;
                                    let mut max_z: f64;
                                    let mut max_u: f64;
                                    let mut max_v: f64;
                                    let mut min_x: f64;
                                    let mut min_y: f64;
                                    let mut min_z: f64;
                                    let mut min_u: f64;
                                    let mut min_v: f64;
                                    // Maxx = Maxy = Maxz = Maxu = Maxv =
                                    //   -RealLast();
                                    max_x = -f64::MAX;
                                    max_y = -f64::MAX;
                                    max_z = -f64::MAX;
                                    max_u = -f64::MAX;
                                    max_v = -f64::MAX;
                                    // Minx = Miny = Minz = Minu = Minv =
                                    //   RealLast();
                                    min_x = f64::MAX;
                                    min_y = f64::MAX;
                                    min_z = f64::MAX;
                                    min_u = f64::MAX;
                                    min_v = f64::MAX;

                                    for i in 1..=nbp {
                                        knots.push(i as f64);
                                        mults.push(1);
                                        // const gp_Pnt& P =
                                        //   Line.Point(i + ipF - 1).Value();
                                        let p = line.point(i + ip_f - 1).value();
                                        if p.x < min_x {
                                            min_x = p.x;
                                        }
                                        if p.y < min_y {
                                            min_y = p.y;
                                        }
                                        if p.z < min_z {
                                            min_z = p.z;
                                        }
                                        if p.x > max_x {
                                            max_x = p.x;
                                        }
                                        if p.y > max_y {
                                            max_y = p.y;
                                        }
                                        if p.z > max_z {
                                            max_z = p.z;
                                        }
                                        points.push(p);
                                    }
                                    mults[0] = 2;
                                    mults[nbp - 1] = 2;
                                    // occ::handle<Geom_BSplineCurve> AppC;
                                    // occ::handle<Geom2d_BSplineCurve> AppC2d;
                                    // AppC = new Geom_BSplineCurve(Points,
                                    //   knots, mults, 1).
                                    let app_c = Arc::new(BSplineCurve3 {
                                        degree: 1,
                                        knots: expand_knots(&knots, &mults),
                                        control_points: points,
                                        weights: vec![],
                                        is_periodic: false,
                                    });

                                    let mut app_c2d: Option<Arc<BSplineCurve2>> = None;
                                    if with_pcurve {
                                        let mut points2d: Vec<DVec2> =
                                            Vec::with_capacity(nbp);
                                        for i in 1..=nbp {
                                            let (u, v) = line
                                                .point(i + ip_f - 1)
                                                .parameters_on_surface(false);
                                            if u < min_u {
                                                min_u = u;
                                            }
                                            if v < min_v {
                                                min_v = v;
                                            }
                                            if u > max_u {
                                                max_u = u;
                                            }
                                            if v > max_v {
                                                max_v = v;
                                            }
                                            points2d.push(DVec2::new(u, v));
                                        }
                                        // AppC2d = new Geom2d_BSplineCurve(
                                        //   Points2d, knots, mults, 1).
                                        app_c2d = Some(Arc::new(BSplineCurve2 {
                                            degree: 1,
                                            knots: {
                                                let kn = knots.clone();
                                                let mu = mults.clone();
                                                expand_knots(&kn, &mu)
                                            },
                                            control_points: points2d,
                                            weights: vec![],
                                        }));
                                    }
                                    first = 1.0;
                                    last = nbp as f64;

                                    // occ::handle<BRepApprox_ApproxLine> AppLine;
                                    // occ::handle<Geom2d_BSplineCurve> CNull;
                                    // AppLine = new BRepApprox_ApproxLine(
                                    //   AppC, AppC2d, CNull).
                                    let app_line = Arc::new(ApproxLine::new(
                                        Some(app_c.clone()),
                                        app_c2d
                                            .as_ref()
                                            .map(|x| Arc::new((**x).clone())),
                                        None,
                                    ));

                                    // int dmin = 4, dmax = 8, niter = 0;
                                    // bool tg = false;
                                    let dmin: usize = 4;
                                    let dmax: usize = 8;
                                    let niter: i32 = 0;
                                    let tg = false;
                                    let mut approx = BRepApproxApprox::new();
                                    // double TOL3d, TOL2d, TOL = 0.0001;
                                    let tol_param = 0.0001;

                                    // Maxx -= Minx; Maxy -= Miny;
                                    // Maxz -= Minz; Maxu -= Minu;
                                    // Maxv -= Minv;
                                    max_x -= min_x;
                                    max_y -= min_y;
                                    max_z -= min_z;
                                    max_u -= min_u;
                                    max_v -= min_v;
                                    // if (Maxy > Maxx) { Maxx = Maxy; }
                                    if max_y > max_x {
                                        max_x = max_y;
                                    }
                                    // if (Maxz > Maxx) { Maxx = Maxy; } — the
                                    // OCCT literal (Maxy, not Maxz).
                                    if max_z > max_x {
                                        max_x = max_y;
                                    }
                                    // if (Maxv > Maxu) { Maxu = Maxv; }
                                    if max_v > max_u {
                                        max_u = max_v;
                                    }

                                    // TOL3d = TOL * Maxx;
                                    let mut tol3d = tol_param * max_x;
                                    if tol3d < 1e-12 {
                                        tol3d = 1e-12;
                                    } else if tol3d > 0.1 {
                                        tol3d = 0.1;
                                    }
                                    // TOL2d = TOL * Maxu;
                                    let mut tol2d = tol_param * max_u;
                                    if tol2d < 1e-12 {
                                        tol2d = 1e-12;
                                    } else if tol2d > 0.1 {
                                        tol2d = 0.1;
                                    }

                                    //-- std::cout<<"\nHLRTopoBRep_DSFiller :
                                    //-- nbp="<<nbp<<"  Tol3d="<<TOL3d<<"
                                    //-- Tol2d="<<TOL2d<<std::endl; — commented
                                    //-- out.

                                    // Approx.SetParameters(TOL3d, TOL2d, dmin,
                                    //   dmax, niter, 30, tg) — the hxx default
                                    //   Parametrization = Approx_ChordLength
                                    //   (BRepApprox_Approx.hxx L111).
                                    approx.set_parameters(
                                        tol3d,
                                        tol2d,
                                        dmin,
                                        dmax,
                                        niter,
                                        30,
                                        tg,
                                        ApproxParamType::ChordLength,
                                    );
                                    // Approx.Perform(AppLine, true, true,
                                    //   false, 1, nbp).
                                    approx.perform_wline(
                                        app_line, true, true, false, 1, nbp,
                                    );
                                    if !approx.is_done() {
                                        // C = AppC; C2d = AppC2d;
                                        // first = 1; last = nbp.
                                        c = Some(Curve3::BSpline((*app_c).clone()));
                                        c2d = app_c2d
                                            .as_ref()
                                            .map(|x| Curve2d::BSpline((**x).clone()));
                                        first = 1.0;
                                        last = nbp as f64;
                                    } else {
                                        // const AppParCurves_MultiBSpCurve&
                                        //   AppVal = Approx.Value(1) — the
                                        //   Index argument is ignored
                                        //   (ApproxInt_Approx.gxx L487-497):
                                        //   Value(2) returns the same curve,
                                        //   so one rcad value() call stands
                                        //   for both.
                                        let app_val = approx.value();
                                        // NCollection_Array1<gp_Pnt>
                                        //   poles3d(1, AppVal.NbPoles());
                                        // AppVal.Curve(1, poles3d);
                                        let mut poles3d: Vec<DVec3> = Vec::new();
                                        app_val.curve(1, &mut poles3d);
                                        let knots_flat =
                                            expand_knots(&app_val.knots, &app_val.mults);
                                        // C = new Geom_BSplineCurve(poles3d,
                                        //   AppVal.Knots(),
                                        //   AppVal.Multiplicities(),
                                        //   AppVal.Degree()).
                                        c = Some(Curve3::BSpline(BSplineCurve3 {
                                            degree: app_val.degree,
                                            knots: knots_flat.clone(),
                                            control_points: poles3d,
                                            weights: vec![],
                                            is_periodic: false,
                                        }));

                                        // const AppParCurves_MultiBSpCurve&
                                        //   AppVal2 = Approx.Value(2) — the
                                        //   same MultiBSpCurve;
                                        //   AppVal2.Curve(2, poles2d) — the
                                        //   global 2d slot (2 > the single 3d
                                        //   slot) is the first pcurve.
                                        let mut poles2d: Vec<DVec2> = Vec::new();
                                        app_val.curve2d(2, &mut poles2d);
                                        // C2d = new Geom2d_BSplineCurve(
                                        //   poles2d, AppVal2.Knots(),
                                        //   AppVal2.Multiplicities(),
                                        //   AppVal2.Degree()).
                                        let c2d_bs = BSplineCurve2 {
                                            degree: app_val.degree,
                                            knots: knots_flat,
                                            control_points: poles2d,
                                            weights: vec![],
                                        };
                                        // first = C2d->FirstParameter();
                                        // last = C2d->LastParameter();
                                        first = c2d_bs.knots[c2d_bs.degree];
                                        last = c2d_bs.knots
                                            [c2d_bs.knots.len() - 1 - c2d_bs.degree];
                                        c2d = Some(Curve2d::BSpline(c2d_bs));
                                    }
                                }
                            }
                            IType::Restriction => {
                                // throw Standard_ProgramError(
                                //   "HLRTopoBRep_DSFiller::InsertFace :
                                //   Restriction").
                                panic!(
                                    "Standard_ProgramError: HLRTopoBRep_DSFiller::InsertFace : Restriction"
                                );
                            }
                        }

                        // compute the PCurve
                        // make the edge
                        if !insuffisant_number_of_points {
                            // TopoDS_Edge E; BRep_Builder B;
                            // B.MakeEdge(E, C, tol) — the new edge TShape
                            // carrying C with edge tolerance tol; the
                            // B.Add / B.Range below map onto the combined
                            // rcad add_edge (start vertex FORWARD, end
                            // vertex REVERSED, the parameter range).
                            let c = c.expect("the switch arms assign the curve");
                            vf.orientation = Orientation::Forward;
                            vl.orientation = Orientation::Reversed;
                            let e = b.add_edge(
                                brep,
                                Some(c),
                                vf.clone(),
                                vl.clone(),
                                [first, last],
                            );
                            // OCCT MakeEdge(E, C, tol) — the edge tolerance
                            // is tol.
                            brep.edge_mut_inplace(e.clone()).tolerance = tol;
                            // B.Range(E, first, last) — the add_edge range
                            // argument.

                            if let Some(c2d) = c2d {
                                // B.UpdateEdge(E, C2d, F,
                                //   BRep_Tool::Tolerance(F)).
                                b.update_edge_pcurve(
                                    brep,
                                    e.clone(),
                                    c2d,
                                    f.clone(),
                                    brep.tolerance(f),
                                );
                            }

                            // add the edge in the DS
                            // if (!E.IsNull()) — a rcad edge handle is never
                            // null.
                            {
                                ds.add_int_l(f).push(e);
                            }
                        }
                    }
                }
            }
        }
    }

    // Correction of internal outlines: unite coinciding vertices
    // const double SqTol = tol * tol;
    let sq_tol = tol * tol;
    let int_l = ds.add_int_l(f).clone();
    // NCollection_List<TopoDS_Shape>::Iterator itl1(IntL);
    for itl1 in 0..int_l.len() {
        // TopoDS_Edge anIntLine = TopoDS::Edge(itl1.Value());
        let mut an_int_line = int_l[itl1].clone();
        an_int_line.orientation = Orientation::Forward;
        // TopoDS_Vertex aVer[2];
        // TopExp::Vertices(anIntLine, aVer[0], aVer[1]).
        let a_ver = [brep.first_vertex(&an_int_line), brep.last_vertex(&an_int_line)];
        // NCollection_List<TopoDS_Shape>::Iterator itl2 = itl1;
        for itl2 in itl1..int_l.len() {
            // TopoDS_Edge anIntLine2 = TopoDS::Edge(itl2.Value());
            let mut an_int_line2 = int_l[itl2].clone();
            an_int_line2.orientation = Orientation::Forward;
            if an_int_line2.is_same(&an_int_line) {
                continue;
            }
            // TopoDS_Vertex aVer2[2];
            // TopExp::Vertices(anIntLine2, aVer2[0], aVer2[1]).
            let a_ver2 = [
                brep.first_vertex(&an_int_line2),
                brep.last_vertex(&an_int_line2),
            ];
            for i in 0..2 {
                if i == 1 && a_ver[0].is_same(&a_ver[1]) {
                    continue;
                }
                // gp_Pnt Pnt1 = BRep_Tool::Pnt(aVer[i]);
                let pnt1 = brep.vertex_position(&a_ver[i]);
                for j in 0..2 {
                    if a_ver[i].is_same(&a_ver2[j]) {
                        continue;
                    }
                    // gp_Pnt Pnt2 = BRep_Tool::Pnt(aVer2[j]);
                    let pnt2 = brep.vertex_position(&a_ver2[j]);
                    if pnt1.distance_squared(pnt2) <= sq_tol {
                        // BRep_Builder aBB;
                        // aBB.Remove(anIntLine2, aVer2[j]);
                        b.remove_from_edge(brep, an_int_line2.clone(), a_ver2[j].clone());
                        // aVer[i].Orientation(
                        //   (j == 0) ? TopAbs_FORWARD : TopAbs_REVERSED);
                        let mut v = a_ver[i].clone();
                        v.orientation = if j == 0 {
                            Orientation::Forward
                        } else {
                            Orientation::Reversed
                        };
                        // aBB.Add(anIntLine2, aVer[i]);
                        b.add_to_edge(brep, an_int_line2.clone(), v);
                    }
                }
            }
        }
    }
}

//=======================================================================
// function : MakeVertex
// purpose  : private, make a vertex from an intersection point
// OCCT HLRTopoBRep_DSFiller.cxx L599-659
//=======================================================================

/// OCCT HLRTopoBRep_DSFiller::MakeVertex(P, tol, DS).
pub(crate) fn make_vertex(brep: &mut BRep, p: &ContapPoint, tol: f64, ds: &mut Data) -> Shape {
    let mut b = BRepBuilder::new();
    // TopoDS_Vertex V;
    let mut v = Shape::null();
    if p.is_vertex() {
        // V = occ::down_cast<BRepTopAdaptor_HVertex>(P.Vertex())->Vertex().
        v = p
            .vertex()
            .topo_vertex()
            .expect("Standard_TypeMismatch: BRepTopAdaptor_HVertex")
            .clone();
        ds.add_out_v(&v);
    } else {
        // if on arc, insert in the DS
        if p.is_on_arc() {
            // const TopoDS_Edge& E =
            //   (*(BRepAdaptor_Curve2d*)(P.Arc().get())).Edge().
            let e = p
                .arc()
                .as_any()
                .downcast_ref::<BRepCurve2d>()
                .expect("Standard_TypeMismatch: BRepAdaptor_Curve2d")
                .edge()
                .clone();
            let par = p.parameter_on_arc();
            let p3d = p.value();

            // for (DS.InitVertex(E); DS.MoreVertex(); DS.NextVertex()).
            ds.init_vertex(&e);
            while ds.more_vertex() {
                // TopoDS_Vertex curV = DS.Vertex();
                let cur_v = ds.vertex();
                // double curP = DS.Parameter();
                let cur_p = ds.parameter();
                // const gp_Pnt& PPP = BRep_Tool::Pnt(curV);
                let ppp = brep.vertex_position(&cur_v);
                // double TTT = BRep_Tool::Tolerance(curV);
                let ttt = brep.vertex_tolerance(&cur_v);
                // if (P3d.IsEqual(PPP, TTT)).
                if p3d.distance(ppp) <= ttt {
                    v = cur_v;
                    break;
                } else if par < cur_p {
                    // B.MakeVertex(V, P.Value(), tol).
                    v = b.add_vertex(brep, p.value(), tol);
                    ds.insert_before(&v, par);
                    break;
                }
                ds.next_vertex();
            }
            if !ds.more_vertex() {
                // B.MakeVertex(V, P.Value(), tol).
                v = b.add_vertex(brep, p.value(), tol);
                ds.append(&v, par);
            }
            ds.add_out_v(&v);
        }
        // if internal create a vertex and insert in the DS
        else {
            // B.MakeVertex(V, P.Value(), tol).
            v = b.add_vertex(brep, p.value(), tol);
            if p.is_internal() {
                ds.add_int_v(&v);
            } else {
                ds.add_out_v(&v);
            }
        }
    }
    v
}

//=======================================================================
// function : InsertVertex
// purpose  : private, insert a vertex from an internal intersection point
//           on resctriction
// OCCT HLRTopoBRep_DSFiller.cxx L667-706
//=======================================================================

/// OCCT HLRTopoBRep_DSFiller::InsertVertex(P, tol, E, DS).
pub(crate) fn insert_vertex(brep: &mut BRep, p: &ContapPoint, tol: f64, e: &Shape, ds: &mut Data) {
    let mut b = BRepBuilder::new();
    // TopoDS_Vertex V;
    let mut v = Shape::null();

    if p.is_vertex() {
        // V = occ::down_cast<BRepTopAdaptor_HVertex>(P.Vertex())->Vertex().
        v = p
            .vertex()
            .topo_vertex()
            .expect("Standard_TypeMismatch: BRepTopAdaptor_HVertex")
            .clone();
    } else {
        // double Par = P.ParameterOnLine();
        let par = p.parameter_on_line();

        // for (DS.InitVertex(E); DS.MoreVertex(); DS.NextVertex()).
        ds.init_vertex(e);
        while ds.more_vertex() {
            let cur_v = ds.vertex();
            let cur_p = ds.parameter();
            // if (P.Value().IsEqual(BRep_Tool::Pnt(curV),
            //   BRep_Tool::Tolerance(curV))).
            if p.value().distance(brep.vertex_position(&cur_v))
                <= brep.vertex_tolerance(&cur_v)
            {
                v = cur_v;
                break;
            } else if par < cur_p {
                // B.MakeVertex(V, P.Value(), tol).
                v = b.add_vertex(brep, p.value(), tol);
                ds.insert_before(&v, par);
                break;
            }
            ds.next_vertex();
        }
        if !ds.more_vertex() {
            // B.MakeVertex(V, P.Value(), tol).
            v = b.add_vertex(brep, p.value(), tol);
            ds.append(&v, par);
        }
    }
    ds.add_int_v(&v);
}

//=======================================================================
// function : ProcessEdges
// purpose  : private, split edges with outline vertices
// OCCT HLRTopoBRep_DSFiller.cxx L713-758
//=======================================================================

/// OCCT HLRTopoBRep_DSFiller::ProcessEdges(DS).
pub(crate) fn process_edges(brep: &mut BRep, ds: &mut Data) {
    let mut b = BRepBuilder::new();
    // TopoDS_Edge newE; TopoDS_Vertex VF, VL, VI;
    // double PF, PL, PI; — the loop-carried locals are (re)assigned per
    // iteration below.

    // for (DS.InitEdge(); DS.MoreEdge(); DS.NextEdge()).
    ds.init_edge();
    while ds.more_edge() {
        // TopoDS_Edge E = DS.Edge();
        let e = ds.edge();
        // NCollection_List<TopoDS_Shape>& SplE = DS.AddSplE(E); — the rcad
        // fetches re-address the same bound list (add_spl_e is idempotent).
        ds.add_spl_e(&e);
        // TopExp::FirstVertex(E) / TopExp::LastVertex(E).
        let mut vf = brep.first_vertex(&e);
        let mut vl = brep.last_vertex(&e);
        // BRep_Tool::Range(E, PF, PL).
        let range = brep.edge_range(&e);
        let mut pf = range[0];
        let pl = range[1];
        vf.orientation = Orientation::Forward;
        vl.orientation = Orientation::Reversed;

        // for (DS.InitVertex(E); DS.MoreVertex(); DS.NextVertex()).
        ds.init_vertex(&e);
        while ds.more_vertex() {
            // VI = DS.Vertex(); PI = DS.Parameter();
            let mut vi = ds.vertex();
            let pi = ds.parameter();
            vi.orientation = Orientation::Reversed;
            // newE = E; newE.EmptyCopy();
            let mut new_e = brep.empty_copy(e.clone());
            new_e.orientation = Orientation::Forward;
            // B.Add(newE, VF);
            b.add_to_edge(brep, new_e.clone(), vf.clone());
            // B.UpdateVertex(VF, PF, newE, BRep_Tool::Tolerance(VF)).
            b.update_vertex_on_edge(
                brep,
                vf.clone(),
                pf,
                new_e.clone(),
                brep.vertex_tolerance(&vf),
            );
            // B.Add(newE, VI);
            b.add_to_edge(brep, new_e.clone(), vi.clone());
            // B.UpdateVertex(VI, PI, newE, BRep_Tool::Tolerance(VI)).
            b.update_vertex_on_edge(
                brep,
                vi.clone(),
                pi,
                new_e.clone(),
                brep.vertex_tolerance(&vi),
            );
            // newE.Orientation(E.Orientation()).
            new_e.orientation = e.orientation;
            // SplE.Append(newE).
            ds.add_spl_e(&e).push(new_e);
            // VF = VI; PF = PI; VF.Orientation(TopAbs_FORWARD).
            vf = vi;
            pf = pi;
            vf.orientation = Orientation::Forward;
            ds.next_vertex();
        }
        // newE = E; newE.EmptyCopy();
        let mut new_e = brep.empty_copy(e.clone());
        new_e.orientation = Orientation::Forward;
        // B.Add(newE, VF);
        b.add_to_edge(brep, new_e.clone(), vf.clone());
        // B.UpdateVertex(VF, PF, newE, BRep_Tool::Tolerance(VF)).
        b.update_vertex_on_edge(
            brep,
            vf.clone(),
            pf,
            new_e.clone(),
            brep.vertex_tolerance(&vf),
        );
        // B.Add(newE, VL);
        b.add_to_edge(brep, new_e.clone(), vl.clone());
        // B.UpdateVertex(VL, PL, newE, BRep_Tool::Tolerance(VL)).
        b.update_vertex_on_edge(
            brep,
            vl.clone(),
            pl,
            new_e.clone(),
            brep.vertex_tolerance(&vl),
        );
        // newE.Orientation(E.Orientation()).
        new_e.orientation = e.orientation;
        // SplE.Append(newE).
        ds.add_spl_e(&e).push(new_e);

        ds.next_edge();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::int_patch::transitions::Transition;
    use crate::geomalgo::int_surf::{LineOn2S, PntOn2S};
    use glam::DVec3;

    /// The test rig: the unit square face of the topol_tool_brep fixture
    /// (vertices (0,0), (1,0), (1,1), (0,1) on z = 0; pcurves = XY) with
    /// its first wire edge e1 (the bottom edge, range [0, 1]) exposed.
    fn fixture() -> (BRep, Shape, Shape, Vec<Shape>) {
        let (brep, face) =
            crate::topalgo::brep_top_adaptor::topol_tool_brep::tests::square_face();
        let wire = brep.face(face.clone()).outer_wire.clone();
        let edges = brep.wire(wire).edges.clone();
        (brep, face, edges[0].clone(), edges)
    }

    /// A Contap point standing on the restriction arc of (e1, face) at the
    /// given 3d value and arc/line parameter.
    fn on_arc_point(brep: &BRep, e: &Shape, f: &Shape, value: DVec3, par: f64) -> ContapPoint {
        let bc2d = BRepCurve2d::new_edge_face(Arc::new(brep.clone()), e, f);
        let mut p = ContapPoint::with_uv(value, 0.0, 0.0);
        if par >= 0.0 {
            p.set_arc(Arc::new(bc2d), par, Transition::new(), Transition::new());
        }
        p
    }

    /// OCCT anchor: ProcessEdges (cxx L713-758) — an edge carrying two
    /// middle vertices splits into three segments; every segment is an
    /// EmptyCopy (a distinct TShape), carries the segment range, the
    /// FORWARD/REVERSED endpoint vertices and the original orientation.
    #[test]
    fn process_edges_splits_at_vertices() {
        let (mut brep, _face, _e1, es) = fixture();
        let e1 = es[0].clone();
        let v1 = brep.first_vertex(&e1);
        let vlast = brep.last_vertex(&e1);
        let mut ds = Data::new();

        // Two outline vertices on e1 at t = 0.25 and t = 0.75 (the Data
        // vertex list of an edge carries the middle outline vertices).
        let mut b = BRepBuilder::new();
        let mid1 = b.add_vertex(&mut brep, DVec3::new(0.25, 0.0, 0.0), 1e-7);
        let mid2 = b.add_vertex(&mut brep, DVec3::new(0.75, 0.0, 0.0), 1e-7);
        ds.init_vertex(&e1);
        ds.append(&mid1, 0.25);
        ds.append(&mid2, 0.75);

        process_edges(&mut brep, &mut ds);

        let spl_e = ds.edge_spl_e(&e1).clone();
        assert_eq!(spl_e.len(), 3);
        // Segment 0: v1 (0.0) -> mid1 (0.25).
        let seg0 = &spl_e[0];
        assert_eq!(brep.edge_range(seg0), [0.0, 0.25]);
        assert!(brep.first_vertex(seg0).is_same(&v1));
        assert!(brep.last_vertex(seg0).is_same(&mid1));
        assert_eq!(seg0.orientation, Orientation::Forward);
        // The range carries the vertex parameters (the FORWARD child updated
        // the First, the REVERSED child the Last).
        let seg1 = &spl_e[1];
        assert_eq!(brep.edge_range(seg1), [0.25, 0.75]);
        assert!(brep.first_vertex(seg1).is_same(&mid1));
        assert!(brep.last_vertex(seg1).is_same(&mid2));
        let seg2 = &spl_e[2];
        assert_eq!(brep.edge_range(seg2), [0.75, 1.0]);
        assert!(brep.first_vertex(seg2).is_same(&mid2));
        assert!(brep.last_vertex(seg2).is_same(&vlast));
        // The segments are new TShapes (EmptyCopy), not the original edge.
        assert!(!seg0.is_same(&e1));
        // The vertex orientations: start FORWARD, end REVERSED.
        let sh = brep.edge(seg0.clone()).my_shapes.clone();
        assert_eq!(sh.len(), 2);
        assert_eq!(sh[0].orientation, Orientation::Forward);
        assert_eq!(sh[1].orientation, Orientation::Reversed);
    }

    /// OCCT anchor: MakeVertex (cxx L599-659) — a point on an arc is merged
    /// with an existing vertex at the same position (IsEqual over the stored
    /// tolerance), otherwise inserted at the parameter-ordered position
    /// (InsertBefore) or appended after the last.
    #[test]
    fn make_vertex_dedup_and_parameter_order() {
        let (mut brep, face, e1, _es) = fixture();
        let v1 = brep.first_vertex(&e1);
        let v2 = brep.last_vertex(&e1);
        let mut ds = Data::new();

        // The stored vertex list of e1: v1 @ 0.0, v2 @ 1.0.
        ds.init_vertex(&e1);
        ds.append(&v1, 0.0);
        ds.append(&v2, 1.0);

        // 1) The dedup path: the point sits exactly on v2, its parameter is
        //    smaller than v2's — the IsEqual branch wins.
        let p_dup = on_arc_point(&brep, &e1, &face, DVec3::new(1.0, 0.0, 0.0), 0.9);
        let v = make_vertex(&mut brep, &p_dup, 1e-7, &mut ds);
        assert!(v.is_same(&v2));
        assert!(ds.is_out_v(&v));

        // 2) The InsertBefore path: a new point at (0.25, 0, 0) with
        //    parameter 0.5 — larger than v1's 0.0, smaller than v2's 1.0.
        let p_mid = on_arc_point(&brep, &e1, &face, DVec3::new(0.25, 0.0, 0.0), 0.5);
        let v = make_vertex(&mut brep, &p_mid, 1e-7, &mut ds);
        assert!(!v.is_same(&v1) && !v.is_same(&v2));
        ds.init_vertex(&e1);
        assert!(ds.vertex().is_same(&v1));
        assert_eq!(ds.parameter(), 0.0);
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v));
        assert_eq!(ds.parameter(), 0.5);
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v2));

        // 3) The Append path: a point past every stored parameter.
        let p_end = on_arc_point(&brep, &e1, &face, DVec3::new(0.9, 0.0, 0.0), 1.5);
        let v = make_vertex(&mut brep, &p_end, 1e-7, &mut ds);
        ds.init_vertex(&e1);
        let mut last_v = Shape::null();
        let mut last_p = -1.0;
        while ds.more_vertex() {
            last_v = ds.vertex();
            last_p = ds.parameter();
            ds.next_vertex();
        }
        assert!(last_v.is_same(&v));
        assert_eq!(last_p, 1.5);
        assert!(ds.is_out_v(&v));
    }

    /// OCCT anchor: MakeVertex — an internal (not on any arc) point creates
    /// a fresh vertex registered in the IntV set; a non-internal one in the
    /// OutV set.
    #[test]
    fn make_vertex_internal_and_outline() {
        let (mut brep, _face, _vs, _es) = fixture();
        let mut ds = Data::new();

        let p_int = {
            let mut q = ContapPoint::with_uv(DVec3::new(0.5, 0.5, 0.0), 0.0, 0.0);
            q.set_internal();
            q
        };
        let v = make_vertex(&mut brep, &p_int, 1e-7, &mut ds);
        assert!(!v.is_null());
        assert_eq!(brep.vertex_position(&v), DVec3::new(0.5, 0.5, 0.0));
        assert!(ds.is_int_v(&v));
        assert!(!ds.is_out_v(&v));

        let p_out = ContapPoint::with_uv(DVec3::new(0.25, 0.25, 0.0), 0.0, 0.0);
        let v = make_vertex(&mut brep, &p_out, 1e-7, &mut ds);
        assert!(ds.is_out_v(&v));
        assert!(!ds.is_int_v(&v));
    }

    /// OCCT anchor: InsertVertex (cxx L667-706) — the parameter-ordered
    /// insertion into the restriction edge's vertex list and the IntV
    /// registration.
    #[test]
    fn insert_vertex_parameter_order_and_int_v() {
        let (mut brep, _face, e1, _es) = fixture();
        let v1 = brep.first_vertex(&e1);
        let v2 = brep.last_vertex(&e1);
        let mut ds = Data::new();
        ds.init_vertex(&e1);
        ds.append(&v1, 0.0);
        ds.append(&v2, 1.0);

        // A new point between v1 and v2 at line parameter 0.5 (no arc —
        // InsertVertex reads ParameterOnLine for the non-vertex case).
        let p = {
            let mut q = ContapPoint::with_uv(DVec3::new(0.25, 0.0, 0.0), 0.0, 0.0);
            q.set_parameter(0.5);
            q
        };
        insert_vertex(&mut brep, &p, 1e-7, &e1, &mut ds);

        ds.init_vertex(&e1);
        assert!(ds.vertex().is_same(&v1));
        ds.next_vertex();
        let inserted = ds.vertex();
        assert_eq!(ds.parameter(), 0.5);
        ds.next_vertex();
        assert!(ds.vertex().is_same(&v2));
        // InsertVertex registers the vertex as an internal outline vertex.
        assert!(ds.is_int_v(&inserted));
    }

    /// OCCT anchor: InsertFace's Restriction arm (cxx L154-176) — the arc
    /// edge lands in OutL and the internal points are gated by the double
    /// IsEqual(VF) && IsEqual(VL) check; and the Walking <5-point arm
    /// (cxx L321-350) — a degree-1 BSpline over the walked points with
    /// knots 1..nbp, end multiplicities 2, range [1, nbp].
    #[test]
    fn insert_face_restriction_and_walking_degree1() {
        let (mut brep, face, e1, _es) = fixture();
        let mut fo = Contour::new();
        fo.set_done(true);
        let mut ds = Data::new();

        // --- Restriction line on the face's bottom edge ---
        let bc2d = BRepCurve2d::new_edge_face(Arc::new(brep.clone()), &e1, &face);
        let mut rline = crate::hlr::contap::line::Line::new();
        rline.set_value_arc(Arc::new(bc2d));
        // An internal point not coincident with both endpoints: the gate
        // keeps it out of the vertex list (the OutL append still happens).
        let mut p = ContapPoint::with_uv(DVec3::new(0.25, 0.0, 0.0), 0.0, 0.0);
        p.set_internal();
        p.set_parameter(0.5);
        rline.add(p);
        fo.slin_mut().push(rline);

        insert_face(&mut brep, 1, &face, &mut fo, &mut ds, true);
        assert_eq!(ds.face_out_l(&face).len(), 1);
        assert!(ds.face_out_l(&face)[0].is_same(&e1));
        assert!(ds.face_int_l(&face).is_empty());
        assert!(ds.face_has_out_l(&face));

        // --- Walking line with ipL - ipF = 3 (< 5): the degree-1 branch ---
        let mut fo2 = Contour::new();
        fo2.set_done(true);
        let mut ds2 = Data::new();

        let mut walking = crate::hlr::contap::line::Line::new();
        let mut pts = LineOn2S::new();
        // Six walking points on z = 0, UV = XY.
        let walk3d = [
            DVec3::new(0.1, 0.1, 0.0),
            DVec3::new(0.2, 0.1, 0.0),
            DVec3::new(0.3, 0.2, 0.0),
            DVec3::new(0.4, 0.3, 0.0),
            DVec3::new(0.5, 0.4, 0.0),
            DVec3::new(0.6, 0.5, 0.0),
        ];
        for (k, w) in walk3d.iter().enumerate() {
            let mut p2 = PntOn2S::new();
            p2.set_value_all(*w, w.x, w.y, w.x, w.y);
            pts.add(&p2);
        }
        walking.set_line_on_2s(&pts);
        // The two Contap vertices delimit the walked range [1, 4]:
        // ipF = 1, ipL = 4 -> ipL - ipF = 3 < 5.
        let mut pf = ContapPoint::with_uv(walk3d[1], 0.0, 0.0);
        pf.set_parameter(1.0);
        let mut pl = ContapPoint::with_uv(walk3d[4], 0.0, 0.0);
        pl.set_parameter(4.0);
        walking.add(pf);
        walking.add(pl);
        fo2.slin_mut().push(walking);

        insert_face(&mut brep, 1, &face, &mut fo2, &mut ds2, true);

        let int_l = ds2.face_int_l(&face);
        assert_eq!(int_l.len(), 1);
        let new_e = int_l[0].clone();
        // The walked range: parF = 1, parL = 4.
        assert_eq!(brep.edge_range(&new_e), [1.0, 4.0]);
        // The curve: degree 1 over the points ipF..ipL (four poles).
        let curve = brep.edge_curve_data(&new_e).expect("edge 3d curve");
        match curve {
            Curve3::BSpline(bs) => {
                assert_eq!(bs.degree, 1);
                assert_eq!(bs.control_points.len(), 4);
                assert_eq!(bs.knots, vec![1.0, 1.0, 2.0, 3.0, 4.0, 4.0]);
                // The poles are the walking points (the same accessor).
                let line_ref = fo2.line(1);
                for i in 0..4 {
                    assert_eq!(bs.control_points[i], line_ref.point(i + 1).value());
                }
            }
            _ => panic!("expected a BSpline walking edge"),
        }
        // The two Contap vertices became the edge endpoints.
        let seg_vs = brep.edge(new_e.clone()).my_shapes.clone();
        assert_eq!(seg_vs.len(), 2);
        // The pcurve: the 2d BSpline over the walked points (the kernel
        // pc_parameter_range stores the normalized [0, 1] range for a
        // BSpline pcurve — the OCCT UpdateEdge would carry the pcurve knot
        // domain [1, 4]; noted as the kernel approximation).
        let (pc, _pc_first, _pc_last) = brep
            .curve_on_surface(&new_e, &face)
            .expect("walking pcurve");
        match pc {
            Curve2d::BSpline(bs2) => {
                assert_eq!(bs2.degree, 1);
                assert_eq!(bs2.control_points.len(), 4);
                assert_eq!(bs2.knots, vec![1.0, 1.0, 2.0, 3.0, 4.0, 4.0]);
            }
            _ => panic!("expected a BSpline walking pcurve"),
        }
    }

    /// OCCT anchor: the Walking >=5-point arm (cxx L351-507) — the
    /// interpolation BSpline feeds BRepApprox_Approx; either the
    /// approximation result or the AppC fallback lands as a BSpline edge in
    /// IntL with a pcurve.
    #[test]
    fn insert_face_walking_approx_arm() {
        let (mut brep, face, _e1, _es) = fixture();
        let mut fo = Contour::new();
        fo.set_done(true);
        let mut ds = Data::new();

        let mut walking = crate::hlr::contap::line::Line::new();
        let mut pts = LineOn2S::new();
        let n = 8; // ipF = 1, ipL = 8 -> 7 >= 5 points
        for k in 0..n {
            let x = 0.1 + 0.1 * k as f64;
            let y = 0.1 + 0.05 * (k as f64) * (k as f64);
            let mut p2 = PntOn2S::new();
            p2.set_value_all(DVec3::new(x, y, 0.0), x, y, x, y);
            pts.add(&p2);
        }
        walking.set_line_on_2s(&pts);
        let mut pf = ContapPoint::with_uv(DVec3::new(0.1, 0.1, 0.0), 0.0, 0.0);
        pf.set_parameter(1.0);
        let mut pl = ContapPoint::with_uv(
            DVec3::new(0.1 + 0.1 * (n - 1) as f64, 0.1 + 0.05 * ((n - 1) as f64) * ((n - 1) as f64), 0.0),
            0.0,
            0.0,
        );
        pl.set_parameter((n - 1) as f64);
        walking.add(pf);
        walking.add(pl);
        fo.slin_mut().push(walking);

        insert_face(&mut brep, 1, &face, &mut fo, &mut ds, true);

        let int_l = ds.face_int_l(&face);
        assert_eq!(int_l.len(), 1);
        let new_e = int_l[0].clone();
        let curve = brep.edge_curve_data(&new_e).expect("edge 3d curve");
        match curve {
            Curve3::BSpline(bs) => {
                assert!(bs.control_points.len() >= 2);
                assert!(bs.knots.len() >= 2);
            }
            _ => panic!("expected a BSpline walking edge"),
        }
        // The pcurve exists (withPCurve).
        assert!(brep.curve_on_surface(&new_e, &face).is_some());
        // The endpoints carry the Contap vertex positions.
        let seg_vs = brep.edge(new_e.clone()).my_shapes.clone();
        assert_eq!(seg_vs.len(), 2);
    }

    /// OCCT anchor: Insert's MST cache (cxx L84-96) — a pre-bound face tool
    /// is reused (no second Bind) and the tool's surface drives
    /// Contap_Contour::Perform; the DS is cleared on entry and
    /// ProcessEdges closes the pass.
    #[test]
    fn insert_mst_cache_reuse() {
        let (mut brep, face, _vs, _es) = fixture();
        let mut ds = Data::new();
        let mut mst: Vec<(Shape, BRepTopAdaptorTool)> = Vec::new();

        // Pre-bind the tool for the face (a silhouette direction along the
        // plane normal produces no contour lines: the InsertFace step is
        // skipped and the MST bind path is what the test observes).
        let pre_bound = BRepTopAdaptorTool::new_face(Arc::new(brep.clone()), &face, 1e-7);
        mst_bind(&mut mst, &face, pre_bound);

        let mut fo = Contour::with_direction(DVec3::new(0.0, 0.0, 1.0));
        insert(&mut brep, &face, &mut fo, &mut ds, &mut mst, 0);

        // The pre-bound tool survived: no second Bind for the same face.
        assert_eq!(mst.len(), 1);
        assert!(mst_is_bound(&mst, &face));
        assert!(fo.is_done());
        assert!(fo.is_empty());
        assert!(!ds.face_has_int_l(&face));
        assert!(!ds.face_has_out_l(&face));
    }

    /// OCCT anchor: Insert (cxx L62-113) — an unbound face gets a fresh
    /// BRepTopAdaptor_Tool (Bind) and nbIso == 0 skips FaceIsoLiner.
    #[test]
    fn insert_binds_unbound_face() {
        let (mut brep, face, _vs, _es) = fixture();
        let mut ds = Data::new();
        let mut mst: Vec<(Shape, BRepTopAdaptorTool)> = Vec::new();
        let mut fo = Contour::with_direction(DVec3::new(0.0, 0.0, 1.0));

        insert(&mut brep, &face, &mut fo, &mut ds, &mut mst, 0);
        assert_eq!(mst.len(), 1);
        assert!(mst_is_bound(&mst, &face));
        assert!(fo.is_done());
    }
}
