// OCCT BRepFeat.cxx statics + vehicles consumed by BRepFeat_Form — 1:1
// translation (split from brep_feat_form.rs per the 2000-line rule; the OCCT
// line anchors are the same files as the parent module header).
//
// Content: the re-hosted BRepFeat package statics (ParametricMinMax /
// IsInside / FaceUntil / Tool — BRepFeat.cxx), the BRepAlgo::IsValid gap
// marker (BRepAlgo_1.cxx L39-43), the BRepAlgoAPI_Cut vehicle (CutVehicle)
// and the result-pool helpers (pool_top_shapes / builder_result_shape). The
// architecture-difference numbering of the parent module header applies.

use crate::bop::algo::builder::BooleanOpType;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::history::BRepToolsHistory;
use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder, OcctShapeMap};
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, TrimmedSurface};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};

/// OCCT NECHANTBARYC (BRepFeat.cxx L57).
pub(crate) const NECHANTBARYC: i32 = 11;

/// OCCT BRep_Tool::Surface(fac) (architecture difference #3).
fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::IsSame(S).
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT TopAbs::Reverse (TopAbs.hxx) — FORWARD<->REVERSED, INTERNAL/
/// EXTERNAL unchanged.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// OCCT NCollection_Map::Add — returns true when newly added.
fn map_add(m: &mut OcctShapeMap, s: &Shape) -> bool {
    m.add(shape_key(s), s.clone())
}

/// OCCT BRep_Tool::Curve(edg, Loc, f, l) — the 3D curve and parameter range
/// of an edge (architecture difference #3: the location transformation is
/// the identity-location reduction).
fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> glam::DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// OCCT TopExp::MapShapes over a BRep result pool: the pool is the
/// TopoDS_Compound root; its top-level (unreferenced) TShapes are the
/// compound children. Re-hosted from brep_feat_builder.rs (pool_top_shapes —
/// private there) and bop/brep_algo_api::brep_top_shapes.
pub(crate) fn pool_top_shapes(brep: &BRep, t: ShapeType) -> Vec<Shape> {
    let mut referenced = vec![false; brep.tshapes.len()];
    fn mark(sr: &Shape, referenced: &mut Vec<bool>) {
        let i = sr.index;
        if i >= referenced.len() || referenced[i] {
            return;
        }
        referenced[i] = true;
        match &*sr.data {
            TShape::Solid(sd) => {
                for x in &sd.shells {
                    mark(x, referenced);
                }
                for x in &sd.internal_vertices {
                    mark(x, referenced);
                }
                for x in &sd.internal_edges {
                    mark(x, referenced);
                }
            }
            TShape::Shell(sd) => {
                for x in &sd.faces {
                    mark(x, referenced);
                }
            }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, referenced);
                for w in &fd.inner_wires {
                    mark(w, referenced);
                }
                for v in &fd.internal_vertices {
                    mark(v, referenced);
                }
            }
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    mark(e, referenced);
                }
            }
            TShape::Edge(ed) => {
                mark(&ed.first, referenced);
                mark(&ed.last, referenced);
            }
            TShape::CompSolid(cs) => {
                for x in cs {
                    mark(x, referenced);
                }
            }
            TShape::Compound(cd) => {
                for x in cd {
                    mark(x, referenced);
                }
            }
            _ => {}
        }
    }
    for ts in &brep.tshapes {
        match ts.as_ref() {
            TShape::Solid(sd) => {
                for sr in &sd.shells {
                    mark(sr, &mut referenced);
                }
                for sr in &sd.internal_vertices {
                    mark(sr, &mut referenced);
                }
                for sr in &sd.internal_edges {
                    mark(sr, &mut referenced);
                }
            }
            TShape::Shell(sd) => {
                for sr in &sd.faces {
                    mark(sr, &mut referenced);
                }
            }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, &mut referenced);
                for w in &fd.inner_wires {
                    mark(w, &mut referenced);
                }
                for v in &fd.internal_vertices {
                    mark(v, &mut referenced);
                }
            }
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    mark(e, &mut referenced);
                }
            }
            TShape::Edge(ed) => {
                mark(&ed.first, &mut referenced);
                mark(&ed.last, &mut referenced);
            }
            TShape::CompSolid(shapes) => {
                for sr in shapes {
                    mark(sr, &mut referenced);
                }
            }
            TShape::Compound(shapes) => {
                for sr in shapes {
                    mark(sr, &mut referenced);
                }
            }
            _ => {}
        }
    }
    fn type_of(ts: &TShape) -> ShapeType {
        match ts {
            TShape::Vertex(_) => ShapeType::Vertex,
            TShape::Edge(_) => ShapeType::Edge,
            TShape::Wire(_) => ShapeType::Wire,
            TShape::Face(_) => ShapeType::Face,
            TShape::Shell(_) => ShapeType::Shell,
            TShape::Solid(_) => ShapeType::Solid,
            TShape::CompSolid(_) => ShapeType::CompSolid,
            TShape::Compound(_) => ShapeType::Compound,
        }
    }
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(i, ts)| !referenced[*i] && type_of(ts) == t)
        .map(|(i, ts)| {
            Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// BRepFeat statics (re-hosted from BRepFeat.cxx — architecture difference 3).
// ---------------------------------------------------------------------------

/// OCCT BRepFeat::ParametricMinMax (BRepFeat.cxx L202-333) — computes the
/// parameter range of the shape S on the curve CC. The OCCT out-parameters
/// (prmin, prmax, prbmin, prbmax, flag) map to the returned tuple; theOri
/// keeps the OCCT name.
pub fn brep_feat_parametric_min_max(
    the_s: &Shape,
    the_cc: &Curve3,
    the_ori: bool,
) -> (f64, f64, f64, f64, bool) {
    // OCCT L211-214: LocOpe_CSIntersector ASI(S); scur = {CC}; ASI.Perform.
    let mut a_si = LocOpeCSIntersector::new();
    a_si.init(the_s);
    let scur = vec![Some(the_cc.clone())];
    a_si.perform_cur(&scur);
    // OCCT L204-209 declares the out-parameters uninitialized; both branches
    // below assign them.
    let prmin: f64;
    let prmax: f64;
    let flag: bool;
    // OCCT L215-243.
    if a_si.is_done() && a_si.nb_points(1) >= 1 {
        if !the_ori {
            prmin = a_si
                .point(1, 1)
                .parameter()
                .min(a_si.point(1, a_si.nb_points(1)).parameter());
            prmax = a_si
                .point(1, 1)
                .parameter()
                .max(a_si.point(1, a_si.nb_points(1)).parameter());
        } else {
            let ori = a_si.point(1, 1).orientation();
            if ori == Orientation::Forward {
                prmin = a_si.point(1, 1).parameter();
                prmax = a_si.point(1, a_si.nb_points(1)).parameter();
            } else {
                prmax = a_si.point(1, 1).parameter();
                prmin = a_si.point(1, a_si.nb_points(1)).parameter();
            }
        }
        flag = true;
    } else {
        // OCCT L240-242: prmax = RealFirst(); prmin = RealLast().
        prmax = f64::MIN;
        prmin = f64::MAX;
        flag = false;
    }
    //
    // OCCT L245-254: theMap + the per-edge sampler (the Initialize/Perform
    // pair of Extrema_ExtPC collapses into the rcad constructor, which
    // carries no separate Initialize).
    let mut the_map = OcctShapeMap::new();
    let mut prbmin = f64::MAX;
    let mut prbmax = f64::MIN;
    for edg in explorer(the_s, ShapeType::Edge, ShapeType::Shape) {
        if !map_add(&mut the_map, &edg) {
            continue;
        }
        if !brep_tool_degenerated(&edg) {
            let Some((c, f, l)) = brep_tool_curve(&edg) else {
                continue;
            };
            for i in 1..NECHANTBARYC {
                let prm = ((NECHANTBARYC - i) as f64 * f + i as f64 * l) / NECHANTBARYC as f64;
                let pone = c.point_at(prm);
                // OCCT L273: extpc.Perform(pone).
                let extpc = rcad_kernel::base::extrema::ExtPC::new(
                    pone,
                    the_cc,
                    CONFUSION,
                    the_cc.default_domain()[0],
                    the_cc.default_domain()[1],
                );
                if extpc.is_done() && extpc.nb_ext() >= 1 {
                    let mut dist2_min = extpc.square_distance(1);
                    let mut kmin = 1usize;
                    for k in 2..=extpc.nb_ext() {
                        let dist2 = extpc.square_distance(k);
                        if dist2 < dist2_min {
                            dist2_min = dist2;
                            kmin = k;
                        }
                    }
                    let prmp = extpc.point(kmin).param;
                    if prmp <= prbmin {
                        prbmin = prmp;
                    }
                    if prmp >= prbmax {
                        prbmax = prmp;
                    }
                }
            }
        }
    }
    // OCCT L301-332: adds every vertex.
    for vtx in explorer(the_s, ShapeType::Vertex, ShapeType::Shape) {
        if !map_add(&mut the_map, &vtx) {
            continue;
        }
        let pone = brep_tool_pnt(&vtx);
        let extpc = rcad_kernel::base::extrema::ExtPC::new(
            pone,
            the_cc,
            CONFUSION,
            the_cc.default_domain()[0],
            the_cc.default_domain()[1],
        );
        if extpc.is_done() && extpc.nb_ext() >= 1 {
            let mut dist2_min = extpc.square_distance(1);
            let mut kmin = 1usize;
            for k in 2..=extpc.nb_ext() {
                let dist2 = extpc.square_distance(k);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                    kmin = k;
                }
            }
            let prmp = extpc.point(kmin).param;
            if prmp <= prbmin {
                prbmin = prmp;
            }
            if prmp >= prbmax {
                prbmax = prmp;
            }
        }
    }
    (prmin, prmax, prbmin, prbmax, flag)
}

/// OCCT BRepFeat::ParametricBarycenter (BRepFeat.cxx L119-195) — the
/// "parametric" barycentre of the shape S on the curve CC (consumed by the
/// SensOfPrism/SensOfRevol statics of the form-feature subclasses, 3b).
pub fn brep_feat_parametric_barycenter(the_s: &Shape, the_cc: &Curve3) -> f64 {
    // OCCT L121-130: theMap; extpc over [FirstParameter, LastParameter].
    let mut the_map = OcctShapeMap::new();
    let mut nbp = 0i32;
    let mut parbar = 0f64;
    for edg in explorer(the_s, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L134-138.
        if !map_add(&mut the_map, &edg) {
            continue;
        }
        if !brep_tool_degenerated(&edg) {
            // OCCT L141-142: BRep_Tool::Curve + Transformed (the
            // identity-location reduction).
            let Some((c, f, l)) = brep_tool_curve(&edg) else {
                continue;
            };
            for i in 1..NECHANTBARYC {
                let prm = ((NECHANTBARYC - i) as f64 * f + i as f64 * l) / NECHANTBARYC as f64;
                let pone = c.point_at(prm);
                // OCCT L148: extpc.Perform(pone) — projection on CC.
                let extpc = rcad_kernel::base::extrema::ExtPC::new(
                    pone,
                    the_cc,
                    CONFUSION,
                    the_cc.default_domain()[0],
                    the_cc.default_domain()[1],
                );
                if extpc.is_done() && extpc.nb_ext() >= 1 {
                    let mut dist2_min = extpc.square_distance(1);
                    let mut kmin = 1usize;
                    for k in 2..=extpc.nb_ext() {
                        let dist2 = extpc.square_distance(k);
                        if dist2 < dist2_min {
                            dist2_min = dist2;
                            kmin = k;
                        }
                    }
                    nbp += 1;
                    let prmp = extpc.point(kmin).param;
                    parbar += prmp;
                }
            }
        }
    }
    // OCCT L169-191: adds every vertex (the k-loop of the OCCT source only
    // refreshes Dist2Min; nbp increments regardless — source quirk kept).
    for vtx in explorer(the_s, ShapeType::Vertex, ShapeType::Shape) {
        if !map_add(&mut the_map, &vtx) {
            continue;
        }
        let pone = brep_tool_pnt(&vtx);
        let extpc = rcad_kernel::base::extrema::ExtPC::new(
            pone,
            the_cc,
            CONFUSION,
            the_cc.default_domain()[0],
            the_cc.default_domain()[1],
        );
        if extpc.is_done() && extpc.nb_ext() >= 1 {
            let mut dist2_min = extpc.square_distance(1);
            for k in 2..=extpc.nb_ext() {
                let dist2 = extpc.square_distance(k);
                if dist2 < dist2_min {
                    dist2_min = dist2;
                }
            }
            nbp += 1;
        }
    }
    // OCCT L193-194.
    parbar /= nbp as f64;
    parbar
}

/// OCCT BRepFeat::IsInside(F1, F2) (BRepFeat.cxx L467-520) — GAP
/// (architecture difference #4): the body needs BRepTopAdaptor_FClass2d +
/// GCPnts_QuasiUniformDeflection + GeomProjLib::Curve2d over
/// BRepTools::UVBounds; the classifier is not translated yet. The OCCT
/// structure: every edge of F1 is sampled (QuasiUniformDeflection at
/// 100*Confusion) and classified on F2 (forward-oriented); a single OUT
/// sample rejects.
pub fn brep_feat_is_inside(_the_f1: &Shape, _the_f2: &Shape) -> bool {
    panic!("GAP(BRepFeat_Form): BRepFeat::IsInside needs BRepTopAdaptor_FClass2d + GCPnts_QuasiUniformDeflection (pending translation)");
}

/// OCCT BRepAlgo::IsValid(S) (BRepAlgo_1.cxx L39-43) — GAP (architecture
/// difference #5): BRepCheck_Analyzer on a standalone rcad Shape is pending
/// (rcad's brep_check works on a BRep pool).
pub(crate) fn brep_algo_is_valid(_the_s: &Shape) -> bool {
    panic!("GAP(BRepFeat_Form): BRepAlgo::IsValid needs BRepCheck_Analyzer on a standalone Shape (pending translation)");
}

/// OCCT BRepFeat::FaceUntil(Sbase, FUntil) (BRepFeat.cxx L524-638) —
/// enlarges the face FUntil to cover the bounding box of Sbase ("Limitation
/// of the shape until the case of infinite faces"). The face is rebuilt
/// through the rcad BRepBuilder pool (architecture difference #9).
pub(crate) fn brep_feat_face_until(the_sbase: &Shape, the_f_until: &mut Shape) {
    // OCCT L526-530: Bnd_Box B; BRepBndLib::Add(Sbase, B); diam.
    let Some((pmin, pmax)) = BRepFeatBuilder::shape_box(the_sbase, &[]) else {
        // OCCT: the box always exists; rcad reports None only for shapes
        // without geometry — the OCCT Get would return an empty box and the
        // diam computation degenerates; keep the null-face outcome.
        *the_f_until = Shape::null();
        return;
    };
    let x = [pmin.x, pmax.x];
    let y = [pmin.y, pmax.y];
    let z = [pmin.z, pmax.z];
    let sq_extent = (pmax.x - pmin.x).powi(2)
        + (pmax.y - pmin.y).powi(2)
        + (pmax.z - pmin.z).powi(2);
    let diam = 10.0 * sq_extent.sqrt();
    //
    // OCCT L532-538: the surface of FUntil, basis of a trimmed surface.
    let Some(mut s) = brep_tool_surface(the_f_until) else {
        *the_f_until = Shape::null();
        return;
    };
    if let Surface3::Trimmed(t) = &s {
        s = (*t.basis).clone();
    }
    //
    // OCCT L539-635.
    let str_opt: Option<TrimmedSurface> = match &s {
        Surface3::Plane(pln) => {
            // OCCT L542-576: the u/v box of the 8 box corners on the plane.
            let mut umin = f64::MAX;
            let mut umax = f64::MIN;
            let mut vmin = f64::MAX;
            let mut vmax = f64::MIN;
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let a_p = glam::DVec3::new(x[i], y[j], z[k]);
                        let (u, v) = rcad_kernel::math::el::elslib_plane_parameters(
                            a_p, pln.origin, pln.u_dir, pln.v_dir,
                        );
                        if u < umin {
                            umin = u;
                        }
                        if u > umax {
                            umax = u;
                        }
                        if v < vmin {
                            vmin = v;
                        }
                        if v > vmax {
                            vmax = v;
                        }
                    }
                }
            }
            umin -= diam;
            umax += diam;
            vmin -= diam;
            vmax += diam;
            Some(TrimmedSurface::new(s.clone(), umin, umax, vmin, vmax))
        }
        Surface3::Cylinder(a_cyl) => {
            // OCCT L577-603: the v range of the box corners on the cylinder.
            let mut vmin = f64::MAX;
            let mut vmax = f64::MIN;
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let a_p = glam::DVec3::new(x[i], y[j], z[k]);
                        let x_dir = a_cyl.ref_dir;
                        let y_dir = a_cyl
                            .y_dir
                            .unwrap_or_else(|| a_cyl.axis.cross(x_dir));
                        let (u, v) = rcad_kernel::math::el::elslib_cylinder_parameters(
                            a_p,
                            a_cyl.origin,
                            x_dir,
                            y_dir,
                            a_cyl.axis,
                            a_cyl.radius,
                        );
                        let _ = u;
                        if v < vmin {
                            vmin = v;
                        }
                        if v > vmax {
                            vmax = v;
                        }
                    }
                }
            }
            vmin -= diam;
            vmax += diam;
            let dom = rcad_kernel::SurfaceEval::default_domain(&s);
            Some(TrimmedSurface::new(s.clone(), dom[0], dom[1], vmin, vmax))
        }
        Surface3::Cone(a_con) => {
            // OCCT L604-630: the v range of the box corners on the cone.
            let mut vmin = f64::MAX;
            let mut vmax = f64::MIN;
            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let a_p = glam::DVec3::new(x[i], y[j], z[k]);
                        let x_dir = a_con.ref_dir;
                        let y_dir = a_con.axis.cross(a_con.ref_dir);
                        let (u, v) = rcad_kernel::math::el::elslib_cone_parameters(
                            a_p,
                            a_con.apex,
                            x_dir,
                            y_dir,
                            a_con.axis,
                            a_con.radius,
                            a_con.half_angle_rad,
                        );
                        let _ = u;
                        if v < vmin {
                            vmin = v;
                        }
                        if v > vmax {
                            vmax = v;
                        }
                    }
                }
            }
            vmin -= diam;
            vmax += diam;
            let dom = rcad_kernel::SurfaceEval::default_domain(&s);
            Some(TrimmedSurface::new(s.clone(), dom[0], dom[1], vmin, vmax))
        }
        _ => {
            // OCCT L631-635: FUntil.Nullify(); return.
            *the_f_until = Shape::null();
            return;
        }
    };
    // OCCT L637: FUntil = BRepLib_MakeFace(str, Precision::Confusion()).
    if let Some(str) = str_opt {
        let mut pool = BRep::new();
        let mut b = BRepBuilder::new();
        *the_f_until = b.make_face(&mut pool, Some(Surface3::Trimmed(str)), Shape::null());
    }
}

/// OCCT BRepFeat::Tool(SRef, Fac, Orf) (BRepFeat.cxx L642-712) — builds the
/// tool solid from the faces of SRef, oriented by Fac/Orf. Returns None for
/// the OCCT null solid.
pub(crate) fn brep_feat_tool(the_s_ref: &Shape, the_fac: &Shape, the_orf: Orientation) -> Option<Shape> {
    // OCCT L646-655: lfaces = the faces of SRef.
    let lfaces = explorer(the_s_ref, ShapeType::Face, ShapeType::Shape);
    // OCCT L657-658: LocOpe_BuildShape bs(lfaces); Res = bs.Shape().
    let bs = LocOpeBuildShape::with_faces(&lfaces);
    let Some(res) = bs.shape().cloned() else {
        return None;
    };
    // OCCT L659-674: Sh extraction.
    let mut sh: Option<Shape> = None;
    if res.shape_type() == ShapeType::Shell {
        sh = Some(res.clone());
    } else if res.shape_type() == ShapeType::Solid {
        let mut exp = explorer(&res, ShapeType::Shell, ShapeType::Shape).into_iter();
        sh = exp.next();
        if exp.next().is_some() {
            sh = None;
        }
    }
    let Some(mut sh) = sh else {
        // OCCT L676-680: return the null solid.
        return None;
    };
    // OCCT L682: Sh.Orientation(TopAbs_FORWARD).
    sh.orientation = Orientation::Forward;
    //
    // OCCT L684-693.
    let mut orient = Orientation::Forward;
    for cur in explorer(&sh, ShapeType::Face, ShapeType::Shape) {
        if shape_is_same(&cur, the_fac) {
            orient = cur.orientation;
            break;
        }
    }
    //
    // OCCT L695-705.
    let reverse = (orient == the_fac.orientation && the_orf == Orientation::Reversed)
        || (orient != the_fac.orientation && the_orf == Orientation::Forward);
    if reverse {
        sh.orientation = top_abs_reverse(sh.orientation);
    }
    //
    // OCCT L707-711: B.MakeSolid(Soc); B.Add(Soc, Sh).
    let mut pool = BRep::new();
    let mut b = BRepBuilder::new();
    let soc = b.make_solid(&mut pool, vec![sh]);
    Some(soc)
}

/// OCCT static Descendants(S, theFB, mapF) (BRepFeat_Form.cxx L1490-1508).
/// The theFB.Modified reads (architecture difference #6) return the empty
/// list — the OCCT myFillHistory=false fallback of
/// BOPAlgo_BuilderShape::Modified — so mapF stays empty exactly as in that
/// fallback.
/// OCCT BRepAlgoAPI_Cut re-host (architecture difference #2) — the
/// (PaveFiller, Builder) vehicle with the history knob on; exposes the
/// BRepAlgoAPI surface consumed by BRepFeat_Form::GlobalPerform
/// (Shape/Modified/IsDeleted; BRepAlgoAPI_BuilderAlgo.cxx L203-235).
pub(crate) struct CutVehicle {
    // BRepAlgoAPI_BuilderShape::myShape — the root result shape.
    my_shape: Option<Shape>,
    // BRepAlgoAPI_Algo::myHistory.
    my_history: BRepToolsHistory,
}

impl CutVehicle {
    /// OCCT BRepAlgoAPI_Cut(S1, S2) + Perform (BRepAlgoAPI_Cut.cxx L44-66:
    /// myBuilder myOperation = BOPAlgo_CUT; Perform with the two arguments).
    pub(crate) fn new(the_s1: &Shape, the_s2: &Shape) -> Self {
        Self::with_operation(the_s1, the_s2, BooleanOpType::Cut)
    }

    /// OCCT BRepAlgoAPI_Fuse/Cut(S1, S2) with an explicit operation — the
    /// subclass constructors (BRepAlgoAPI_Fuse.cxx / BRepAlgoAPI_Cut.cxx
    /// L44-66) differ only in myOperation; consumed by the form-feature
    /// subclasses (stage 3b).
    pub(crate) fn with_operation(the_s1: &Shape, the_s2: &Shape, op: BooleanOpType) -> Self {
        let mut vehicle = CutVehicle {
            my_shape: None,
            my_history: BRepToolsHistory::new(),
        };
        // BOPAlgo_BOP::Perform on the rcad vehicle (brep_algo_api::run_build
        // model): a fresh PaveFiller over [S1, S2], a method-scoped Builder.
        let mut p_pf = PaveFiller::new();
        p_pf.set_arguments(vec![the_s1.clone(), the_s2.clone()]);
        let a_prog = rcad_kernel::message::NoopProgress;
        let a_ps = rcad_kernel::message::ProgressScope::new(&a_prog, "intersect", 100);
        p_pf.perform(&a_ps);
        let mut a_builder = crate::bop::algo::builder::Builder::new(
            p_pf.ds(),
            op,
            p_pf.fuzzy_value(),
        );
        a_builder.my_arguments = p_pf.ds().arguments.clone();
        a_builder.my_tools = a_builder.my_arguments[1..].to_vec();
        // OCCT: BRepAlgoAPI_Algo::SetFillHistory(true) — the API default;
        // rcad's Builder defaults to false, set the knob explicitly.
        a_builder.my_fill_history = true;
        match a_builder.build_with_history_topods() {
            Ok((brep, _)) => {
                // The root result shape (the run_build convention: the last
                // Solid/Shell TShape of the pool).
                let root = brep
                    .tshapes
                    .iter()
                    .enumerate()
                    .rev()
                    .find(|(_, ts)| {
                        matches!(
                            ts.as_ref(),
                            TShape::Solid(_) | TShape::Shell(_)
                        )
                    })
                    .map(|(i, ts)| {
                        Shape::from_parts(
                            ts.clone(),
                            i,
                            0,
                            rcad_kernel::topods::Orientation::Forward,
                        )
                    });
                vehicle.my_shape = root;
            }
            Err(_) => {
                vehicle.my_shape = None;
            }
        }
        vehicle.my_history = a_builder
            .my_history
            .take()
            .unwrap_or_else(BRepToolsHistory::new);
        vehicle
    }

    /// OCCT BRepAlgoAPI_Cut default ctor + SetArguments(aLObj) +
    /// SetTools(aLTools) + Build() — the N-arguments/N-tools form consumed
    /// by the form-feature subclasses (stage 3b). The PaveFiller receives
    /// the concatenated list; the Builder keeps the tools tail (the
    /// run_build vehicle model).
    pub(crate) fn with_args_tools(
        arguments: Vec<Shape>,
        tools: Vec<Shape>,
        op: BooleanOpType,
    ) -> Self {
        let mut vehicle = CutVehicle {
            my_shape: None,
            my_history: BRepToolsHistory::new(),
        };
        // BOPAlgo_BOP with N arguments + N tools on the rcad vehicle.
        let mut all_args = arguments;
        all_args.extend(tools.iter().cloned());
        let mut p_pf = PaveFiller::new();
        p_pf.set_arguments(all_args);
        let a_prog = rcad_kernel::message::NoopProgress;
        let a_ps = rcad_kernel::message::ProgressScope::new(&a_prog, "intersect", 100);
        p_pf.perform(&a_ps);
        let mut a_builder = crate::bop::algo::builder::Builder::new(
            p_pf.ds(),
            op,
            p_pf.fuzzy_value(),
        );
        a_builder.my_arguments = p_pf.ds().arguments.clone();
        a_builder.my_tools = tools;
        // OCCT: BRepAlgoAPI_Algo::SetFillHistory(true) — the API default.
        a_builder.my_fill_history = true;
        match a_builder.build_with_history_topods() {
            Ok((brep, _)) => {
                let root = brep
                    .tshapes
                    .iter()
                    .enumerate()
                    .rev()
                    .find(|(_, ts)| {
                        matches!(
                            ts.as_ref(),
                            TShape::Solid(_) | TShape::Shell(_)
                        )
                    })
                    .map(|(i, ts)| {
                        Shape::from_parts(
                            ts.clone(),
                            i,
                            0,
                            rcad_kernel::topods::Orientation::Forward,
                        )
                    });
                vehicle.my_shape = root;
            }
            Err(_) => {
                vehicle.my_shape = None;
            }
        }
        vehicle.my_history = a_builder
            .my_history
            .take()
            .unwrap_or_else(BRepToolsHistory::new);
        vehicle
    }

    /// OCCT BRepAlgoAPI_BuilderShape::Shape().
    pub(crate) fn shape(&self) -> Option<&Shape> {
        self.my_shape.as_ref()
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::Modified(theS) (BuilderAlgo.cxx
    /// L203-210) — the history list; returned by value (rcad history API).
    pub(crate) fn modified(&self, the_s: &Shape) -> Vec<Shape> {
        self.my_history.modified(the_s)
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::IsDeleted(theS) (BuilderAlgo.cxx
    /// L228-231) — the history IsRemoved.
    pub(crate) fn is_deleted(&self, the_s: &Shape) -> bool {
        self.my_history.is_removed(the_s)
    }
}

/// OCCT BRepFeat_Builder::Shape() — the result shape of the builder
/// (architecture difference #9: the rcad result is a BRep pool; the root
/// compound is rebuilt over its top-level shapes).
pub(crate) fn builder_result_shape(the_builder: &mut BRepFeatBuilder) -> Option<Shape> {
    let brep = the_builder.shape_brep()?;
    let tops = pool_top_shapes(&brep, ShapeType::Shape);
    let mut pool = BRep::new();
    let mut b = BRepBuilder::new();
    let comp = b.make_compound(&mut pool, tops);
    Some(comp)
}

