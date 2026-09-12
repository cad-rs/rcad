// OCCT BRepFeat.cxx statics + vehicles consumed by BRepFeat_Form — 1:1
// translation (split from brep_feat_form.rs per the 2000-line rule; the OCCT
// line anchors are the same files as the parent module header).
//
// Content: the re-hosted BRepFeat package statics (ParametricMinMax /
// IsInside / FaceUntil / Tool — BRepFeat.cxx), the BRepAlgo::IsValid gap
// marker (BRepAlgo_1.cxx L39-43), the BRepAlgoAPI_Cut vehicle (CutVehicle)
// and the result-pool helpers (pool_top_shapes / builder_result_shape).
// IsInside is the full 1:1 body (BRepFeat.cxx L337-520: IsIn /
// PutInBoundsU / PutInBoundsV / IsInside over the topalgo FClass2d +
// GCPnts_QuasiUniformDeflection + GeomProjLib::Curve2d bodies).
// The architecture-difference numbering of the parent module header applies.

use crate::bop::algo::builder::BooleanOpType;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::history::BRepToolsHistory;
use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder, OcctShapeMap};
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::Curve2dEval;
use rcad_kernel::geom::{Curve2d, Curve3, CurveEval, Surface3, SurfaceEval, TrimmedSurface};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, State};

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
    // OCCT L245-254: theMap; GeomAdaptor_Curve TheCurve(CC); the default
    // Extrema_ExtPC and its Initialize(TheCurve, CC->FirstParameter(),
    // CC->LastParameter()) — the default theTolF is 1.0e-10.
    let mut the_map = OcctShapeMap::new();
    let a_adaptor = GeomCurveAdaptor::new(the_cc.clone());
    let a_tool = CurveToolHandle::for_curve3(the_cc, &a_adaptor, &a_adaptor);
    let mut extpc = ExtremaExtPC::new();
    extpc.initialize(
        &a_tool,
        the_cc.default_domain()[0],
        the_cc.default_domain()[1],
        1.0e-10,
    );
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
                // OCCT L271-273: gp_Pnt pone = C->Value(prm);
                // extpc.Perform(pone) — the projection on CC.
                extpc.perform(pone);
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
        // OCCT L309-311: extpc.Perform(pone) — the projection on CC.
        extpc.perform(pone);
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
    // OCCT L121-130: theMap; GeomAdaptor_Curve TheCurve(CC); the default
    // Extrema_ExtPC and its Initialize(TheCurve, CC->FirstParameter(),
    // CC->LastParameter()) — the default theTolF is 1.0e-10.
    let mut the_map = OcctShapeMap::new();
    let mut nbp = 0i32;
    let mut parbar = 0f64;
    let a_adaptor = GeomCurveAdaptor::new(the_cc.clone());
    let a_tool = CurveToolHandle::for_curve3(the_cc, &a_adaptor, &a_adaptor);
    let mut extpc = ExtremaExtPC::new();
    extpc.initialize(
        &a_tool,
        the_cc.default_domain()[0],
        the_cc.default_domain()[1],
        1.0e-10,
    );
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
                // OCCT L146-148: gp_Pnt pone = C->Value(prm);
                // extpc.Perform(pone) — the projection on CC.
                extpc.perform(pone);
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
        // OCCT L183-185: extpc.Perform(pone) — the projection on CC.
        extpc.perform(pone);
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

/// OCCT Geom_Surface::UPeriod() — the U period of the elementary periodic
/// surfaces (pure-math re-host over the rcad surface value).
fn surface_u_period(s: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match s {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => TAU,
        _ => 0.0,
    }
}

/// OCCT Geom_SphericalSurface/ToroidalSurface::VPeriod (2*PI sphere,
/// 2*minorRadius torus; pure-math re-host).
fn surface_v_period(s: &Surface3) -> f64 {
    use std::f64::consts::TAU;
    match s {
        Surface3::Sphere(_) => TAU,
        Surface3::Torus(t) => 2.0 * t.minor_radius,
        _ => 0.0,
    }
}

/// OCCT static IsIn(FC, AC) (BRepFeat.cxx L337-351): every deflection sample
/// of the pcurve is classified on F2; a single OUT sample rejects.
fn brep_feat_is_in(
    fc: &crate::topalgo::brep_top_adaptor::fclass2d::FClass2d,
    src: &crate::topalgo::shape_source::FaceShapeSource,
    c2d: &Curve2d,
    f1: f64,
    l1: f64,
) -> bool {
    // OCCT L339: Def = 100 * Precision::Confusion().
    let def = 100.0 * rcad_kernel::precision::CONFUSION;
    // OCCT L340: GCPnts_QuasiUniformDeflection QU(AC, Def) — the adaptor
    // carries the (f1, l1) range.
    let q_u = crate::topalgo::gcpnts::QuasiUniformDeflection::new(c2d, def, f1, l1);
    for i in 1..=q_u.nb_points() {
        // OCCT L344: P = AC.Value(QU.Parameter(i)).
        let p = c2d.point_at(q_u.parameter(i));
        // OCCT L345: FC.Perform(P, false) == TopAbs_OUT.
        if fc.perform(src, p, false) == State::Out {
            return false;
        }
    }
    true
}

/// OCCT static PutInBoundsU (BRepFeat.cxx L361-407) — recadre la courbe 2d
/// dans les bounds de la face (U direction).
fn put_in_bounds_u(
    umin: f64,
    umax: f64,
    eps: f64,
    period: f64,
    f: f64,
    l: f64,
    c2d: &mut Curve2d,
) {
    let pf = c2d.point_at(f);
    let pl = c2d.point_at(l);
    let pm = c2d.point_at(0.34 * f + 0.66 * l);
    let mut min_c = pf.x.min(pl.x);
    min_c = min_c.min(pm.x);
    let mut max_c = pf.x.max(pl.x);
    max_c = max_c.max(pm.x);
    let mut du = 0.0;
    if min_c < umin - eps {
        du = (((umin - min_c) / period) as i32 as f64 + 1.0) * period;
    }
    if min_c > umax + eps {
        du = -(((min_c - umax) / period) as i32 as f64 + 1.0) * period;
    }
    if du != 0.0 {
        // OCCT: C2d->Translate(gp_Vec2d(du, 0.)).
        *c2d = rcad_kernel::geom::translate_curve2d(c2d, glam::DVec2::new(du, 0.0));
        min_c += du;
        max_c += du;
    }
    // Ajuste au mieux la courbe dans le domaine.
    if max_c > umax + 100.0 * eps {
        let d1 = max_c - umax;
        let d2 = umin - min_c + period;
        if d2 < d1 {
            du = -period;
        }
        if du != 0.0 {
            *c2d = rcad_kernel::geom::translate_curve2d(c2d, glam::DVec2::new(du, 0.0));
        }
    }
}

/// OCCT static PutInBoundsV (BRepFeat.cxx L417-463) — the V-direction twin
/// of PutInBoundsU.
fn put_in_bounds_v(
    vmin: f64,
    vmax: f64,
    eps: f64,
    period: f64,
    f: f64,
    l: f64,
    c2d: &mut Curve2d,
) {
    let pf = c2d.point_at(f);
    let pl = c2d.point_at(l);
    let pm = c2d.point_at(0.34 * f + 0.66 * l);
    let mut min_c = pf.y.min(pl.y);
    min_c = min_c.min(pm.y);
    let mut max_c = pf.y.max(pl.y);
    max_c = max_c.max(pm.y);
    let mut dv = 0.0;
    if min_c < vmin - eps {
        dv = (((vmin - min_c) / period) as i32 as f64 + 1.0) * period;
    }
    if min_c > vmax + eps {
        dv = -(((min_c - vmax) / period) as i32 as f64 + 1.0) * period;
    }
    if dv != 0.0 {
        // OCCT: C2d->Translate(gp_Vec2d(0., dv)).
        *c2d = rcad_kernel::geom::translate_curve2d(c2d, glam::DVec2::new(0.0, dv));
        min_c += dv;
        max_c += dv;
    }
    // Ajuste au mieux la courbe dans le domaine.
    if max_c > vmax + 100.0 * eps {
        let d1 = max_c - vmax;
        let d2 = vmin - min_c + period;
        if d2 < d1 {
            dv = -period;
        }
        if dv != 0.0 {
            *c2d = rcad_kernel::geom::translate_curve2d(c2d, glam::DVec2::new(0.0, dv));
        }
    }
}

/// OCCT BRepFeat::IsInside(F1, F2) (BRepFeat.cxx L467-520): every edge of F1
/// is projected on F2's surface (GeomProjLib::Curve2d), re-fitted into the
/// surface bounds for the periodic directions (PutInBoundsU/V) and sampled
/// (QuasiUniformDeflection at 100*Confusion) against the forward-oriented
/// F2 classifier (BRepTopAdaptor_FClass2d); a single OUT sample rejects.
pub fn brep_feat_is_inside(the_f1: &Shape, the_f2: &Shape) -> bool {
    // OCCT L475: S = BRep_Tool::Surface(F2).
    let Some(s) = brep_tool_surface(the_f2) else {
        // OCCT dereferences the null surface handle (S->IsUPeriodic()).
        panic!("Standard_NoSuchObject: BRepFeat::IsInside - null face surface");
    };
    // OCCT L477: BRepTools::UVBounds(F2, umin, umax, vmin, vmax).
    let Some(bounds) = crate::feat::loc_ope_wires_on_shape_b::brep_tools_uv_bounds(the_f2)
    else {
        panic!("Standard_NoSuchObject: BRepFeat::IsInside - empty UV bounds");
    };
    let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
    // OCCT L479-489: the periodic flags.
    let (mut flagu, mut flagv) = (0i32, 0i32);
    let (mut uperiod, mut vperiod) = (0.0f64, 0.0f64);
    if s.is_u_periodic() {
        flagu = 1;
        uperiod = surface_u_period(&s);
    }
    if s.is_v_periodic() {
        flagv = 1;
        vperiod = surface_v_period(&s);
    }
    // OCCT L490-491: BRepTopAdaptor_FClass2d FC(F2.Oriented(FORWARD),
    // Precision::Confusion()).
    let mut f2_forward = the_f2.clone();
    f2_forward.orientation = Orientation::Forward;
    let locations = [glam::DAffine3::IDENTITY];
    let src = crate::topalgo::shape_source::FaceShapeSource::new(
        &f2_forward,
        s.clone(),
        &locations,
    );
    let fc = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
        &src,
        0,
        rcad_kernel::precision::CONFUSION,
    );
    // OCCT L494-518: the F1 edge loop.
    for exp in explorer(the_f1, ShapeType::Edge, ShapeType::Shape) {
        let Some((c0, f1, l1)) = crate::feat::loc_ope_wires_on_shape_b::brep_tool_curve(&exp)
        else {
            // OCCT: a null 3D curve makes GeomProjLib::Curve2d return a null
            // handle and Geom2dAdaptor_Curve AC(null) fails on Value — the
            // rcad curve engine has no null curve value; the edge is skipped
            // (documented arch. diff.).
            continue;
        };
        // OCCT L498: C = GeomProjLib::Curve2d(C0, f1, l1, S).
        let Some(mut c) =
            rcad_kernel::base::geom_proj_lib::curve2d_simple(&c0, f1, l1, &s)
        else {
            // OCCT: the null pcurve handle would fail the adaptor; the rcad
            // Option is the failure outcome — treated as OUT (the sample
            // cannot be classified inside).
            return false;
        };
        // OCCT L500-512: the periodic re-fitting.
        if flagu == 1 || flagv == 1 {
            let eps = crate::feat::loc_ope_wires_on_shape_b::brep_tool_tolerance(&exp);
            // OCCT L503: BRep_Tool::Range(E, f1, l1) — the same range.
            let _ = crate::feat::loc_ope_wires_on_shape_b::brep_tool_range(&exp);
            if flagu == 1 {
                put_in_bounds_u(umin, umax, eps, uperiod, f1, l1, &mut c);
            }
            if flagv == 1 {
                put_in_bounds_v(vmin, vmax, eps, vperiod, f1, l1, &mut c);
            }
        }
        // OCCT L513-517: Geom2dAdaptor_Curve AC(C, f1, l1); IsIn(FC, AC).
        if !brep_feat_is_in(&fc, &src, &c, f1, l1) {
            return false;
        }
    }
    true
}

/// OCCT BRepAlgo::IsValid(S) (BRepAlgo_1.cxx L39-43):
/// `BRepCheck_Analyzer ana(S); return ana.IsValid();` — the default ctor
/// enables GeomControls.
///
/// Architecture difference: OCCT's `TopoDS_Shape` carries its `TShape` graph
/// by pointer, so a shape assembled from several `BRep_Builder` scopes is
/// still one coherent graph. The rcad `BRep` is a flat pool and a `Shape`'s
/// `index` field is only meaningful inside the pool it was created in — the
/// feat callers routinely mix shapes from several pools (the profile pool,
/// the BndFace/CutVehicle result pool, cloned input-wire edges), so equal
/// indices can belong to different TShapes. A scratch pool naively keyed by
/// the source indices therefore aliases unrelated TShapes and makes the
/// analyzer read the wrong shape ("Shape N is not a Face").
///
/// The subgraph is therefore collected into a scratch pool under FRESH
/// contiguous indices, and every copied TShape has its child reference
/// re-pointed at the scratch slot AND the scratch `TShape` handle (a bare
/// index re-point would leave `data` on the source pool and re-open the
/// collision one level deeper). Identity inside the analyzer is still by
/// `TShape` pointer (`ShapeKey::of`), so the renumbering is transparent to
/// every consumer of the scratch pool — the shapes simply all belong to the
/// scratch pool now.
///
/// The walk is a POST-order DFS: a sub-graph is a DAG (a sub-shape can be
/// shared by several parents), so only post-order guarantees that a node is
/// laid down after every child it references.
pub(crate) fn brep_algo_is_valid(the_s: &Shape) -> bool {
    let mut order: Vec<Shape> = Vec::new();
    let mut visited: std::collections::HashSet<u64> = std::collections::HashSet::new();
    collect_subgraph_indexed(the_s, &mut order, &mut visited);
    // ptr_id -> scratch slot.
    let remap: std::collections::HashMap<u64, usize> = order
        .iter()
        .enumerate()
        .map(|(i, s)| (s.ptr_id(), i))
        .collect();
    // Children first (post-order), each node laid down exactly once.
    let mut final_arcs: Vec<Option<std::sync::Arc<TShape>>> = vec![None; order.len()];
    for i in 0..order.len() {
        let mut ts = (*order[i].data).clone();
        remap_tshape_children(&mut ts, &remap, &final_arcs);
        final_arcs[i] = Some(std::sync::Arc::new(ts));
    }
    let root_index = order.len() - 1; // the root closes the post-order walk
    let mut scratch = rcad_kernel::topods::BRep::new();
    scratch.tshapes = final_arcs
        .into_iter()
        .map(|a| a.expect("every slot is built by the post-order pass"))
        .collect();
    let mut root = Shape::from_parts(
        scratch.tshapes[root_index].clone(),
        root_index,
        the_s.location,
        the_s.orientation,
    );
    let ana = crate::topalgo::brep_check::brep_check_analyzer::BRepCheckAnalyzer::new(
        &scratch,
        &root,
        true,
    );
    let ok = ana.is_valid(&scratch, &root);
    ok
}

/// Re-point every child reference of `ts` at its scratch slot (the pool
/// renumbering of `brep_algo_is_valid`). References whose TShape is not in
/// the scratch pool (a null child) are left untouched.
fn remap_tshape_children(
    ts: &mut TShape,
    remap: &std::collections::HashMap<u64, usize>,
    final_arcs: &[Option<std::sync::Arc<TShape>>],
) {
    let fix = |s: &mut Shape| {
        if let Some(&i) = remap.get(&s.ptr_id()) {
            s.index = i;
            s.data = final_arcs[i]
                .as_ref()
                .expect("children are built before their parent")
                .clone();
        }
    };
    match ts {
        TShape::Vertex(vd) => {
            vd.my_shapes.iter_mut().for_each(fix);
        }
        TShape::Edge(ed) => {
            ed.my_shapes.iter_mut().for_each(fix);
            fix(&mut ed.first);
            fix(&mut ed.last);
        }
        TShape::Wire(wd) => {
            wd.my_shapes.iter_mut().for_each(fix);
            wd.edges.iter_mut().for_each(fix);
        }
        TShape::Face(fd) => {
            fd.my_shapes.iter_mut().for_each(fix);
            fix(&mut fd.outer_wire);
            fd.inner_wires.iter_mut().for_each(fix);
            fd.internal_vertices.iter_mut().for_each(fix);
        }
        TShape::Shell(sd) => {
            sd.my_shapes.iter_mut().for_each(fix);
            sd.faces.iter_mut().for_each(fix);
        }
        TShape::Solid(sd) => {
            sd.my_shapes.iter_mut().for_each(fix);
            sd.shells.iter_mut().for_each(fix);
            sd.internal_edges.iter_mut().for_each(fix);
            sd.internal_vertices.iter_mut().for_each(fix);
        }
        TShape::CompSolid(children) | TShape::Compound(children) => {
            children.iter_mut().for_each(fix);
        }
    }
}

/// The subgraph walk of brep_builderapi_transform::collect_subgraph — the
/// POST-order dfs of the distinct sub-shapes (deduplicated by TShape
/// pointer), so a node is always emitted after every child it references.
fn collect_subgraph_indexed(
    s: &Shape,
    out: &mut Vec<Shape>,
    visited: &mut std::collections::HashSet<u64>,
) {
    if s.is_null() || !visited.insert(s.ptr_id()) {
        return;
    }
    match s.data.as_ref() {
        TShape::Vertex(_) => {}
        TShape::Edge(ed) => {
            collect_subgraph_indexed(&ed.first, out, visited);
            collect_subgraph_indexed(&ed.last, out, visited);
        }
        TShape::Wire(wd) => {
            for e in &wd.edges {
                collect_subgraph_indexed(e, out, visited);
            }
        }
        TShape::Face(fd) => {
            collect_subgraph_indexed(&fd.outer_wire, out, visited);
            for w in &fd.inner_wires {
                collect_subgraph_indexed(w, out, visited);
            }
            for v in &fd.internal_vertices {
                collect_subgraph_indexed(v, out, visited);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                collect_subgraph_indexed(f, out, visited);
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                collect_subgraph_indexed(sh, out, visited);
            }
            for e in &sd.internal_edges {
                collect_subgraph_indexed(e, out, visited);
            }
            for v in &sd.internal_vertices {
                collect_subgraph_indexed(v, out, visited);
            }
        }
        TShape::CompSolid(children) | TShape::Compound(children) => {
            for c in children {
                collect_subgraph_indexed(c, out, visited);
            }
        }
    }
    // POST-order: the node is emitted once all of its children are.
    out.push(s.clone());
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
                // Solid/Shell TShape of the pool).  A face-level BOP result
                // (e.g. BRepAlgoAPI_Common(Solid, Face) — BOPAlgo_BOP::
                // BuildShape L1092 result compound) has no container: the
                // compound members are the result faces, carried flat in the
                // rcad pool, so fall back to the last result Face.
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
                    .or_else(|| {
                        brep.tshapes.iter().enumerate().rev().find(|(_, ts)| {
                            matches!(ts.as_ref(), TShape::Face(_))
                        })
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
                // The root result shape — the Solid/Shell convention with the
                // face-level fallback (see with_operation).
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
                    .or_else(|| {
                        brep.tshapes.iter().enumerate().rev().find(|(_, ts)| {
                            matches!(ts.as_ref(), TShape::Face(_))
                        })
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

