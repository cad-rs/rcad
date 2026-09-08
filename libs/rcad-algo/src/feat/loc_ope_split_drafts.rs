// OCCT LocOpe_SplitDrafts.hxx L17-113 + LocOpe_SplitDrafts.cxx L17-1976 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitDrafts.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_SplitDrafts.cxx
//
// OCCT inheritance chain: none (standalone value class).
//
// The four OCCT file-local statics (NewPlane L1461-1538, MakeFace
// L1542-1712, Contains L1716-1727, NewEdge L1731-1976) live in the sibling
// submodule loc_ope_split_drafts_b (single-file <2000-line rule).
//
// Architecture differences (referenced from the affected functions):
// 1. LocOpe_Spliter / LocOpe_WiresOnShape — deferred translations (Stage 3d;
//    the loc_ope_spliter.rs / loc_ope_wires_on_shape.rs modules are still
//    placeholders). The API-surface carriers already defined in
//    feat::brep_feat_split_shape are reused here; until the Stage 3d bodies
//    land, Spls.Perform is a no-op with IsDone=false and Perform returns
//    through the OCCT L367-370 guard.
// 2. LocOpe_SplitShape — not translated yet anywhere in rcad; carried below
//    as the API-surface carrier LocOpeSplitShape (hxx L38-91 member set plus
//    the ctor/Add(V,P,E)/DescendantShapes surface consumed by Perform). The
//    carrier keeps the OCCT L298-308 structure; until its body lands
//    DescendantShapes is unbound and Perform returns at the OCCT
//    `lres.Extent() != 1` guard.
// 3. GeomFill_Pipe (TKGeomAlgo) — not translated (rcad geomfill carries the
//    laws only); carried below as GeomFillPipe with the
//    GenerateParticularCase/Init/Perform/IsDone/Surface surface consumed by
//    Perform. IsDone stays false until the body lands, so the OCCT
//    `throw Standard_ConstructionError("GeomFill_Pipe : ...")` raise
//    structure is preserved.
// 4. GeomInt_IntSS (TKBool) — not translated; the reduced kernel vehicle
//    rcad_kernel base::geom_api::int_ss::IntSS is the GeomAPI_IntSS wrapper
//    and lacks the AppS1/AppS2 Perform overloads and the LineOnS1/LineOnS2
//    p-curves both consumed here. Carried below as GeomIntIntSS over the
//    OCCT member surface (myDone, myLines, myLinesOnS1, myLinesOnS2); the
//    intersection body is deferred (IsDone=false), so Perform takes the
//    OCCT L332 else path (wire not split).
// 5. BRepTools_Substitution — carried below as BRepToolsSubstitution;
//    Substitute/Copy/IsCopied are the real OCCT map operations
//    (BRepTools_Substitution.cxx L30-58); Build(S) (the ancestor rebuild
//    walk) is deferred — it needs the TopoDS rebuild vehicle — and is a
//    no-op until that lands.
// 6. BRep_Builder — rcad's BRepBuilder over a local topods::BRep pool (the
//    LocOpe_BuildShape model). MakeVertex is pool.add_tvertex_unique +
//    tolerance (OCCT MakeVertex never dedups by position). MakeEdge(C) +
//    Add(V)/UpdateVertex(V,P,E) collapse into BRepBuilder::add_edge +
//    set_vertex_param (the OCCT vertex parameters are stored explicitly).
//    MakeFace(F, S, Tol) is pool.add_tface_tol over a fresh marker wire —
//    the kernel face cache keys on the outer wire, and the OCCT face is
//    wire-less. BRep_Builder::Add(F, W) first-wire-becomes-outer semantics
//    are re-hosted as builder_add_face_wire.
// 7. TopExp_Explorer / TopExp::Vertices / TopExp::MapShapesAndAncestors —
//    the re-hosts are TopExpExplorer (over feat::brep_feat_builder::
//    explorer), fillet::chfi3d_builder_0::topexp_vertices and
//    feat::loc_ope_glued_shape::map_shapes_and_ancestors (same crate).
//    TopoDS::Face/Edge/Wire/Vertex casts are identity on rcad Shape
//    handles.
// 8. GeomAdaptor_Curve / GeomAdaptor_Surface — Load(x) is the construction
//    of bop::int_tools::bean_face_intersector::{BRepAdaptorCurve,
//    BRepAdaptorSurface}; IntCurveSurface_HInter is
//    geomalgo::int_curve_surface::HInter driven by the matching Tool
//    markers (bop::int_tools::hinter_adaptor).
// 9. IntAna_QuadQuadGeo(Pln, Pln, tol) — the plane-plane case is
//    rcad_kernel base::int_ana::intersect_plane_plane_intana (PlnPlnResult;
//    TypeInter()==IntAna_Line maps to PlnPlnResult::Line).
// 10. Extrema_ExtPC — rcad_kernel base::extrema::ExtPC; the OCCT default
//     constructor runs over the adaptor parameter range with
//     Precision::Confusion(). TrimmedSquareDistances (the square distances
//     at FirstParameter/LastParameter) is re-hosted below
//     (ext_trimmed_square_distances).
// 11. BRepGProp::SurfaceProperties(NewFace, GP) — GAP: the rcad gprop
//     vehicle computes whole-pool areas; the single-Shape-face overload is
//     deferred (brep_gprop_surface_properties_mass returns 0.0), so the
//     OCCT `GP.Mass() < 0` wire-reversal branch keeps its structure but
//     stays inactive until the BRepGProp translation lands.
// 12. Geom_Curve::Transformed(Loc.Transformation()) — rcad carries the
//     location as an id and the feat flows are identity-location; the
//     re-application of the location transformation is an identity
//     (GAP with the kernel location plumbing).
// 13. TopoDS_Shape null handle — Option/Shape::null() carries it; the OCCT
//     map raises (NoSuchObject/StdFail_NotDone) map to expect()/panic! at
//     the same points.
//
// Draft package note: this file does not consume Draft_Modification (the
// Draft package, 2d batch not started) — SplitDrafts carries its own
// NewPlane/draft-plane math; no Draft gap arises here.

use crate::bop::int_tools::bean_face_intersector::{
    BRepAdaptorCurve, BRepAdaptorSurface, IntCurveSurfaceHInter,
};
use crate::feat::brep_feat_builder::explorer;
use crate::feat::brep_feat_split_shape::{LocOpeSpliter, LocOpeWiresOnShape};
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::fillet::chfi3d_builder_0::{brep_tool_parameter, topexp_vertices};
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::base::extrema::ExtPC;
use rcad_kernel::geom::{Curve2d, Curve3, CurveEval, Line3, Plane, Surface3};
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::{HashMap, HashSet};

use crate::geomalgo::geomfill::corrected_frenet::Trihedron;

use crate::feat::loc_ope_split_drafts_b::{
    contains, make_face, new_edge, new_plane, BRepToolsSubstitution, GeomFillPipe, GeomIntIntSS,
    LocOpeSplitShape,
};

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
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
pub(crate) fn with_orientation(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT TopExp_Explorer (TopExp_Explorer.cxx) — the Init/More/Current/Next/
/// ReInit surface over the rcad explorer walk.
pub(crate) struct TopExpExplorer {
    the_items: Vec<Shape>,
    the_cur: usize,
}

impl TopExpExplorer {
    /// OCCT TopExp_Explorer::Init(S, ToFind, ToAvoid).
    pub(crate) fn init(&mut self, the_s: &Shape, the_to_find: ShapeType, the_to_avoid: ShapeType) {
        self.the_items = explorer(the_s, the_to_find, the_to_avoid);
        self.the_cur = 0;
    }

    /// OCCT TopExp_Explorer::More().
    pub(crate) fn more(&self) -> bool {
        self.the_cur < self.the_items.len()
    }

    /// OCCT TopExp_Explorer::Current().
    pub(crate) fn current(&self) -> Shape {
        self.the_items[self.the_cur].clone()
    }

    /// OCCT TopExp_Explorer::Next().
    pub(crate) fn next(&mut self) {
        self.the_cur += 1;
    }

    /// OCCT TopExp_Explorer::ReInit() — re-explores the same target.
    pub(crate) fn re_init(&mut self) {
        self.the_cur = 0;
    }
}

impl Default for TopExpExplorer {
    fn default() -> Self {
        TopExpExplorer {
            the_items: Vec::new(),
            the_cur: 0,
        }
    }
}

/// OCCT BRep_Tool::Pnt(vtx).
pub(crate) fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(shape) — vertex/edge/face tolerance.
pub(crate) fn brep_tool_tolerance(the_shape: &Shape) -> f64 {
    match the_shape.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        TShape::Face(fd) => fd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Curve(edg, f, l) — the 3D curve and parameter range of
/// the edge (architecture difference #6: the OCCT location overload and the
/// Transformed re-application are an identity in the rcad feat flows —
/// architecture difference #12).
pub(crate) fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Range(edg, f, l).
pub(crate) fn brep_tool_range(edg: &Shape) -> (f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (ed.range[0], ed.range[1]),
        _ => (0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Surface(face).
pub(crate) fn brep_tool_surface(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edg, face, f, l) — the pcurve of the edge
/// on the face with its range.
pub(crate) fn brep_tool_curve_on_surface(edg: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT Geom_TrimmedCurve::BasisCurve() dispatch — the DynamicType ==
/// STANDARD_TYPE(Geom_TrimmedCurve) test of the OCCT source.
pub(crate) fn geom_trimmed_basis_curve(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Trimmed(tc) => (*tc.curve).clone(),
        _ => c.clone(),
    }
}

/// OCCT Geom_RectangularTrimmedSurface::BasisSurface() dispatch — the
/// DynamicType == STANDARD_TYPE(Geom_RectangularTrimmedSurface) test of the
/// OCCT source.
pub(crate) fn geom_rectangular_trimmed_basis_surface(s: &Surface3) -> Surface3 {
    match s {
        Surface3::Trimmed(ts) => (*ts.basis).clone(),
        _ => s.clone(),
    }
}

/// OCCT gp_Pln::Direct() (gp_Ax3::Direct — (X ^ Y) * Z > 0).
pub(crate) fn gp_pln_direct(p: &Plane) -> bool {
    p.u_dir.cross(p.v_dir).dot(p.normal) > 0.0
}

/// OCCT gp_Pln::Axis() — the position axis (location + normal).
pub(crate) fn gp_pln_axis(p: &Plane) -> Ax1 {
    Ax1::new(p.origin, p.normal)
}

/// OCCT gp_Pln::Rotated(A1, Ang) — the plane position rotated about the
/// axis (gp_Trsf::SetRotation; the quaternion form of OCCT gp_Trsf.cxx
/// L100-140 and the Rodrigues form used here are the same rotation up to
/// floating-point evaluation order).
pub(crate) fn gp_pln_rotated(p: &Plane, the_axe: &Ax1, the_theta: f64) -> Plane {
    let (s, c) = the_theta.sin_cos();
    let rot = |v: DVec3| -> DVec3 {
        // Rodrigues: v*cos(T) + (K x v)*sin(T) + K*(K.v)*(1-cos(T)).
        let k = the_axe.direction;
        v * c + k.cross(v) * s + k * (k.dot(v)) * (1.0 - c)
    };
    let d = the_axe.location;
    Plane {
        origin: d + rot(p.origin - d),
        normal: rot(p.normal),
        u_dir: rot(p.u_dir),
        v_dir: rot(p.v_dir),
    }
}

/// OCCT Extrema_ExtPC::TrimmedSquareDistances(DistMin, Dist, PB1, PB2) —
/// the square distances from the constructor point to the curve points of
/// parameter FirstParameter / LastParameter (architecture difference #10).
pub(crate) fn ext_trimmed_square_distances(
    the_p: DVec3,
    the_curve: &Curve3,
) -> (f64, f64, DVec3, DVec3) {
    let dom = the_curve.default_domain();
    let p1 = the_curve.point_at(dom[0]);
    let p2 = the_curve.point_at(dom[1]);
    (
        (the_p - p1).length_squared(),
        (the_p - p2).length_squared(),
        p1,
        p2,
    )
}


/// OCCT BRepGProp::SurfaceProperties(NewFace, GP) over the single-face
/// shape — GAP: the rcad gprop vehicle computes whole-pool areas
/// (architecture difference #11); until the BRepGProp translation lands the
/// mass is 0.0 and the OCCT `GP.Mass() < 0` branch keeps its structure but
/// stays inactive.
pub(crate) fn brep_gprop_surface_properties_mass(_pool: &BRep, _the_face: &Shape) -> f64 {
    0.0
}

/// OCCT BRep_Builder::Add(F, W) — the first wire added to a wire-less face
/// becomes the outer wire (BRep_Builder.cxx L491-505); later wires go to
/// the inner list. The rcad add_tface marker wire (an empty wire) is
/// replaced.
pub(crate) fn builder_add_face_wire(pool: &mut BRep, the_face: &mut Shape, the_wire: Shape) {
    let fd = pool.face_mut(the_face.clone());
    let outer_is_empty = match fd.outer_wire.data.as_ref() {
        TShape::Wire(wd) => wd.edges.is_empty(),
        _ => true,
    };
    if outer_is_empty {
        fd.outer_wire = the_wire.clone();
        fd.my_shapes.push(the_wire);
    } else {
        fd.inner_wires.push(the_wire.clone());
        fd.my_shapes.push(the_wire);
    }
}

/// OCCT LocOpe_SplitDrafts (LocOpe_SplitDrafts.hxx L38-111) — splits a face
/// of a shape with a wire and puts a draft angle on both sides of the wire.
pub struct LocOpeSplitDrafts {
    my_shape: Shape,   // OCCT: myShape (TopoDS_Shape; null = Shape::null())
    my_result: Shape,  // OCCT: myResult
    my_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)>, // OCCT: myMap
}

impl Default for LocOpeSplitDrafts {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeSplitDrafts {
    /// OCCT LocOpe_SplitDrafts::LocOpe_SplitDrafts() (hxx L44).
    pub fn new() -> Self {
        LocOpeSplitDrafts {
            my_shape: Shape::null(),
            my_result: Shape::null(),
            my_map: IndexMap::new(),
        }
    }

    /// OCCT LocOpe_SplitDrafts::LocOpe_SplitDrafts(S) (hxx L47-50).
    pub fn new_with_shape(the_s: &Shape) -> Self {
        LocOpeSplitDrafts {
            my_shape: the_s.clone(),
            my_result: Shape::null(),
            my_map: IndexMap::new(),
        }
    }

    /// OCCT LocOpe_SplitDrafts::Init(S) (cxx L86-91).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_shape = the_s.clone();
        self.my_result = Shape::null(); // OCCT L89: myResult.Nullify()
        self.my_map.clear(); // OCCT L90: myMap.Clear()
    }

    /// OCCT LocOpe_SplitDrafts::Perform(F, W, Extr, NPl, Angle)
    /// (cxx L95-102) — the left-only overload.
    pub fn perform(
        &mut self,
        the_f: &Shape,
        the_w: &Shape,
        the_extr: DVec3,
        the_npl: &Plane,
        the_angle: f64,
    ) {
        self.perform_full(
            the_f, the_w, the_extr, the_npl, the_angle, the_extr, the_npl, the_angle, true, false,
        );
    }

    /// OCCT LocOpe_SplitDrafts::Perform(F, W, Extrg, NPlg, Angleg, Extrd,
    /// NPld, Angled, ModLeft, ModRight) (cxx L106-1434).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_full(
        &mut self,
        the_f: &Shape,
        the_w: &Shape,
        extrg: DVec3,
        nplg: &Plane,
        angleg: f64,
        extrd: DVec3,
        npld: &Plane,
        angled: f64,
        mod_left: bool,
        mod_right: bool,
    ) {
        // OCCT L118: int j;
        let mut j: i32;

        // OCCT L120-121.
        self.my_result = Shape::null();
        self.my_map.clear();

        // OCCT L122-125.
        if self.my_shape.is_null() || the_f.is_null() || the_w.is_null() {
            panic!("Standard_NullObject");
        }

        // OCCT L127-130.
        if !mod_left && !mod_right {
            panic!("Standard_ConstructionError");
        }

        // OCCT L132: TopAbs_Orientation OriF = TopAbs_FORWARD;
        let mut ori_f = Orientation::Forward;

        // OCCT L134: bool FinS = false; TopExp_Explorer exp, exp2;
        let mut fin_s = false;
        let mut exp = TopExpExplorer::default();
        let mut exp2 = TopExpExplorer::default();

        // OCCT L136-146.
        exp.init(&self.my_shape, ShapeType::Face, ShapeType::Shape);
        while exp.more() {
            let fac = exp.current();
            let thelist: Vec<Shape> = Vec::new();
            self.my_map.insert(shape_key(&fac), (fac.clone(), thelist));
            if fac.is_same(the_f) {
                ori_f = fac.orientation;
                fin_s = true;
            }
            exp.next();
        }

        // OCCT L148-152.
        if !fin_s {
            println!("LocOpe_SplitDrafts:!Fins throw Standard_ConstructionError()");
            panic!("Standard_ConstructionError");
        }

        // OCCT L154: gp_Pln NewPlg, NewPld; gp_Ax1 NormalFg, NormalFd;
        let mut new_plg = Plane::new(DVec3::ZERO, DVec3::Z);
        let mut new_pld = Plane::new(DVec3::ZERO, DVec3::Z);
        let mut normal_fg = Ax1::new(DVec3::ZERO, DVec3::Z);
        let mut normal_fd = Ax1::new(DVec3::ZERO, DVec3::Z);

        // OCCT L156: TopoDS_Shape aLocalFace = F.Oriented(OriF);
        let a_local_face = with_orientation(the_f, ori_f);

        // OCCT L158-166.
        if !new_plane(
            &a_local_face, extrg, nplg, angleg, &mut new_plg, &mut normal_fg, mod_left,
        ) || !new_plane(
            &a_local_face, extrd, npld, angled, &mut new_pld, &mut normal_fd, mod_right,
        ) {
            return;
        }

        // OCCT L168: NCollection_List<TopoDS_Shape>::Iterator itl;
        // OCCT L169: BRep_Builder B; (architecture difference #6)
        let mut pool = BRep::new();
        let mut b = BRepBuilder::new();

        // OCCT L171-174.
        let new_sg = Surface3::Plane(new_plg);
        let new_sd = Surface3::Plane(new_pld);
        let the_line_pipe = Curve3::Line(Line3::new(normal_fg.location, normal_fg.direction));
        let mut i2s = GeomIntIntSS::new_with_surfaces(&new_sg, &new_sd, CONFUSION);

        // OCCT L176-179.
        let mut the_map: HashSet<(u64, u32)> = HashSet::new();
        let mut has = BRepAdaptorSurface::new(new_sg.clone()); // HAS->Load
        let mut intcs = IntCurveSurfaceHInter::new();

        // OCCT L181: TopoDS_Wire theW = W;
        let mut the_w = the_w.clone();

        // OCCT L182-331.
        if i2s.is_done() && i2s.nb_lines() > 0 {
            // on split le wire

            // OCCT L186-193.
            let mut the_pipe = GeomFillPipe::new();
            let pipe_profile = i2s
                .line(1)
                .expect("GeomInt_IntSS::Line(1) out of range");
            the_pipe.set_generate_particular_case(true);
            // OCCT L188: thePipe.Init(TheLinePipe, i2s.Line(1)) — the
            // Init(Path, FirstSect, Option = IsCorrectedFrenet) form.
            the_pipe.init_with_trihedron(
                &the_line_pipe.clone(),
                &pipe_profile,
                Trihedron::IsCorrectedFrenet,
            );
            the_pipe.perform(true, false);
            if !the_pipe.is_done() {
                panic!("Standard_ConstructionError: GeomFill_Pipe : Cannot make a surface");
            }

            // OCCT L195-196: Spl = thePipe.Surface(); HAS->Load(Spl).
            let spl = the_pipe.surface().expect("GeomFill_Pipe::Surface").clone();
            has = BRepAdaptorSurface::new(spl);

            // OCCT L198: LocOpe_SplitShape splw(W);
            let mut splw = LocOpeSplitShape::new_with_shape(&the_w);

            // OCCT L200-296.
            exp.init(&the_w, ShapeType::Edge, ShapeType::Shape);
            while exp.more() {
                let edg = exp.current(); // TopoDS::Edge(exp.Current())
                if the_map.insert(shape_key(&edg)) {
                    // OCCT L205-207.
                    let (c, f, l) = brep_tool_curve(&edg).expect("BRep_Tool::Curve");
                    let hac = BRepAdaptorCurve::new(c);
                    intcs.perform(&hac, &has);
                    if !intcs.is_done() {
                        exp.next();
                        continue; // voir ce qu`on peut faire de mieux
                    }

                    // OCCT L215-218.
                    if intcs.nb_segments() >= 2 {
                        exp.next();
                        continue; // Not yet implemented...and probably never
                    }

                    // OCCT L220-260.
                    if intcs.nb_segments() == 1 {
                        // OCCT L222-225: the first/second points of the
                        // segment (the base Intersection stores the 3D
                        // points; the bop wrapper tuple carries (w,u,v)).
                        let p1 = intcs.base.segment(1).first_point();
                        let p2 = intcs.base.segment(1).second_point();
                        let pf = p1.pnt();
                        let pl = p2.pnt();
                        // OCCT L226-227.
                        let (vf, vl) = topexp_vertices(&edg);
                        let pf_v = brep_tool_pnt(&vf);
                        let pl_v = brep_tool_pnt(&vl);
                        let mut tolf = brep_tool_tolerance(&vf);
                        let mut toll = brep_tool_tolerance(&vl);
                        tolf *= tolf;
                        toll *= toll;

                        // OCCT L235-238.
                        let dff = pf.distance_squared(pf_v);
                        let dfl = pf.distance_squared(pl_v);
                        let dlf = pl.distance_squared(pf_v);
                        let dll = pl.distance_squared(pl_v);

                        // OCCT L240-259.
                        if (dff <= tolf && dll <= toll) || (dlf <= tolf && dfl <= toll) {
                            exp.next();
                            continue;
                        } else {
                            // on segmente edg en pf et pl
                            let vnewf = make_vertex_confusion(&mut pool, pf);
                            let vnewl = make_vertex_confusion(&mut pool, pl);
                            if p1.w() >= f && p1.w() <= l && p2.w() >= f && p2.w() <= l {
                                splw.add(&vnewf, p1.w(), &edg);
                                splw.add(&vnewl, p2.w(), &edg);
                            } else {
                                exp.next();
                                continue;
                            }
                        }
                    }
                    // OCCT L261-294.
                    else if intcs.nb_points() != 0 {
                        // OCCT L263-270.
                        let (vf, vl) = topexp_vertices(&edg);
                        let pf = brep_tool_pnt(&vf);
                        let pl = brep_tool_pnt(&vl);
                        let mut tolf = brep_tool_tolerance(&vf);
                        let mut toll = brep_tool_tolerance(&vl);
                        tolf *= tolf;
                        toll *= toll;

                        // OCCT L272-293.
                        for i in 1..=intcs.nb_points() {
                            let pi = intcs.base.point(i);
                            let pi_pt = pi.pnt();
                            let dif = pi_pt.distance_squared(pf);
                            let dil = pi_pt.distance_squared(pl);
                            if dif <= tolf {
                            } else if dil <= toll {
                            } else {
                                if pi.w() >= f && pi.w() <= l {
                                    let vnew = make_vertex_confusion(&mut pool, pi_pt);
                                    splw.add(&vnew, pi.w(), &edg);
                                }
                            }
                        }
                    }
                    let _ = (f, l);
                }
                exp.next();
            }

            // OCCT L298-308.
            let lres = splw
                .descendant_shapes(&the_w)
                .cloned()
                .unwrap_or_default();
            if lres.len() != 1 {
                return;
            }

            if !the_w.is_same(&lres[0]) {
                the_w = Shape::null(); // OCCT L306: theW.Nullify()
                the_w = lres[0].clone(); // OCCT L307: TopoDS::Wire(lres.First())
            }

            // OCCT L310-330.
            exp.re_init();
            while exp.more() {
                let cur = exp.current();
                if !self.my_map.contains_key(&shape_key(&cur)) {
                    self.my_map.insert(shape_key(&cur), (cur.clone(), Vec::new()));
                    // OCCT L316-319.
                    for itl in splw
                        .descendant_shapes(&cur)
                        .cloned()
                        .unwrap_or_default()
                    {
                        if let Some(entry) = self.my_map.get_mut(&shape_key(&cur)) {
                            entry.1.push(itl);
                        }
                    }
                    // OCCT L320-328.
                    exp2.init(&cur, ShapeType::Vertex, ShapeType::Shape);
                    while exp2.more() {
                        let vtx = exp2.current();
                        if !self.my_map.contains_key(&shape_key(&vtx)) {
                            self.my_map
                                .insert(shape_key(&vtx), (vtx.clone(), vec![vtx.clone()]));
                        }
                        exp2.next();
                    }
                }
                exp.next();
            }
        }
        // OCCT L332-352.
        else {
            exp.init(&the_w, ShapeType::Edge, ShapeType::Shape);
            while exp.more() {
                let cur = exp.current();
                if !self.my_map.contains_key(&shape_key(&cur)) {
                    self.my_map
                        .insert(shape_key(&cur), (cur.clone(), vec![cur.clone()]));
                    // OCCT L341-349.
                    exp2.init(&cur, ShapeType::Vertex, ShapeType::Shape);
                    while exp2.more() {
                        let vtx = exp2.current();
                        if !self.my_map.contains_key(&shape_key(&vtx)) {
                            self.my_map
                                .insert(shape_key(&vtx), (vtx.clone(), vec![vtx.clone()]));
                        }
                        exp2.next();
                    }
                }
                exp.next();
            }
        }

        // On split la face par le wire

        // OCCT L356-358.
        let mut won_s = LocOpeWiresOnShape::new(&self.my_shape);
        let mut spls = LocOpeSpliter::new_with_shape(&self.my_shape);
        won_s.bind_wire_face(&the_w, the_f);

        // OCCT L366-370.
        spls.perform(&mut won_s);
        if !spls.is_done() {
            return;
        }

        // OCCT L372-373.
        let res = spls
            .resulting_shape()
            .cloned()
            .expect("StdFail_NotDone: ResultingShape");
        let the_left: Vec<Shape> = spls.direct_left().clone();

        // Descendants
        // OCCT L376-383.
        exp.init(&self.my_shape, ShapeType::Face, ShapeType::Shape);
        while exp.more() {
            let fac = exp.current();
            let desc = spls
                .descendant_shapes(&fac)
                .cloned()
                .unwrap_or_default();
            for itl in desc {
                if let Some(entry) = self.my_map.get_mut(&shape_key(&fac)) {
                    entry.1.push(itl);
                }
            }
            exp.next();
        }

        // OCCT L385-399.
        let mut map_w: HashMap<(u64, u32), Shape> = HashMap::new();
        exp.init(&the_w, ShapeType::Edge, ShapeType::Shape);
        while exp.more() {
            let cur = exp.current();
            if !map_w.contains_key(&shape_key(&cur)) {
                map_w.insert(shape_key(&cur), Shape::null());
                // OCCT L391-397.
                exp2.init(&cur, ShapeType::Vertex, ShapeType::Shape);
                while exp2.more() {
                    let vtx = exp2.current();
                    if !map_w.contains_key(&shape_key(&vtx)) {
                        map_w.insert(shape_key(&vtx), Shape::null());
                    }
                    exp2.next();
                }
            }
            exp.next();
        }

        // OCCT L401-403.
        let mut the_map_ef: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        map_shapes_and_ancestors(&res, ShapeType::Edge, ShapeType::Face, &mut the_map_ef);

        // On stocke les geometries potentiellement generees par les edges
        // OCCT L406-408: MapEV (genere); MapSg, MapSd (image a gauche et a
        // droite).
        let mut map_ev: IndexMap<(u64, u32), (Shape, Shape)> = IndexMap::new();
        let mut map_sg: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        let mut map_sd: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();

        // OCCT L410: int Nbedges, index;
        let nbedges: usize;
        let mut index: usize;

        // OCCT L411-469.
        let fkey = shape_key(the_f);
        let my_map_f = self
            .my_map
            .get(&fkey)
            .expect("NoSuchObject: myMap(F)")
            .1
            .clone();
        for itl in &my_map_f {
            let fac = itl.clone(); // TopoDS::Face(itl.Value())
            exp.init(&fac, ShapeType::Edge, ShapeType::Shape);
            while exp.more() {
                let edg = exp.current();
                // OCCT L417-420: MapEV.FindIndex(edg) != 0 -> continue.
                if map_ev.get_index_of(&shape_key(&edg)).is_some() {
                    exp.next();
                    continue;
                }
                // OCCT L421-450.
                if map_w.contains_key(&shape_key(&edg)) {
                    // edge du wire initial
                    // OCCT L423-430.
                    let (c0, f, l) = brep_tool_curve(&edg).expect("BRep_Tool::Curve");
                    let mut c = geom_trimmed_basis_curve(&c0);
                    c = geom_transformed_identity(&c);
                    let _ = (f, l);
                    // OCCT L432-439.
                    let mut the_pipe = GeomFillPipe::new();
                    the_pipe.set_generate_particular_case(true);
                    // OCCT L433: thePipe.Init(TheLinePipe, C) — the
                    // Init(Path, FirstSect, Option = IsCorrectedFrenet) form.
                    the_pipe.init_with_trihedron(
                        &the_line_pipe.clone(),
                        &c,
                        Trihedron::IsCorrectedFrenet,
                    );
                    the_pipe.perform(true, false);
                    if !the_pipe.is_done() {
                        panic!("Standard_ConstructionError: GeomFill_Pipe : Cannot make a surface");
                    }

                    // OCCT L441-445.
                    let the_ps = geom_rectangular_trimmed_basis_surface(
                        the_pipe.surface().expect("GeomFill_Pipe::Surface"),
                    );

                    // OCCT L447-449.
                    let new_face = make_face_from_surface(&mut pool, Some(the_ps), CONFUSION);
                    map_ev.insert(shape_key(&edg), (edg.clone(), new_face));
                }
                // OCCT L451-467.
                else {
                    // on recupere la face.
                    let index0 = the_map_ef
                        .get_index_of(&shape_key(&edg))
                        .expect("NoSuchObject: theMapEF.FindIndex");
                    if the_map_ef.get_index(index0).expect("theMapEF entry").1.1.len() != 2 {
                        return; // NotDone
                    }
                    let entry = the_map_ef.get_index(index0).expect("theMapEF entry").1.clone();
                    let the_face = Shape::null();
                    let _ = the_face;
                    // OCCT L459: theMapEF(index).First() / .Last() — the
                    // FIRST/LAST of the ancestor face list.
                    if entry.1[0].is_same(&fac) {
                        map_ev.insert(shape_key(&edg), (edg.clone(), entry.1[1].clone()));
                    } else {
                        map_ev.insert(shape_key(&edg), (edg.clone(), entry.1[0].clone()));
                    }
                }
                exp.next();
            }
        }

        // OCCT L471: NCollection_DataMap<...> MapSonS;
        let mut map_son_s: IndexMap<(u64, u32), (Shape, Shape)> = IndexMap::new();

        // OCCT L473-474.
        nbedges = map_ev.len();
        index = 1;
        while index <= nbedges {
            // OCCT L476-482.
            let index_key = map_ev.get_index(index - 1).expect("MapEV entry").1.0.clone();
            exp.init(&index_key, ShapeType::Vertex, ShapeType::Shape);
            while exp.more() {
                let vtx = exp.current(); // TopoDS::Vertex(exp.Current())
                // OCCT L479-482.
                if map_ev.get_index_of(&shape_key(&vtx)).is_some() {
                    exp.next();
                    continue;
                }

                // Localisation du vertex :
                //    - entre 2 edges d`origine : on recupere l`edge qui n`est
                //                                pas dans F
                //    - entre 2 edges du wire   : droite
                //    - mixte                   : intersection de surfaces
                // OCCT L489-507.
                j = 1;
                while j <= nbedges as i32 {
                    if j == index as i32 {
                        j += 1;
                        continue;
                    }
                    let jkey = map_ev
                        .get_index((j - 1) as usize)
                        .expect("MapEV entry")
                        .1.0
                        .clone();
                    let mut found = false;
                    exp2.init(&jkey, ShapeType::Vertex, ShapeType::Shape);
                    while exp2.more() {
                        let vtx2 = exp2.current();
                        if vtx2.is_same(&vtx) {
                            found = true;
                            break;
                        }
                        exp2.next();
                    }
                    if found {
                        break;
                    }
                    j += 1;
                }

                // OCCT L508-518.
                let mut choice: i32 = 0;
                let edg1 = index_key.clone();
                let mut edg2 = Shape::null();
                if j <= nbedges as i32 {
                    edg2 = map_ev
                        .get_index((j - 1) as usize)
                        .expect("MapEV entry")
                        .1.0
                        .clone();
                } else {
                    edg2 = edg1.clone();
                }

                // OCCT L519-548.
                if map_w.contains_key(&shape_key(&edg1)) {
                    if j > nbedges as i32 {
                        choice = 2; // doit correspondre a edge ferme — droite
                    } else if map_w.contains_key(&shape_key(&edg2)) {
                        choice = 2; // droite
                    } else {
                        choice = 3; // mixte
                    }
                } else if j > nbedges as i32 {
                    choice = 1; // doit correspondre a edge ferme — edge a retrouver
                } else if !map_w.contains_key(&shape_key(&edg2)) {
                    choice = 1; // edge a retrouver
                } else {
                    choice = 3; // mixte
                }

                // OCCT L549-552.
                let mut newc: Option<Curve3> = None;
                let mut new_cs1: Option<Curve2d> = None;
                let mut new_cs2: Option<Curve2d> = None;
                let mut knownp = 0.0f64;
                let mut e_bind = Shape::null();

                // OCCT L553-679.
                match choice {
                    1 => {
                        // OCCT L555-594.
                        exp2.init(&res, ShapeType::Edge, ShapeType::Shape);
                        while exp2.more() {
                            let cur = exp2.current();
                            if cur.is_same(&edg1) || cur.is_same(&edg2) {
                                exp2.next();
                                continue;
                            }
                            // OCCT L562-574: exp3 explores the
                            // FORWARD-oriented edge vertices.
                            let fwd = with_orientation(&cur, Orientation::Forward);
                            let mut exp3 = TopExpExplorer::default();
                            exp3.init(&fwd, ShapeType::Vertex, ShapeType::Shape);
                            while exp3.more() {
                                if exp3.current().is_same(&vtx) {
                                    break;
                                }
                                exp3.next();
                            }
                            if exp3.more() {
                                break;
                            }
                            exp2.next();
                        }
                        if exp2.more() {
                            // OCCT L577-584.
                            let cur = exp2.current();
                            let (c0, _f, _l) = brep_tool_curve(&cur).expect("BRep_Tool::Curve");
                            let c = geom_transformed_identity(&c0);
                            newc = Some(c);
                            e_bind = cur;
                            knownp = brep_tool_parameter(&pool, &vtx, &e_bind);
                        } else {
                            // droite ??? il vaudrait mieux sortir
                            return;
                        }
                    }
                    2 => {
                        // OCCT L596-601.
                        let mut the_line = Line3::new(normal_fg.location, normal_fg.direction);
                        // OCCT: theLine.Translate(NormalFg.Location(),
                        // BRep_Tool::Pnt(vtx)) — the position moves by
                        // P2 - P1 (gp_Ax1::Translate).
                        the_line.origin = brep_tool_pnt(&vtx);
                        newc = Some(Curve3::Line(the_line));
                        knownp = 0.0;
                    }
                    3 => {
                        // OCCT L603-678.
                        let f1 = map_ev
                            .get(&shape_key(&edg1))
                            .expect("MapEV.FindFromKey")
                            .1
                            .clone();
                        let f2 = map_ev
                            .get(&shape_key(&edg2))
                            .expect("MapEV.FindFromKey")
                            .1
                            .clone();
                        let s1 = brep_tool_surface(&f1).expect("BRep_Tool::Surface");
                        let s2 = brep_tool_surface(&f2).expect("BRep_Tool::Surface");
                        let mut app_s1 = false;
                        let mut app_s2 = false;
                        if !matches!(s1, Surface3::Plane(_)) {
                            app_s1 = true;
                        }
                        if !matches!(s2, Surface3::Plane(_)) {
                            app_s2 = true;
                        }
                        i2s.perform(&s1, &s2, CONFUSION, true, app_s1, app_s2);
                        if !i2s.is_done() || i2s.nb_lines() == 0 {
                            return;
                        }

                        // OCCT L624-627.
                        let mut pmin = 0.0f64;
                        let mut dist2;
                        let mut dist2_min;
                        let mut glob2_min = f64::MAX; // RealLast()

                        // OCCT L628-629.
                        let mut imin: i32;
                        let pv = brep_tool_pnt(&vtx);
                        imin = 0;
                        for i in 1..=i2s.nb_lines() {
                            let crv = i2s.line(i).expect("GeomInt_IntSS::Line");
                            // OCCT L632-633: TheCurve.Load(i2s.Line(i));
                            // Extrema_ExtPC myExtPC(pv, TheCurve).
                            let dom = crv.default_domain();
                            let my_ext_pc = ExtPC::new(pv, &crv, CONFUSION, dom[0], dom[1]);

                            if my_ext_pc.is_done() {
                                // OCCT L637-639.
                                let mut thepmin = dom[0]; // TheCurve.FirstParameter()
                                let (d2min0, d2last, _p1b, _p2b) =
                                    ext_trimmed_square_distances(pv, &crv);
                                dist2_min = d2min0;
                                dist2 = d2last;
                                if dist2 < dist2_min {
                                    thepmin = dom[1]; // TheCurve.LastParameter()
                                }
                                // OCCT L644-652.
                                for k in 1..=my_ext_pc.nb_ext() {
                                    dist2 = my_ext_pc.square_distance(k);
                                    if dist2 < dist2_min {
                                        dist2_min = dist2;
                                        thepmin = my_ext_pc.point(k).param;
                                    }
                                }

                                // OCCT L654-659.
                                if dist2_min < glob2_min {
                                    glob2_min = dist2_min;
                                    pmin = thepmin;
                                    imin = i as i32;
                                }
                            }
                        }
                        // OCCT L662-665.
                        if imin == 0 {
                            return;
                        }

                        // OCCT L667-676.
                        newc = i2s.line(imin as usize);
                        knownp = pmin;
                        if app_s1 {
                            new_cs1 = i2s.line_on_s1(imin as usize);
                        }
                        if app_s2 {
                            new_cs2 = i2s.line_on_s2(imin as usize);
                        }
                    }
                    _ => {}
                }

                // Determination des vertex par intersection sur Plg ou/et Pld
                // OCCT L683-724.
                let hac = BRepAdaptorCurve::new(newc.clone().expect("Newc"));
                let mut nbfois = 2i32;
                let mut vtx1 = Shape::null();
                let mut vtx2 = Shape::null();
                let mut p1 = 0.0f64;
                let mut p2 = 0.0f64;
                let mut is_left = false;
                if choice == 1 {
                    // edge retrouve : on ne fait qu`une seule intersection
                    // il faut utiliser Plg ou Pld

                    // OCCT L693-719.
                    let indedgf = the_map_ef
                        .get_index_of(&shape_key(&edg1))
                        .expect("NoSuchObject: theMapEF.FindIndex");
                    let ancestor_list = the_map_ef
                        .get_index(indedgf)
                        .expect("theMapEF entry")
                        .1
                        .1
                        .clone();
                    let mut hit = false;
                    for itl in &ancestor_list {
                        if contains(self.my_map_f_list(&fkey), itl) {
                            if contains(&the_left, itl) {
                                has = BRepAdaptorSurface::new(new_sg.clone());
                                is_left = true;
                            } else {
                                has = BRepAdaptorSurface::new(new_sd.clone());
                                is_left = false;
                            }

                            nbfois = 1;
                            vtx2 = vtx.clone();
                            p2 = knownp;
                            hit = true;
                            break;
                        }
                    }
                    if !hit {
                        println!("LocOpe_SplitDrafts: betite probleme ");
                        return;
                    }
                } else {
                    has = BRepAdaptorSurface::new(new_sg.clone());
                }

                // OCCT L726-760.
                for it in 1..=nbfois {
                    if it == 2 {
                        has = BRepAdaptorSurface::new(new_sd.clone());
                    }

                    intcs.perform(&hac, &has);
                    if !intcs.is_done() {
                        // voir ce qu`on peut faire de mieux
                        return;
                    }
                    let mut imin = 1usize;
                    let mut delta = (knownp - intcs.base.point(1).w()).abs();
                    for i in 2..=intcs.nb_points() {
                        let newdelta = (knownp - intcs.base.point(i).w()).abs();
                        if newdelta < delta {
                            imin = i;
                            delta = newdelta;
                        }
                    }
                    if it == 1 {
                        vtx1 = make_vertex_confusion(&mut pool, intcs.base.point(imin).pnt());
                        p1 = intcs.base.point(imin).w();
                        knownp = p1;
                    } else {
                        vtx2 = make_vertex_confusion(&mut pool, intcs.base.point(imin).pnt());
                        p2 = intcs.base.point(imin).w();
                    }
                }

                // OCCT L761-837.
                if (p1 - p2).abs() > PCONFUSION {
                    // OCCT L763-776: MakeEdge(Newc) + Add(vtx1/vtx2 with the
                    // OCCT orientations) + UpdateVertex(p1/p2).
                    let (first, last) = if p1 < p2 {
                        (vtx1.clone(), vtx2.clone())
                    } else {
                        (vtx2.clone(), vtx1.clone())
                    };
                    let new_edge_s = b.add_edge(
                        &mut pool,
                        newc.clone(),
                        first,
                        last,
                        [p1.min(p2), p1.max(p2)],
                    );
                    b.set_vertex_param(&mut pool, new_edge_s.clone(), vtx1.clone(), p1);
                    b.set_vertex_param(&mut pool, new_edge_s.clone(), vtx2.clone(), p2);

                    // OCCT L777-791.
                    if let Some(cs1) = &new_cs1 {
                        b.update_edge_pcurve(
                            &mut pool,
                            new_edge_s.clone(),
                            cs1.clone(),
                            map_ev
                                .get(&shape_key(&edg1))
                                .expect("MapEV.FindFromKey")
                                .1
                                .clone(),
                            CONFUSION,
                        );
                    }

                    if let Some(cs2) = &new_cs2 {
                        b.update_edge_pcurve(
                            &mut pool,
                            new_edge_s.clone(),
                            cs2.clone(),
                            map_ev
                                .get(&shape_key(&edg2))
                                .expect("MapEV.FindFromKey")
                                .1
                                .clone(),
                            CONFUSION,
                        );
                    }

                    // OCCT L793: MapEV.Add(vtx, NewEdge);
                    map_ev.insert(shape_key(&vtx), (vtx.clone(), new_edge_s));

                    // OCCT L795-828.
                    if choice == 1 {
                        // OCCT L797-798.
                        let ne = pool.empty_copied(&e_bind);
                        // OCCT L800-814.
                        exp2.init(&e_bind, ShapeType::Vertex, ShapeType::Shape);
                        while exp2.more() {
                            let thevtx = exp2.current(); // TopoDS::Vertex
                            if thevtx.is_same(&vtx) {
                                b.add_to_edge(
                                    &mut pool,
                                    ne.clone(),
                                    with_orientation(&vtx1, thevtx.orientation),
                                );
                                b.update_vertex_on_edge(
                                    &mut pool,
                                    vtx1.clone(),
                                    p1,
                                    ne.clone(),
                                    CONFUSION,
                                );
                            } else {
                                b.add_to_edge(&mut pool, ne.clone(), thevtx.clone());
                                let theprm = brep_tool_parameter(&pool, &thevtx, &e_bind);
                                b.update_vertex_on_edge(
                                    &mut pool,
                                    thevtx.clone(),
                                    theprm,
                                    ne.clone(),
                                    brep_tool_tolerance(&thevtx),
                                );
                            }
                            exp2.next();
                        }
                        map_son_s.insert(
                            shape_key(&e_bind),
                            (e_bind.clone(), with_orientation(&ne, Orientation::Forward)),
                        );
                        if is_left {
                            map_sg.insert(shape_key(&vtx), vec![vtx1.clone()]);
                        } else {
                            map_sd.insert(shape_key(&vtx), vec![vtx1.clone()]);
                        }
                    } else {
                        // OCCT L829-836.
                        map_sg.insert(shape_key(&vtx), vec![vtx1.clone()]);
                        map_sd.insert(shape_key(&vtx), vec![vtx2.clone()]);
                    }
                }
                // OCCT L838-864.
                else {
                    map_ev.insert(shape_key(&vtx), (vtx.clone(), vtx2.clone()));
                    if choice == 1 {
                        if is_left {
                            map_sg.insert(shape_key(&vtx), vec![vtx.clone()]);
                        } else {
                            map_sd.insert(shape_key(&vtx), vec![vtx.clone()]);
                        }
                    } else {
                        map_sg.insert(shape_key(&vtx), vec![vtx2.clone()]);
                        map_sd.insert(shape_key(&vtx), vec![vtx2.clone()]);
                    }
                }
                exp.next();
            }
            index += 1;
        }

        // OCCT L868-1100.
        the_map.clear();
        exp.init(&the_w, ShapeType::Edge, ShapeType::Shape);
        while exp.more() {
            let edg = exp.current(); // TopoDS::Edge(exp.Current())
            if !the_map.insert(shape_key(&edg)) {
                // precaution sans doute inutile...
                exp.next();
                continue;
            }
            // OCCT L876-877.
            let indedg = map_ev
                .get_index_of(&shape_key(&edg))
                .expect("NoSuchObject: MapEV.FindIndex");
            let mut gen_f = map_ev.get_index(indedg).expect("MapEV entry").1.1.clone();

            // OCCT L878-880.
            map_sg.insert(shape_key(&edg), Vec::new()); // genere a gauche
            map_sd.insert(shape_key(&edg), Vec::new()); // genere a droite

            // OCCT L881-886.
            let edg_fwd = with_orientation(&edg, Orientation::Forward);
            let (vf, vl) = topexp_vertices(&edg_fwd);
            let gvf = map_ev
                .get(&shape_key(&vf))
                .expect("MapEV.FindFromKey")
                .1
                .clone();
            let gvl = map_ev
                .get(&shape_key(&vl))
                .expect("MapEV.FindFromKey")
                .1
                .clone();

            // OCCT L997-1019: nouveau code
            let mut vfg = Shape::null();
            let mut vfd = Shape::null();
            let mut vlg = Shape::null();
            let mut vld = Shape::null();
            if gvf.shape_type() == ShapeType::Vertex {
                vfg = gvf.clone();
                vfd = vfg.clone();
            } else {
                vfg = map_sg
                    .get(&shape_key(&vf))
                    .expect("MapSg(Vf)")
                    .first()
                    .cloned()
                    .expect("MapSg(Vf).First()");
                vfd = map_sd
                    .get(&shape_key(&vf))
                    .expect("MapSd(Vf)")
                    .first()
                    .cloned()
                    .expect("MapSd(Vf).First()");
            }
            if gvl.shape_type() == ShapeType::Vertex {
                vlg = gvl.clone();
                vld = vlg.clone();
            } else {
                vlg = map_sg
                    .get(&shape_key(&vl))
                    .expect("MapSg(Vl)")
                    .first()
                    .cloned()
                    .expect("MapSg(Vl).First()");
                vld = map_sd
                    .get(&shape_key(&vl))
                    .expect("MapSd(Vl)")
                    .first()
                    .cloned()
                    .expect("MapSd(Vl).First()");
            }

            // OCCT L1021-1031.
            let new_edgg = new_edge(
                &edg, &gen_f, &new_sg, &vfg, &vlg, &mut pool, &mut b,
            );
            if new_edgg.is_null() {
                return;
            }

            let new_edgd = new_edge(
                &edg, &gen_f, &new_sd, &vfd, &vld, &mut pool, &mut b,
            );
            if new_edgg.is_null() {
                // OCCT L1028: the source tests NewEdgg a second time — kept
                // as-is.
                return;
            }

            // OCCT L1033-1078.
            let mut isedg = false;
            if gvf.shape_type() == ShapeType::Vertex && gvl.shape_type() == ShapeType::Vertex {
                // edg ou face a 2 cotes

                // Comparaison NewEdgg et NewEdgd
                // OCCT L1039-1045.
                let (cg, fg, lg) = brep_tool_curve(&new_edgg).expect("BRep_Tool::Curve");
                let (cd, fd, ld) = brep_tool_curve(&new_edgd).expect("BRep_Tool::Curve");
                let prmg = (fg + lg) / 2.0;
                let prmd = (fd + ld) / 2.0;
                let pg = cg.point_at(prmg);
                let pd = cd.point_at(prmd);
                // OCCT L1046: the source tests NewEdgg twice — kept as-is.
                let tol = brep_tool_tolerance(&new_edgg).max(brep_tool_tolerance(&new_edgg));
                if pg.distance_squared(pd) <= tol * tol {
                    isedg = true;
                    // raffinement pour essayer de partager l`edge de depart...
                    // OCCT L1051-1063.
                    let mut modified = true;
                    if gvf.is_same(&vf) && gvl.is_same(&vl) {
                        // Comparaison avec l`edge de depart
                        let (cd0, fd0, ld0) = brep_tool_curve(&edg).expect("BRep_Tool::Curve");
                        let cd_start = cd0;
                        let prmd_start = (fd0 + ld0) / 2.0;
                        let pd_start = cd_start.point_at(prmd_start);
                        let tol_start = brep_tool_tolerance(&new_edgg)
                            .max(brep_tool_tolerance(&edg));
                        if pg.distance_squared(pd_start) <= tol_start * tol_start {
                            modified = false;
                        }
                    }

                    // OCCT L1065-1076.
                    if !modified {
                        map_w.insert(shape_key(&edg), edg.clone());
                        map_sg.get_mut(&shape_key(&edg)).expect("MapSg(edg)").push(edg.clone());
                        map_sd.get_mut(&shape_key(&edg)).expect("MapSd(edg)").push(edg.clone());
                    } else {
                        map_w.insert(shape_key(&edg), new_edgg.clone());
                        map_sg
                            .get_mut(&shape_key(&edg))
                            .expect("MapSg(edg)")
                            .push(new_edgg.clone());
                        map_sd
                            .get_mut(&shape_key(&edg))
                            .expect("MapSd(edg)")
                            .push(new_edgg.clone());
                    }
                }
                let _ = (cg, cd);
            }

            // OCCT L1080-1099.
            if !isedg {
                // face a 2 ou 3 ou 4 cotes
                map_sg
                    .get_mut(&shape_key(&edg))
                    .expect("MapSg(edg)")
                    .push(new_edgg.clone());
                map_sd
                    .get_mut(&shape_key(&edg))
                    .expect("MapSd(edg)")
                    .push(new_edgd.clone());

                let mut theedges: Vec<Shape> = Vec::new();
                theedges.push(new_edgg.clone());
                theedges.push(new_edgd.clone());
                if gvf.shape_type() == ShapeType::Edge {
                    theedges.push(gvf.clone());
                }
                if gvl.shape_type() == ShapeType::Edge {
                    theedges.push(gvl.clone());
                }
                make_face(&mut gen_f, &mut theedges, &mut pool, &mut b);
                map_w.insert(shape_key(&edg), gen_f.clone());
            }

            // OCCT L877 reference write-back: GenF is MapEV(indedg).
            map_ev.get_index_mut(indedg).expect("MapEV entry").1.1 = gen_f;
            exp.next();
        }

        // OCCT L1102-1336.
        let mut mapedgadded: HashSet<(u64, u32)> = HashSet::new();
        let mut thefaces: Vec<Shape> = Vec::new();

        let my_map_f2 = self
            .my_map
            .get(&fkey)
            .expect("NoSuchObject: myMap(F)")
            .1
            .clone();
        for itl in &my_map_f2 {
            let mut fac = itl.clone(); // TopoDS::Face(itl.Value())
            the_map.clear();
            // OCCT L1109-1120: elle est FORWARD
            let mut drft_face = Shape::null();
            let is_left;
            if contains(&the_left, &fac) {
                drft_face = make_face_from_surface(
                    &mut pool,
                    Some(new_sg.clone()),
                    brep_tool_tolerance(&fac),
                );
                is_left = true;
            } else {
                drft_face = make_face_from_surface(
                    &mut pool,
                    Some(new_sd.clone()),
                    brep_tool_tolerance(&fac),
                );
                is_left = false;
            }

            // OCCT L1122-1334.
            let mut exp3 = TopExpExplorer::default();
            let fac_fwd = with_orientation(&fac, Orientation::Forward);
            exp3.init(&fac_fwd, ShapeType::Wire, ShapeType::Shape);
            while exp3.more() {
                let wir = exp3.current();
                let new_wire_on_f = b.make_wire(&mut pool);
                let wir_fwd = with_orientation(&wir, Orientation::Forward);
                exp.init(&wir_fwd, ShapeType::Edge, ShapeType::Shape);
                while exp.more() {
                    let edg = exp.current(); // TopoDS::Edge(exp.Current())
                    if !the_map.insert(shape_key(&edg)) {
                        // precaution sans doute inutile...
                        exp.next();
                        continue;
                    }
                    // OCCT L1135-1157.
                    if map_w.contains_key(&shape_key(&edg)) {
                        // edge du wire d`origine
                        let ored0 = edg.orientation;
                        let mut ored = ored0;
                        let imgs = if is_left {
                            map_sg.get(&shape_key(&edg)).expect("MapSg(edg)").clone()
                        } else {
                            map_sd.get(&shape_key(&edg)).expect("MapSd(edg)").clone()
                        };
                        for itld in &imgs {
                            if itld.orientation == Orientation::Reversed {
                                ored = top_abs_reverse(ored);
                            }
                            let a_local_edge = with_orientation(itld, ored);
                            b.add_to_wire(&mut pool, new_wire_on_f.clone(), a_local_edge);
                        }
                    }
                    // OCCT L1158-1330.
                    else {
                        let new_s = if is_left {
                            new_sg.clone()
                        } else {
                            new_sd.clone()
                        };
                        // OCCT L1169-1176.
                        let indedg = map_ev
                            .get_index_of(&shape_key(&edg))
                            .expect("NoSuchObject: MapEV.FindIndex");
                        let gen_f = map_ev.get_index(indedg).expect("MapEV entry").1.1.clone();
                        let edg_fwd = with_orientation(&edg, Orientation::Forward);
                        let (vf, vl) = topexp_vertices(&edg_fwd);
                        let gvf = map_ev
                            .get(&shape_key(&vf))
                            .expect("MapEV.FindFromKey")
                            .1
                            .clone();
                        let gvl = map_ev
                            .get(&shape_key(&vl))
                            .expect("MapEV.FindFromKey")
                            .1
                            .clone();

                        // OCCT L1177-1204.
                        if gvf.shape_type() == ShapeType::Vertex
                            && gvl.shape_type() == ShapeType::Vertex
                        {
                            if !gvf.is_same(&vf) || !gvl.is_same(&vl) {
                                let mut new_edg = new_edge(
                                    &edg, &gen_f, &new_s, &gvf, &gvl, &mut pool, &mut b,
                                );
                                if new_edg.is_null() {
                                    return;
                                }

                                map_son_s.insert(
                                    shape_key(&edg),
                                    (edg.clone(), new_edg.clone()),
                                );

                                if new_edg.orientation == Orientation::Reversed {
                                    new_edg.orientation = top_abs_reverse(edg.orientation);
                                } else {
                                    new_edg.orientation = edg.orientation;
                                }
                                b.add_to_wire(&mut pool, new_wire_on_f.clone(), new_edg);
                            } else {
                                // Frozen???
                                b.add_to_wire(&mut pool, new_wire_on_f.clone(), edg.clone());
                            }
                        }
                        // OCCT L1205-1330.
                        else {
                            let mut vff = Shape::null();
                            let mut vll = Shape::null();
                            if gvf.shape_type() == ShapeType::Vertex {
                                vff = gvf.clone();
                            } else if is_left {
                                vff = map_sg
                                    .get(&shape_key(&vf))
                                    .expect("MapSg(Vf)")
                                    .first()
                                    .cloned()
                                    .expect("MapSg(Vf).First()");
                            } else {
                                vff = map_sd
                                    .get(&shape_key(&vf))
                                    .expect("MapSd(Vf)")
                                    .first()
                                    .cloned()
                                    .expect("MapSd(Vf).First()");
                            }
                            if gvl.shape_type() == ShapeType::Vertex {
                                vll = gvl.clone();
                            } else if is_left {
                                vll = map_sg
                                    .get(&shape_key(&vl))
                                    .expect("MapSg(Vl)")
                                    .first()
                                    .cloned()
                                    .expect("MapSg(Vl).First()");
                            } else {
                                vll = map_sd
                                    .get(&shape_key(&vl))
                                    .expect("MapSd(Vl)")
                                    .first()
                                    .cloned()
                                    .expect("MapSd(Vl).First()");
                            }

                            let mut new_edg = new_edge(
                                &edg, &gen_f, &new_s, &vff, &vll, &mut pool, &mut b,
                            );
                            if new_edg.is_null() {
                                return;
                            }

                            // OCCT L1245-1320.
                            if !map_w.contains_key(&shape_key(&vf))
                                && !map_w.contains_key(&shape_key(&vl))
                            {
                                map_son_s.insert(
                                    shape_key(&edg),
                                    (edg.clone(), new_edg.clone()),
                                );
                            } else if map_w.contains_key(&shape_key(&vf)) {
                                if gvf.shape_type() != ShapeType::Edge
                                    || mapedgadded.contains(&shape_key(&gvf))
                                {
                                    map_son_s.insert(
                                        shape_key(&edg),
                                        (edg.clone(), new_edg.clone()),
                                    );
                                } else {
                                    // OCCT L1262-1284.
                                    let new_wir = b.make_wire(&mut pool);
                                    b.add_to_wire(&mut pool, new_wir.clone(), new_edg.clone());

                                    let (vf2, vl2) = topexp_vertices(&gvf);

                                    // ici bug orientation : voir tspdrft6
                                    if vl2.is_same(&vff) {
                                        b.add_to_wire(
                                            &mut pool,
                                            new_wir.clone(),
                                            with_orientation(&gvf, Orientation::Forward),
                                        );
                                    } else {
                                        b.add_to_wire(
                                            &mut pool,
                                            new_wir.clone(),
                                            with_orientation(&gvf, Orientation::Reversed),
                                        );
                                    }
                                    mapedgadded.insert(shape_key(&gvf));
                                    map_son_s.insert(shape_key(&edg), (edg.clone(), new_wir));
                                    // NewWire est FORWARD
                                }
                            } else if gvl.shape_type() != ShapeType::Edge
                                || mapedgadded.contains(&shape_key(&gvl))
                            {
                                map_son_s
                                    .insert(shape_key(&edg), (edg.clone(), new_edg.clone()));
                            } else {
                                // OCCT L1295-1317.
                                let new_wir = b.make_wire(&mut pool);
                                b.add_to_wire(&mut pool, new_wir.clone(), new_edg.clone());

                                let (vf2, vl2) = topexp_vertices(&gvl);

                                // ici bug orientation : voir tspdrft6
                                if vf2.is_same(&vll) {
                                    b.add_to_wire(
                                        &mut pool,
                                        new_wir.clone(),
                                        with_orientation(&gvl, Orientation::Forward),
                                    );
                                } else {
                                    b.add_to_wire(
                                        &mut pool,
                                        new_wir.clone(),
                                        with_orientation(&gvl, Orientation::Reversed),
                                    );
                                }
                                mapedgadded.insert(shape_key(&gvl));
                                map_son_s.insert(shape_key(&edg), (edg.clone(), new_wir));
                                // NewWire est FORWARD
                            }

                            // OCCT L1321-1329.
                            if new_edg.orientation == Orientation::Reversed {
                                new_edg.orientation = top_abs_reverse(edg.orientation);
                            } else {
                                new_edg.orientation = edg.orientation;
                            }
                            b.add_to_wire(&mut pool, new_wire_on_f.clone(), new_edg);
                        }
                    }
                    exp.next();
                }
                // OCCT L1333.
                let oriented = with_orientation(&new_wire_on_f, wir.orientation);
                builder_add_face_wire(&mut pool, &mut drft_face, oriented);
                exp3.next();
            }
            thefaces.push(drft_face);
            let _ = &mut fac;
        }

        // OCCT L1338-1348.
        let mut the_subs = BRepToolsSubstitution::new();
        for (key, (key_shape, value)) in map_son_s.iter() {
            let _ = key;
            let mut lsubs: Vec<Shape> = Vec::new();
            exp.init(value, ShapeType::Edge, ShapeType::Shape);
            while exp.more() {
                lsubs.push(exp.current());
                exp.next();
            }
            the_subs.substitute(key_shape, lsubs);
        }

        // on reconstruit les faces
        // OCCT L1350-1358.
        exp.init(&res, ShapeType::Face, ShapeType::Shape);
        while exp.more() {
            let cur = exp.current();
            if contains(self.my_map_f_list(&fkey), &cur) {
                exp.next();
                continue;
            }
            the_subs.build(&cur);
            exp.next();
        }

        // Stockage des descendants
        // OCCT L1360-1418.
        for mi in 0..self.my_map.len() {
            let (key_shape, value) = {
                let entry = self.my_map.get_index(mi).expect("myMap entry");
                (entry.1.0.clone(), entry.1.1.clone())
            };
            if key_shape.shape_type() == ShapeType::Edge {
                // OCCT L1369-1379.
                let mut thedesc: Vec<Shape> = Vec::new();
                the_map.clear();
                for itl in &value {
                    let mapped = map_w
                        .get(&shape_key(itl))
                        .expect("NoSuchObject: MapW")
                        .clone();
                    if the_map.insert(shape_key(&mapped)) {
                        thedesc.push(mapped);
                    }
                }
                self.my_map.get_index_mut(mi).expect("myMap entry").1.1 = thedesc;
            } else if key_shape.is_same(the_f) {
                // OCCT L1380-1387.
                self.my_map.get_index_mut(mi).expect("myMap entry").1.1.clear();
                for itl in &thefaces {
                    self.my_map
                        .get_index_mut(mi)
                        .expect("myMap entry")
                        .1
                        .1
                        .push(itl.clone());
                }
            } else {
                // OCCT L1388-1417.
                let mut thedesc: Vec<Shape> = Vec::new();
                the_map.clear();
                for itl in &value {
                    if the_subs.is_copied(itl) {
                        let copy = the_subs.copy(itl).expect("NoSuchObject: theSubs.Copy");
                        if copy.len() != 1 {
                            println!("Invalid number of descendant");
                            return;
                        } else {
                            let first = copy[0].clone();
                            if the_map.insert(shape_key(&first)) {
                                thedesc.push(first);
                            }
                        }
                    } else if the_map.insert(shape_key(itl)) {
                        thedesc.push(itl.clone());
                    }
                }
                self.my_map.get_index_mut(mi).expect("myMap entry").1.1 = thedesc;
            }
        }

        // OCCT L1420-1433.
        the_map.clear();
        thefaces.clear();
        for (_, (_, value)) in self.my_map.iter() {
            for itl in value {
                if itl.shape_type() == ShapeType::Face && the_map.insert(shape_key(itl)) {
                    thefaces.push(itl.clone());
                }
            }
        }
        let bs = crate::feat::loc_ope_build_shape::LocOpeBuildShape::with_faces(&thefaces);
        self.my_result = bs.shape().cloned().unwrap_or_else(Shape::null);
    }

    /// OCCT myMap(F) — the bound descendant list of the face argument (the
    /// OCCT DataMap operator() raises NoSuchObject when absent).
    fn my_map_f_list(&self, fkey: &(u64, u32)) -> &[Shape] {
        &self
            .my_map
            .get(fkey)
            .expect("NoSuchObject: myMap(F)")
            .1
    }

    /// OCCT LocOpe_SplitDrafts::IsDone() (hxx L96).
    pub fn is_done(&self) -> bool {
        !self.my_result.is_null()
    }

    /// OCCT LocOpe_SplitDrafts::OriginalShape() (hxx L98).
    pub fn original_shape(&self) -> &Shape {
        &self.my_shape
    }

    /// OCCT LocOpe_SplitDrafts::Shape() (cxx L1438-1445).
    pub fn shape(&self) -> &Shape {
        if self.my_result.is_null() {
            panic!("StdFail_NotDone");
        }
        &self.my_result
    }

    /// OCCT LocOpe_SplitDrafts::ShapesFromShape(S) (cxx L1449-1457).
    pub fn shapes_from_shape(&self, the_s: &Shape) -> &[Shape] {
        if self.my_result.is_null() {
            panic!("StdFail_NotDone");
        }
        &self
            .my_map
            .get(&shape_key(the_s))
            .expect("NoSuchObject: myMap(S)")
            .1
    }
}

/// OCCT B.MakeVertex(V, P, Precision::Confusion()) — always a new TShape
/// (architecture difference #6: no position dedup).
fn make_vertex_confusion(pool: &mut BRep, the_p: DVec3) -> Shape {
    let v = pool.add_tvertex_unique(the_p);
    pool.vertex_mut(v.clone()).tolerance = CONFUSION;
    v
}

/// OCCT B.MakeFace(F, S, Tol) — a wire-less face on the surface (the kernel
/// face cache keys on the outer wire, so a fresh marker wire keeps the face
/// unique — architecture difference #6).
fn make_face_from_surface(pool: &mut BRep, the_s: Option<Surface3>, the_tol: f64) -> Shape {
    let marker = pool.add_twire(Vec::new());
    pool.add_tface_tol(the_s, marker, Vec::new(), None, None, Vec::new(), false, the_tol)
}

/// OCCT C->Transformed(Loc.Transformation()) — an identity in the rcad feat
/// flows (architecture difference #12).
fn geom_transformed_identity(the_c: &Curve3) -> Curve3 {
    the_c.clone()
}
