// OCCT BRepOffset_Tool.cxx L4422-L4659 — module d of the 1:1 translation
// (the ExtentFace tail: IsInOut / CorrectOrientation / CheckPlanesNormals /
// PerformPlanes / UpdateVertexTolerances; split from brep_offset_tool_c.rs
// to respect the 2000-line file limit; see the architecture-difference list
// in brep_offset_tool.rs — numbering continues there).
//
// This module also carries the BRepOffset_Analyse / BRepOffset_Interval GAP
// carriers shared with the Inter3d module (the OCCT dependency direction
// Inter3d -> Tool is kept).

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topo::topods::{Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_range, brep_tool_tolerance,
};

use super::brep_offset_tool::*;

/// OCCT BRepTopAdaptor_FClass2d(Face, Tol) — the rcad carrier over the
/// topalgo FClass2d + FaceShapeSource (architecture difference #29; the
/// loc_ope_wires_on_shape_b.rs #8 precedent).  The classifier is rebuilt
/// per Perform call (the OCCT class caches the UV polygons; the carrier
/// form keeps the OCCT call structure without the self-referential borrow).
pub(crate) struct FClass2dCarrier {
    face: Shape,
    surf: Surface3,
}

impl FClass2dCarrier {
    /// OCCT BRepTopAdaptor_FClass2d(F, Tol).
    pub fn new(face: &Shape) -> Self {
        let surf = face_surface_of(face)
            .expect("BRepTopAdaptor_FClass2d: the face carries no surface");
        FClass2dCarrier { face: face.clone(), surf }
    }

    /// OCCT BRepTopAdaptor_FClass2d::Perform(Puv).
    pub fn perform(&self, p_uv: DVec2) -> State {
        let locations = [glam::DAffine3::IDENTITY];
        let fss = crate::topalgo::shape_source::FaceShapeSource::new(
            &self.face,
            self.surf.clone(),
            &locations,
        );
        let fc = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
            &fss,
            0,
            rcad_kernel::precision::CONFUSION,
        );
        fc.perform(&fss, p_uv, false)
    }
}



// ---------------------------------------------------------------------------
// OCCT static IsInOut (cxx L4422-4438).
// ---------------------------------------------------------------------------

/// OCCT static IsInOut(FC, AC, S) (cxx L4422-4438) — the pcurve samples are
/// all classified S.  The GCPnts_QuasiUniformDeflection sampling is the
/// reduced uniform form (architecture difference #29); the classifier is
/// the topalgo FClass2d over FaceShapeSource (the loc_ope #8 precedent).
pub(crate) fn is_in_out(
    fc: &FClass2dCarrier,
    ac: &Curve2d,
    f: f64,
    l: f64,
    state: State,
) -> bool {
    // OCCT L4426-4427: Def = 100 * Confusion; QU(AC, Def).
    let n = 23usize;
    for i in 1..=n {
        let t = f + (l - f) * (i as f64) / (n as f64);
        let p = Curve2dEval::point_at(ac, t);
        // OCCT L4431-4434: if (FC.Perform(P) != S) return false.
        if fc.perform(p) != state {
            return false;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::CorrectOrientation (hxx L195-200; cxx L4442-4505).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::CorrectOrientation(SI, NewEdges, AsDes,
/// InitOffset, Offset) (cxx L4442-4505).
pub fn correct_orientation(
    si: &Shape,
    new_edges: &OcctIndexedShapeMap,
    as_des: &mut crate::brep_algo::as_des::BRepAlgoAsDes,
    init_offset: &mut crate::brep_algo::image::BRepAlgoImage,
    offset: f64,
) {
    // OCCT L4450-4452: exp.Init(SI, FACE); f = 0., l = 0.
    let mut f_par = 0.0f64;
    let mut l_par = 0.0f64;

    for fi in explorer(si, ShapeType::Face, ShapeType::Shape) {
        // OCCT L4457-4458: LOF = InitOffset.Image(FI).
        let lof = init_offset.image(&fi);
        for of in &lof {
            // OCCT L4463: LOE = AsDes->ChangeDescendant(OF).
            let has_loe = as_des.has_descendant(of);
            if !has_loe {
                continue;
            }
            // OCCT L4465-4475: the YaInt probe.
            let loe_snapshot = as_des.descendant(of).to_vec();
            let mut ya_int = false;
            for oe in &loe_snapshot {
                if new_edges.contains(oe) {
                    ya_int = true;
                    break;
                }
            }
            if ya_int {
                // OCCT L4478-4479: FC = BRepTopAdaptor_FClass2d(the
                // FORWARD FI, Confusion).
                let mut fi_fwd = fi.clone();
                fi_fwd.orientation = Orientation::Forward;
                let fc = FClass2dCarrier::new(&fi_fwd);

                // OCCT L4482-4501: the orientation walk.
                let loe_snapshot = as_des.descendant(of).to_vec();
                for (pos, oe) in loe_snapshot.iter().enumerate() {
                    if new_edges.contains(oe) {
                        // OCCT L4487-4488: CO2d = CurveOnSurface(OE, OF);
                        // AC(CO2d, f, l).
                        let (co2d, cf, cl) = match brep_tool_curve_on_surface(oe, of) {
                            Some(t) => t,
                            None => continue,
                        };
                        f_par = cf;
                        l_par = cl;

                        if offset > 0.0 {
                            // OCCT L4492-4494: if (IsInOut(FC, AC,
                            // TopAbs_OUT)) OE.Reverse().
                            if is_in_out(&fc, &co2d, f_par, l_par, State::Out) {
                                // OCCT L4484: TopoDS_Shape& OE =
                                // itE.ChangeValue() — the in-place
                                // orientation flip of the descendant entry.
                                if let Some(loe) = as_des.change_descendant(of) {
                                    if let Some(oe_mut) = loe.get_mut(pos) {
                                        oe_mut.orientation = bat::top_abs_reverse(oe_mut.orientation);
                                    }
                                }
                            }
                        }
                        // OCCT L4497-4499: the commented-out else branch.
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::CheckPlanesNormals (hxx L206-208; cxx L4509-4533).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::CheckPlanesNormals(theFace1, theFace2, theTolAng)
/// (cxx L4509-4533) — compares the normal directions of the planar faces.
pub fn check_planes_normals(the_face1: &Shape, the_face2: &Shape, the_tol_ang: f64) -> bool {
    // OCCT L4513-4517: the planar probes.
    let s1 = face_surface_of(the_face1);
    let s2 = face_surface_of(the_face2);
    let (p1, p2) = match (&s1, &s2) {
        (Some(Surface3::Plane(a)), Some(Surface3::Plane(b))) => (a, b),
        _ => return false,
    };
    //
    // OCCT L4519-4523: aDN1 = the plane direction (reversed for REVERSED).
    let mut a_dn1 = p1.normal;
    if the_face1.orientation == Orientation::Reversed {
        a_dn1 = -a_dn1;
    }
    //
    // OCCT L4525-4529: aDN2 (idem).
    let mut a_dn2 = p2.normal;
    if the_face2.orientation == Orientation::Reversed {
        a_dn2 = -a_dn2;
    }
    //
    // OCCT L4531-4532: anAngle = aDN1.Angle(aDN2).
    let an_angle = a_dn1.angle_between(a_dn2);
    an_angle < the_tol_ang
}

// ---------------------------------------------------------------------------
// OCCT static PerformPlanes (cxx L4537-4602).
// ---------------------------------------------------------------------------

/// OCCT static PerformPlanes(theFace1, theFace2, theSide, theL1, theL2)
/// (cxx L4537-4602) — the IntTools_FaceFace plane/plane intersection.  GAP
/// leaf (architecture difference #30): the rcad IntTools FaceFace is
/// DS-bound; the carrier takes the OCCT !IsDone() path (empty result
/// lists).
pub(crate) fn perform_planes(
    the_face1: &Shape,
    the_face2: &Shape,
    the_side: State,
    the_l1: &mut Vec<Shape>,
    the_l2: &mut Vec<Shape>,
) {
    the_l1.clear();
    the_l2.clear();
    // Intersect the planes using IntTools_FaceFace directly
    // OCCT L4546-4548: aFF.SetParameters(true, true, true, Confusion);
    // aFF.Perform(theFace1, theFace2).
    let is_done: bool = false; // GAP: IntTools_FaceFace::Perform (bare-face form)
    //
    // OCCT L4550-4553: if (!aFF.IsDone()) return.
    if !is_done {
        return;
    }
    // OCCT L4555-4601: the single-curve edge build, OrientSection + the
    // Side reversal, the result appends — behind the OCCT !IsDone() path.
    let _ = (the_face1, the_face2, the_side);
}

// ---------------------------------------------------------------------------
// OCCT static UpdateVertexTolerances (cxx L4613-4659).
// ---------------------------------------------------------------------------

/// OCCT static UpdateVertexTolerances(theFace) (cxx L4613-4659).
pub(crate) fn update_vertex_tolerances(the_face: &Shape) {
    // OCCT L4615-4618: VEmap = MapShapesAndAncestors(VERTEX, EDGE).
    let mut vemak: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
        the_face,
        ShapeType::Vertex,
        ShapeType::Edge,
        &mut vemak,
    );

    for i in 1..=shape_indexed_data_map::extent(&vemak) {
        let a_vertex = shape_indexed_data_map::find_key_1(&vemak, i).clone();
        let elist = shape_indexed_data_map::value_1(&vemak, i).clone();
        let pnt_vtx = match bat::brep_tool_pnt(&a_vertex) {
            Some(t) => t,
            None => continue,
        };
        for it_value in &elist {
            let an_edge = it_value;
            let (v1, v2) = top_exp_vertices(an_edge);
            let (fpar, lpar) = brep_tool_range(an_edge);
            // OCCT L4633: aParam = (V1.IsSame(aVertex)) ? fpar : lpar.
            let a_param = if v1.is_same(&a_vertex) { fpar } else { lpar };
            if !brep_tool_degenerated(an_edge) {
                // OCCT L4636-4644: the 3d-curve distance probes.
                if let Some((c, _, _)) = brep_tool_curve(an_edge) {
                    let a_pnt = CurveEval::point_at(&c, a_param);
                    let a_dist = pnt_vtx.distance(a_pnt);
                    let mut vtx = a_vertex.clone();
                    bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    if v1.is_same(&v2) {
                        let a_pnt = CurveEval::point_at(&c, lpar);
                        let a_dist = pnt_vtx.distance(a_pnt);
                        bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    }
                }
            }
            // OCCT L4647-4656: the pcurve-on-face distance probes (the
            // BRepAdaptor_Curve(E, theFace) form) — the rcad form evaluates
            // the pcurve point through the face surface.
            if let Some((c2d, _, _)) = brep_tool_curve_on_surface(an_edge, the_face) {
                if let Some(surf) = face_surface_of(the_face) {
                    let p2d = Curve2dEval::point_at(&c2d, a_param);
                    let a_pnt = surf.point_at(p2d.x, p2d.y);
                    let a_dist = pnt_vtx.distance(a_pnt);
                    let mut vtx = a_vertex.clone();
                    bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    if v1.is_same(&v2) {
                        let p2d = Curve2dEval::point_at(&c2d, lpar);
                        let a_pnt = surf.point_at(p2d.x, p2d.y);
                        let a_dist = pnt_vtx.distance(a_pnt);
                        bat::builder_update_vertex_tol(&mut vtx, a_dist);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// BRepOffset_Analyse / BRepOffset_Interval GAP carriers (shared with the
// Inter3d module; the OCCT dependency direction Inter3d -> Tool keeps the
// carriers here).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Interval (BRepOffset_Interval.hxx L30-70) — GAP carrier
/// (staged; the Type() accessor is the only consumed form of this batch).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BRepOffsetInterval {
    /// OCCT: myType (ChFiDS_TypeOfConcavity).
    pub my_type: crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity,
}

/// OCCT BRepOffset_Analyse (BRepOffset_Analyse.hxx L40-140) — GAP carrier
/// (staged as its own translation unit of the package; the accessor
/// surface is the one consumed by Inter3d and Tool).
#[derive(Debug, Clone)]
pub(crate) struct BRepOffsetAnalyse;

impl BRepOffsetAnalyse {
    /// OCCT BRepOffset_Analyse::Type(E) — the interval list of the edge.
    pub fn type_(&self, _e: &Shape) -> Vec<BRepOffsetInterval> {
        panic!("GAP: BRepOffset_Analyse::Type (BRepOffset_Analyse not translated)");
    }
    /// OCCT BRepOffset_Analyse::Ancestors(S).
    pub fn ancestors(&self, _s: &Shape) -> Vec<Shape> {
        panic!("GAP: BRepOffset_Analyse::Ancestors (BRepOffset_Analyse not translated)");
    }
    /// OCCT BRepOffset_Analyse::HasAncestor(S).
    pub fn has_ancestor(&self, _s: &Shape) -> bool {
        panic!("GAP: BRepOffset_Analyse::HasAncestor (BRepOffset_Analyse not translated)");
    }
    /// OCCT BRepOffset_Analyse::Descendants(S) — the nullable form.
    pub fn descendants(&self, _s: &Shape) -> Option<Vec<Shape>> {
        panic!("GAP: BRepOffset_Analyse::Descendants (BRepOffset_Analyse not translated)");
    }
    /// OCCT BRepOffset_Analyse::Generated(S).
    pub fn generated(&self, _s: &Shape) -> Shape {
        panic!("GAP: BRepOffset_Analyse::Generated (BRepOffset_Analyse not translated)");
    }
    /// OCCT BRepOffset_Analyse::NewFaces().
    pub fn new_faces(&self) -> Vec<Shape> {
        panic!("GAP: BRepOffset_Analyse::NewFaces (BRepOffset_Analyse not translated)");
    }
    /// OCCT BRepOffset_Analyse::TangentEdges(Edge, Vertex, TangOnV).
    pub fn tangent_edges(&self, _edge: &Shape, _vertex: &Shape, _tang_on_v: &mut Vec<Shape>) {
        panic!("GAP: BRepOffset_Analyse::TangentEdges (BRepOffset_Analyse not translated)");
    }
}
