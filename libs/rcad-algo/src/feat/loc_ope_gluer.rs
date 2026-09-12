// OCCT LocOpe_Gluer.hxx L33-81 + LocOpe_Gluer.cxx L47-467 +
// LocOpe_Gluer.lxx L21-90 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Gluer.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Gluer.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Gluer.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// The Perform (cxx L156-334) / AddEdges (cxx L471-556) bodies are translated
// below; all their OCCT dependencies are already translated in this crate:
// LocOpe_WiresOnShape (loc_ope_wires_on_shape.rs), LocOpe_GluedShape
// (loc_ope_glued_shape.rs), LocOpe_Spliter (loc_ope_spliter.rs),
// LocOpe_Generator (loc_ope_generator.rs) and LocOpe::TgtFaces
// (loc_ope.rs::tgt_faces).
//
// The single remaining stand-in is BRepExtrema_ExtPF inside AddEdges (see the
// note at the AddEdges site).
//
// Architecture differences (referenced from the affected functions):
// 1. BRep_Tool::Surface / CurveOnSurface — the face surface and the edge
//    pcurve are carried on the rcad TShape (TFaceData.surface /
//    TEdgeData.pcurves keyed by (face ptr, face location)); standalone feat
//    shapes carry identity locations (same reduction as
//    loc_ope_find_edges.rs).
// 2. GeomAdaptor_Surface(GAS) + FirstU/LastU/FirstV/LastV map to
//    Surface3::default_domain(); UResolution/VResolution to
//    rcad_kernel::topo::topods::{u_resolution_for_surface,
//    v_resolution_for_surface}.
// 3. Extrema_ExtPS — rcad_kernel::base::extrema::ExtPS; the OCCT
//    Initialize(...)+SetFlag(Extrema_ExtFlag_MIN) pair collapses into the
//    rcad constructor + Perform (the rcad solver is minimization-only).
// 4. NCollection_IndexedDataMap / NCollection_DataMap of shapes map to
//    insertion-ordered IndexMap / HashMap keyed by (TShape ptr, Location)
//    — the TopTools_ShapeMapHasher identity; the myMapEF iteration in the
//    Perform follows the insertion order. myMapEE is a HashMap carrying its
//    key Shape with the value (the OCCT itm.Key() of cxx L193 is read); the
//    OCCT bucket iteration order is not reproduced (the OcctShapeMap
//    reduction of brep_feat_builder.rs).
// 5. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate);
//    TopExp::MapShapes is brep_feat_builder::map_shapes and
//    TopExp::MapShapesAndAncestors is the local re-host below.
// 6. BRep_Builder::Continuity (cxx L292-293 / L323-324) -> the
//    BRepSweepBRepBuilder::continuity re-host (brep_sweep/brep_sweep_builder.rs,
//    the BRep_Builder payload form: a standalone LocOpe_Gluer result is not in
//    a BRep pool); BRep_Tool::Continuity is the local read-back below.
//
// first consumer: BRepFeat_Form family (3b) — BRepFeat_Gluer wraps
// LocOpe_Gluer for the glue operation of BRepFeat_Form.

use crate::brep_sweep::brep_sweep_builder::BRepSweepBRepBuilder;
use crate::feat::brep_feat_builder::{explorer, map_shapes, OcctShapeMap};
use crate::feat::brep_feat_make_linear_form::BRepExtremaExtPFCarrier;
use crate::feat::loc_ope::tgt_faces;
use crate::feat::loc_ope_generated_shape::LocOpeGeneratedShape;
use crate::feat::loc_ope_generator::LocOpeGenerator;
use crate::feat::loc_ope_glued_shape::LocOpeGluedShape;
use crate::feat::loc_ope_spliter::LocOpeSpliter;
use crate::feat::loc_ope_wires_on_shape::LocOpeWiresOnShape;
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_tolerance;
use indexmap::IndexMap;
use rcad_kernel::geom::{ Curve2d, Surface3 };
use rcad_kernel::SurfaceEval;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{
    CurveRepresentation, GeomAbsShape, Orientation, ShapeType, TShape,
};
use std::collections::HashMap;

use crate::feat::loc_ope_operation::LocOpeOperation;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
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

/// OCCT TopoDS_Face::Oriented(theOr) — a copy carrying the orientation.
fn face_oriented(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT BRep_Tool::Surface(F) (architecture difference #1).
fn brep_tool_surface(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edg, face, f, l) (architecture
/// difference #1) — the pcurve of the edge on the face with its range.
fn brep_tool_curve_on_surface(edg: &Shape, face: &Shape) -> Option<(Curve2d, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .pcurves
            .get(&shape_key(face))
            .map(|(c, f, l)| (c.clone(), *f, *l)),
        _ => None,
    }
}

/// OCCT static GetOrientation(Fn, Fb) (cxx L354-452).
fn get_orientation(the_fn: &Shape, the_fb: &Shape) -> Orientation {
    // OCCT cxx L357-359: Sn = BRep_Tool::Surface(Fn); Sb = Surface(Fb).
    let Some(sn) = brep_tool_surface(the_fn) else {
        return Orientation::Internal;
    };
    let Some(sb) = brep_tool_surface(the_fb) else {
        return Orientation::Internal;
    };

    // OCCT cxx L368-380: the extrema projector initialization
    // (architecture differences #2/#3).
    let tol_u = rcad_kernel::topo::topods::u_resolution_for_surface(&sb, rcad_kernel::precision::CONFUSION);
    let tol_v = rcad_kernel::topo::topods::v_resolution_for_surface(&sb, rcad_kernel::precision::CONFUSION);
    let sb_domain = sb.default_domain();

    // OCCT cxx L382-450.
    for edg in explorer(the_fn, ShapeType::Edge, ShapeType::Shape) {
        // OCCT cxx L384-385: C2d = BRep_Tool::CurveOnSurface(edg, Fn, f, l).
        let Some((c2d, mut f, mut l)) = brep_tool_curve_on_surface(&edg, the_fn) else {
            continue;
        };
        // OCCT cxx L386-398: infinite range clamping.
        let inf = rcad_kernel::precision::INFINITE_VALUE;
        if f <= -inf && l >= inf {
            f = -100.0;
            l = 100.0;
        } else if f <= -inf {
            f = l - 200.0;
        } else if l >= inf {
            l = f + 200.0;
        }
        let deltau = (l - f) / 20.0;
        use rcad_kernel::geom::{Curve2dEval, SurfaceEval};
        for i in 1..=21i32 {
            // OCCT cxx L402: C2d->D0(f + (i-1)*deltau, ptvtx).
            let ptvtx = c2d.point_at(f + (i - 1) as f64 * deltau);
            // OCCT cxx L403: Sn->D1(ptvtx.X(), ptvtx.Y(), pvt, d1u, d1v).
            let (pvt, d1u, d1v) = sn.derivatives(ptvtx.x, ptvtx.y);
            let mut n1 = d1u.cross(d1v);
            // OCCT cxx L405-411.
            if n1.length() > rcad_kernel::precision::CONFUSION {
                n1 = n1.normalize();
                if the_fn.orientation == Orientation::Reversed {
                    n1 = -n1;
                }

                // OCCT cxx L413-426: anExtPS.Perform(pvt); jmin search.
                let mut an_ext_ps =
                    rcad_kernel::base::extrema::ExtPS::new(pvt, &sb, tol_u, tol_v);
                an_ext_ps.perform(
                    pvt,
                    &sb,
                    sb_domain[0],
                    sb_domain[1],
                    sb_domain[2],
                    sb_domain[3],
                    tol_u,
                    tol_v,
                );
                if an_ext_ps.is_done() {
                    let mut dist2min = f64::MAX; // RealLast()
                    let mut jmin: usize = 0;
                    for j in 1..=an_ext_ps.nb_ext() {
                        if an_ext_ps.square_distance(j) < dist2min {
                            jmin = j;
                            dist2min = an_ext_ps.square_distance(j);
                        }
                    }
                    // OCCT cxx L427-445.
                    if jmin != 0 {
                        let pt = an_ext_ps.point(jmin);
                        let (uu, vv) = (pt.u, pt.v);
                        let (_, d1u, d1v) = sb.derivatives(uu, vv);
                        let mut n2 = d1u.cross(d1v);
                        if n2.length() > rcad_kernel::precision::CONFUSION {
                            n2 = n2.normalize();
                            if the_fb.orientation == Orientation::Reversed {
                                n2 = -n2;
                            }
                            if n1.dot(n2) > 0.0 {
                                return Orientation::Reversed;
                            }
                            return Orientation::Forward;
                        }
                    }
                }
            }
        }
    }
    // OCCT cxx L451.
    Orientation::Internal
}

/// OCCT static Contains(L, S) (cxx L456-467).
fn contains(the_l: &[Shape], the_s: &Shape) -> bool {
    for it in the_l {
        if shape_key(it) == shape_key(the_s) {
            return true;
        }
    }
    false
}

/// OCCT TopExp::MapShapesAndAncestors(S, TS, TA, M) (TopExp.cxx L80-120):
/// M[TS sub-shape] = the ordered list of its TA ancestors, then the TS
/// sub-shapes with no ancestor get an empty entry (L107-117).
///
/// Architecture difference: the same local re-host shape as
/// bop/algo/section.rs and brep_algo/loop.rs (no shared module).
fn map_shapes_and_ancestors(
    the_s: &Shape,
    the_ts: ShapeType,
    the_ta: ShapeType,
    the_m: &mut IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
) {
    // OCCT L87-105: visit the ancestors — for each TA ancestor, map its TS
    // sub-shapes to it.
    for a_anc in explorer(the_s, the_ta, ShapeType::Shape) {
        for a_exs in explorer(&a_anc, the_ts, ShapeType::Shape) {
            let key = shape_key(&a_exs);
            let entry = the_m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    // OCCT L107-117: visit the TS sub-shapes that have no TA ancestor.
    for a_ex in explorer(the_s, the_ts, the_ta) {
        the_m.entry(shape_key(&a_ex)).or_insert((a_ex, Vec::new()));
    }
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) (BRep_Tool.cxx L1180-1188) ->
/// Continuity(E, S1, S2, L1, L2) (BRep_Tool.cxx L1223-1246): the
/// BRep_CurveOn2Surfaces regularity record of the edge matching the two face
/// surfaces; GeomAbs_C0 when no record matches (L1245).
///
/// Architecture difference: OCCT composes the record locations from
/// L.Predivided(E.Location()); the standalone feat shapes travel with
/// identity locations, so the location argument is the face location — the
/// same reduction as the BRepSweepBRepBuilder::continuity writer
/// (brep_sweep/brep_sweep_builder.rs) that stores the record.
fn brep_tool_continuity(the_edg: &Shape, the_f1: &Shape, the_f2: &Shape) -> GeomAbsShape {
    // OCCT L1184-1186: S1 = Surface(F1, l1); S2 = Surface(F2, l2).
    let (Some(s1), Some(s2)) = (brep_tool_surface(the_f1), brep_tool_surface(the_f2)) else {
        // the OCCT null-surface case has no matching representation.
        return GeomAbsShape::C0;
    };
    let ed = match the_edg.data.as_ref() {
        TShape::Edge(ed) => ed,
        _ => return GeomAbsShape::C0,
    };
    // OCCT L1236-1244: the representation scan.
    for cr in &ed.representations {
        // OCCT L1239: cr->IsRegularity(S1, S2, l1, l2).
        if cr.is_regularity_on(&s1, &s2, the_f1.location, the_f2.location) {
            // OCCT L1241: return cr->Continuity().
            if let CurveRepresentation::CurveOn2Surfaces { continuity, .. } = cr {
                return *continuity;
            }
        }
    }
    // OCCT L1245: the GeomAbs_C0 default.
    GeomAbsShape::C0
}

/// OCCT TopoDS_Shape::IsNull() for an entry taken out of a face list.
///
/// rcad: `Shape::null()` carries the null marker (the Vertex stub TShape).
/// A pool-less builder face also has `index == usize::MAX`, so the index
/// test alone would drop legitimate faces (see the feat/ pitfall recorded on
/// `Shape::is_null()`); the TShape-kind test is exact in this position
/// because a descendant-face entry is never a Vertex stub.
fn shape_is_null(the_s: &Shape) -> bool {
    the_s.is_null() && matches!(the_s.data.as_ref(), TShape::Vertex(_))
}

/// OCCT LocOpe_Gluer (LocOpe_Gluer.hxx L33-81).
pub struct LocOpeGluer {
    my_done: bool,        // OCCT: myDone
    my_sb: Shape,         // OCCT: mySb
    my_sn: Shape,         // OCCT: mySn
    my_res: Option<Shape>, // OCCT: myRes (None = null)
    my_ori: Orientation,  // OCCT: myOri
    my_ope: LocOpeOperation, // OCCT: myOpe
    // OCCT: myMapEF (NCollection_IndexedDataMap<Shape, Shape>) — values are
    // Option to carry the Nullify() of cxx L135 (arch. diff. #4).
    my_map_ef: IndexMap<(u64, u32), (Shape, Option<Shape>)>,
    // OCCT: myMapEE (NCollection_DataMap<Shape, Shape>) — the key Shape is
    // carried with the value (the OCCT itm.Key() of cxx L193 is read;
    // arch. diff. #4).
    my_map_ee: HashMap<(u64, u32), (Shape, Shape)>,
    // OCCT: myDescF (NCollection_DataMap<Shape, NCollection_List<Shape>>).
    my_desc_f: HashMap<(u64, u32), Vec<Shape>>,
    my_edges: Vec<Shape>,    // OCCT: myEdges
    my_tgt_edges: Vec<Shape>, // OCCT: myTgtEdges
    // OCCT static `nullList` of DescendantFaces (cxx L348) — the shared
    // empty list the OCCT static carries (arch. diff. #4 vehicle).
    my_null_list: Vec<Shape>,
}

impl Default for LocOpeGluer {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeGluer {
    /// OCCT LocOpe_Gluer::LocOpe_Gluer() (lxx L21-26).
    pub fn new() -> Self {
        LocOpeGluer {
            my_done: false,
            my_sb: Shape::null(),
            my_sn: Shape::null(),
            my_res: None,
            my_ori: Orientation::Internal,
            my_ope: LocOpeOperation::Invalid,
            my_map_ef: IndexMap::new(),
            my_map_ee: HashMap::new(),
            my_desc_f: HashMap::new(),
            my_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
            my_null_list: Vec::new(),
        }
    }

    /// OCCT LocOpe_Gluer::LocOpe_Gluer(Sbase, Snew) (lxx L30-37).
    pub fn with_shapes(the_sbase: &Shape, the_snew: &Shape) -> Self {
        LocOpeGluer {
            my_done: false,
            my_sb: the_sbase.clone(),
            my_sn: the_snew.clone(),
            my_res: None,
            my_ori: Orientation::Internal,
            my_ope: LocOpeOperation::Invalid,
            my_map_ef: IndexMap::new(),
            my_map_ee: HashMap::new(),
            my_desc_f: HashMap::new(),
            my_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
            my_null_list: Vec::new(),
        }
    }

    /// OCCT LocOpe_Gluer::Init(Sbase, Snew) (cxx L53-63).
    pub fn init(&mut self, the_sbase: &Shape, the_snew: &Shape) {
        self.my_sb = the_sbase.clone();
        self.my_sn = the_snew.clone();
        self.my_map_ef.clear();
        self.my_map_ee.clear();
        self.my_desc_f.clear();
        self.my_done = false;
        self.my_ori = Orientation::Internal;
        self.my_ope = LocOpeOperation::Invalid;
    }

    /// OCCT LocOpe_Gluer::Bind(Fnew, Fbase) (cxx L67-141) — the face-face
    /// overload.
    pub fn bind_face(&mut self, the_fnew: &Shape, the_fbase: &Shape) {
        // OCCT cxx L69-81: find Fnew among the faces of mySn.
        let mut fnew_ori: Option<Orientation> = None;
        for cur in explorer(&self.my_sn, ShapeType::Face, ShapeType::Shape) {
            if shape_key(&cur) == shape_key(the_fnew) {
                fnew_ori = Some(cur.orientation);
                break;
            }
        }
        let Some(fnew_ori) = fnew_ori else {
            // OCCT cxx L78-81: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        };

        // OCCT cxx L83-84: Fnor = TopoDS::Face(Fnew.Oriented(exp.Current().
        // Orientation())).
        let fnor = face_oriented(the_fnew, fnew_ori);

        // OCCT cxx L87-97: find Fbase among the faces of mySb.
        let mut fbase_ori: Option<Orientation> = None;
        for cur in explorer(&self.my_sb, ShapeType::Face, ShapeType::Shape) {
            if shape_key(&cur) == shape_key(the_fbase) {
                fbase_ori = Some(cur.orientation);
                break;
            }
        }
        let Some(fbase_ori) = fbase_ori else {
            // OCCT cxx L94-97.
            panic!("Standard_ConstructionError");
        };

        // OCCT cxx L99-101: Fbor.
        let fbor = face_oriented(the_fbase, fbase_ori);
        // OCCT cxx L102: Ori = GetOrientation(Fnor, Fbor).
        let ori = get_orientation(&fnor, &fbor);

        // OCCT cxx L104-120.
        if self.my_ori == Orientation::Internal {
            self.my_ori = ori;
            if self.my_ori == Orientation::Reversed {
                // OCCT cxx L109: mySn.Reverse().
                self.my_sn.orientation = top_abs_reverse(self.my_sn.orientation);
                self.my_ope = LocOpeOperation::Cut;
            } else {
                self.my_ope = LocOpeOperation::Fuse;
            }
        } else if ori != Orientation::Forward {
            // OCCT cxx L117-120.
            self.my_ope = LocOpeOperation::Invalid;
        }

        // OCCT cxx L122-138: bind the edges of Fnor on Fbor.
        for edg in explorer(&fnor, ShapeType::Edge, ShapeType::Shape) {
            let key = shape_key(&edg);
            match self.my_map_ef.get(&key) {
                None => {
                    // OCCT cxx L127-129: myMapEF.Add(edg, Fbor).
                    self.my_map_ef.insert(key, (edg.clone(), Some(fbor.clone())));
                }
                Some((_, bound)) => {
                    let different = match bound {
                        Some(b) => shape_key(b) != shape_key(&fbor),
                        None => false,
                    };
                    if different {
                        // OCCT cxx L132-137: ChangeFromKey(edg).Nullify() —
                        // "edg sur 2 faces. a binder avec l'edge commun".
                        self.my_map_ef.get_mut(&key).expect("entry").1 = None;
                    }
                }
            }
        }
        // OCCT cxx L140: myMapEF.Add(Fnor, Fbor).
        self.my_map_ef
            .insert(shape_key(&fnor), (fnor, Some(fbor)));
    }

    /// OCCT LocOpe_Gluer::Bind(Enew, Ebase) (cxx L145-152) — the edge-edge
    /// overload.
    pub fn bind_edge(&mut self, the_enew: &Shape, the_ebase: &Shape) {
        let key = shape_key(the_enew);
        // OCCT cxx L147-150.
        if let Some((_, bound)) = self.my_map_ee.get(&key) {
            if shape_key(bound) != shape_key(the_ebase) {
                panic!("Standard_ConstructionError");
            }
        }
        // OCCT cxx L151.
        self.my_map_ee
            .insert(key, (the_enew.clone(), the_ebase.clone()));
    }

    /// OCCT LocOpe_Gluer::OpeType() (lxx L73-76).
    pub fn ope_type(&self) -> LocOpeOperation {
        self.my_ope
    }

    /// OCCT LocOpe_Gluer::Perform() (cxx L156-334).
    pub fn perform(&mut self) {
        // OCCT cxx L158: int ind;
        let lmap: usize;

        // OCCT cxx L159-162.
        if self.my_done {
            return;
        }
        // OCCT cxx L163-166.
        if self.my_sb.is_null()
            || self.my_sn.is_null()
            || self.my_map_ef.is_empty()
            || self.my_ope == LocOpeOperation::Invalid
        {
            panic!("Standard_ConstructionError");
        }

        // OCCT cxx L168-169.
        let mut the_wons = LocOpeWiresOnShape::new(&self.my_sb);
        let mut the_gs = LocOpeGluedShape::with_shape(&self.my_sn);

        // OCCT cxx L171.
        lmap = self.my_map_ef.len();

        // OCCT cxx L173-188.
        for ind in 1..=lmap {
            // OCCT cxx L175: TopoDS_Shape S = myMapEF.FindKey(ind).
            let s = self
                .my_map_ef
                .get_index(ind - 1)
                .expect("myMapEF entry")
                .1
                 .0
                .clone();
            // OCCT cxx L176.
            if s.shape_type() == ShapeType::Edge {
                // OCCT cxx L178: TopoDS_Shape S2 = myMapEF(ind).
                let s2 = self
                    .my_map_ef
                    .get_index(ind - 1)
                    .expect("myMapEF entry")
                    .1
                     .1
                    .clone();
                // OCCT cxx L179-182: if (!S2.IsNull()).
                if let Some(s2) = s2 {
                    the_wons.bind_edge_face(&s, &s2);
                }
            } else {
                // OCCT cxx L184-187: TopAbs_FACE.
                the_gs.glue_on_face(&s);
            }
        }

        // OCCT cxx L190-194: the myMapEE iterator.
        // The key Shape travels with the value (arch. diff. #4).
        let the_map_ee: Vec<(Shape, Shape)> = self.my_map_ee.values().cloned().collect();
        for (e_new, e_base) in &the_map_ee {
            the_wons.bind_edge_edge(e_new, e_base);
        }

        // OCCT cxx L196.
        the_wons.bind_all();

        // OCCT cxx L198-201.
        if !the_wons.is_done() {
            return;
        }

        // OCCT cxx L203-204.
        let mut the_split = LocOpeSpliter::with_shape(&self.my_sb);
        the_split.perform(&mut the_wons);

        // OCCT cxx L205-208.
        if !the_split.is_done() {
            return;
        }

        // OCCT cxx L209-215: mise a jour des descendants.
        for exp in explorer(&self.my_sb, ShapeType::Face, ShapeType::Shape) {
            // OCCT cxx L214: myDescF.Bind(...) — Bind assigns, and the
            // spliter returns the static empty list for an unbound face.
            let desc = the_split
                .descendant_shapes(&exp)
                .cloned()
                .unwrap_or_default();
            self.my_desc_f.insert(shape_key(&exp), desc);
        }

        // OCCT cxx L217-225.
        for exp in explorer(&self.my_sn, ShapeType::Face, ShapeType::Shape) {
            // OCCT cxx L219-220: an empty list.
            self.my_desc_f.insert(shape_key(&exp), Vec::new());
            // OCCT cxx L221-224.
            if contains(the_gs.oriented_faces(), &exp) {
                self.my_desc_f
                    .get_mut(&shape_key(&exp))
                    .expect("myDescF bound")
                    .push(exp.clone());
            }
        }

        // OCCT cxx L227-228.
        let mut the_gen = LocOpeGenerator::with_shape(
            the_split
                .resulting_shape()
                .expect("OCCT myRes is a non-null result"),
        );
        the_gen.perform(&mut the_gs);

        // OCCT cxx L230-233.
        self.my_done = the_gen.is_done();
        if self.my_done {
            self.my_res = the_gen.resulting_shape().cloned();

            // OCCT cxx L235.
            self.add_edges();

            // OCCT cxx L237-258: mise a jour des descendants. The OCCT
            // iterator mutates the value of the key it currently visits; the
            // rcad HashMap cannot be mutated while iterated, so the keys are
            // snapshotted first (same assignment semantics).
            let keys: Vec<(u64, u32)> = self.my_desc_f.keys().copied().collect();
            for key in keys {
                let mut new_desc: Vec<Shape> = Vec::new();
                // OCCT cxx L243-244: the itl iterator over itd.Value().
                let the_list = self.my_desc_f.get(&key).cloned().unwrap_or_default();
                for itl in &the_list {
                    // OCCT cxx L246-247: the itl2 iterator over
                    // theGen.DescendantFace(Face(itl.Value())).
                    for itl2 in the_gen.descendant_face(itl) {
                        let descface = itl2;
                        // OCCT cxx L251-254: if (!descface.IsNull()).
                        if !shape_is_null(descface) {
                            new_desc.push(descface.clone());
                        }
                    }
                }
                // OCCT cxx L257: myDescF(itd.Key()) = newDesc.
                self.my_desc_f.insert(key, new_desc);
            }
        }

        // recodage des regularites (OCCT cxx L261-265).
        let mut the_map_ef1: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        let mut the_map_ef2: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        // OCCT cxx L264.
        map_shapes_and_ancestors(
            &self.my_sn,
            ShapeType::Edge,
            ShapeType::Face,
            &mut the_map_ef1,
        );
        // OCCT cxx L265: on myRes (the null shape when !myDone yields the
        // same empty map as the rcad None).
        if let Some(the_res) = self.my_res.clone() {
            map_shapes_and_ancestors(
                &the_res,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_map_ef2,
            );
        }

        // OCCT cxx L267-299.
        for ind in 1..=the_map_ef1.len() {
            // OCCT cxx L269-270.
            let (_, (edg, ll)) = the_map_ef1.get_index(ind - 1).expect("theMapEF1 entry");
            let edg = edg.clone();
            let ll = ll.clone();
            if ll.len() == 2 {
                // OCCT cxx L273-274: LL.First() / LL.Last().
                let fac1 = ll[0].clone();
                let fac2 = ll[ll.len() - 1].clone();
                // OCCT cxx L275.
                let the_cont = brep_tool_continuity(&edg, &fac1, &fac2);
                // OCCT cxx L276.
                if the_cont >= GeomAbsShape::G1 {
                    // OCCT cxx L279: ind2 = theMapEF2.FindIndex(edg) (0 = not
                    // found).
                    if let Some(ind2) = the_map_ef2.get_index_of(&shape_key(&edg)) {
                        // OCCT cxx L282-283.
                        let (_, (_, ll2)) = the_map_ef2.get_index(ind2).expect("theMapEF2 entry");
                        let ll2 = ll2.clone();
                        if ll2.len() == 2 {
                            // OCCT cxx L285-286.
                            let ff1 = ll2[0].clone();
                            let ff2 = ll2[ll2.len() - 1].clone();
                            // OCCT cxx L287-294.
                            if (shape_key(&ff1) == shape_key(&fac1)
                                && shape_key(&ff2) == shape_key(&fac2))
                                || (shape_key(&ff1) == shape_key(&fac2)
                                    && shape_key(&ff2) == shape_key(&fac1))
                            {
                                // OCCT cxx L288-289: an empty body.
                            } else {
                                // OCCT cxx L292-293: BRep_Builder B;
                                // B.Continuity(edg, ff1, ff2, thecont).
                                BRepSweepBRepBuilder.continuity(
                                    &edg, &ff1, &ff2, the_cont,
                                );
                            }
                        }
                    }
                }
            }
        }

        // creation de la liste d`edge (OCCT cxx L300-331).
        the_wons.init_edge_iterator();
        while the_wons.more_edge() {
            // OCCT cxx L304.
            let edg = the_wons.edge();
            for ind in 1..=the_map_ef2.len() {
                // OCCT cxx L307.
                let (_, (edg1, l)) = the_map_ef2.get_index(ind - 1).expect("theMapEF2 entry");
                let edg1 = edg1.clone();
                let l = l.clone();
                // OCCT cxx L308.
                if shape_key(&edg1) == shape_key(&edg) {
                    // OCCT cxx L310.
                    self.my_edges.push(edg.clone());
                    // OCCT cxx L311-327: recodage eventuel des regularites
                    // sur cet edge.
                    if l.len() == 2 {
                        // OCCT cxx L315-316.
                        let fac1 = l[0].clone();
                        let fac2 = l[l.len() - 1].clone();
                        // OCCT cxx L317.
                        if tgt_faces(&edg, &fac1, &fac2) {
                            // OCCT cxx L319.
                            self.my_tgt_edges.push(edg.clone());
                            // OCCT cxx L320.
                            let the_cont = brep_tool_continuity(&edg, &fac1, &fac2);
                            // OCCT cxx L321-324.
                            if the_cont < GeomAbsShape::G1 {
                                BRepSweepBRepBuilder.continuity(
                                    &edg,
                                    &fac1,
                                    &fac2,
                                    GeomAbsShape::G1,
                                );
                            }
                        }
                    }
                }
            }
            // OCCT cxx L330.
            the_wons.next_edge();
        }

        // recodage eventuel des regularites sur cet edge (OCCT cxx L333).
    }

    /// OCCT LocOpe_Gluer::IsDone() (lxx L41-44).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_Gluer::ResultingShape() (lxx L48-55).
    pub fn resulting_shape(&self) -> Option<&Shape> {
        if !self.my_done {
            // OCCT lxx L50-53: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.as_ref()
    }

    /// OCCT LocOpe_Gluer::DescendantFaces(F) (cxx L338-350).
    pub fn descendant_faces(&self, the_f: &Shape) -> &Vec<Shape> {
        if !self.my_done {
            // OCCT cxx L340-343: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if let Some(list) = self.my_desc_f.get(&shape_key(the_f)) {
            return list;
        }
        // OCCT cxx L348-349: static NCollection_List nullList.
        &self.my_null_list
    }

    /// OCCT LocOpe_Gluer::BasisShape() (lxx L59-62).
    pub fn basis_shape(&self) -> &Shape {
        &self.my_sb
    }

    /// OCCT LocOpe_Gluer::GluedShape() (lxx L66-69).
    pub fn glued_shape(&self) -> &Shape {
        &self.my_sn
    }

    /// OCCT LocOpe_Gluer::Edges() (lxx L80-83).
    pub fn edges(&self) -> &Vec<Shape> {
        &self.my_edges
    }

    /// OCCT LocOpe_Gluer::TgtEdges() (lxx L87-90).
    pub fn tgt_edges(&self) -> &Vec<Shape> {
        &self.my_tgt_edges
    }

    /// OCCT LocOpe_Gluer::AddEdges() (cxx L471-556).
    ///
    /// GAP carrier: the `BRepExtrema_ExtPF ext;` of cxx L512 is carried by
    /// BRepExtremaExtPFCarrier (feat/brep_feat_make_linear_form.rs). Its
    /// IsDone() is false, i.e. the OCCT failure path of cxx L521-525
    /// (`flag = 0; break;`) is the one taken; the OCCT success branch is only
    /// reachable when BRepExtrema_ExtPF is translated (the two callers of
    /// this function — the OCCT `if (flag == 1) {}` bodies — are empty, so
    /// the observable result is the same).
    fn add_edges(&mut self) {
        // OCCT cxx L473-474: TopExp_Explorer exp, expsb; exp.Init(mySn,
        // TopAbs_EDGE); — exp is re-initialized inside the loops below.
        // OCCT cxx L476: TopLoc_Location Loc — a dead local of the OCCT
        // source (never read).
        // OCCT cxx L478-480: MapV, MapFPrism, MapE; int flag, i.
        let mut map_v = OcctShapeMap::new();
        let mut map_f_prism = OcctShapeMap::new();
        // MapE is never cleared in the OCCT source (a quirk of cxx L478): it
        // accumulates across the myRes faces.
        let mut map_e = OcctShapeMap::new();
        let mut flag: i32;

        // OCCT cxx L482.
        map_shapes(&self.my_sn, ShapeType::Face, &mut map_f_prism);

        // OCCT cxx L484: for (expsb.Init(myRes, TopAbs_FACE); ...). myRes is
        // the null shape when the generator produced none (the explorer then
        // yields nothing).
        let the_res = self.my_res.clone().unwrap_or_else(Shape::null);
        for expsb_current in explorer(&the_res, ShapeType::Face, ShapeType::Shape) {
            // OCCT cxx L486.
            if !map_f_prism.contains(shape_key(&expsb_current)) {
                // OCCT cxx L488-490.
                map_v.clear();
                map_shapes(&expsb_current, ShapeType::Vertex, &mut map_v);
                map_shapes(&expsb_current, ShapeType::Edge, &mut map_e);

                // OCCT cxx L491: for (exp.Init(mySn, TopAbs_EDGE); ...).
                for exp_current in explorer(&self.my_sn, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT cxx L493-497.
                    let e = exp_current;
                    if map_e.contains(shape_key(&e)) {
                        continue;
                    }
                    // OCCT cxx L498.
                    flag = 0;
                    // OCCT cxx L499-507.
                    for v in explorer(&e, ShapeType::Vertex, ShapeType::Shape) {
                        if map_v.contains(shape_key(&v)) {
                            flag = 1;
                        }
                    }
                    // OCCT cxx L508.
                    if flag == 1 {
                        // OCCT cxx L511-513.
                        let mut ext = BRepExtremaExtPFCarrier::new_default();
                        ext.initialize(&expsb_current);
                        // OCCT cxx L514.
                        flag = 0;
                        // OCCT cxx L515-548.
                        for v in explorer(&e, ShapeType::Vertex, ShapeType::Shape) {
                            // OCCT cxx L518.
                            if !map_v.contains(shape_key(&v)) {
                                // OCCT cxx L520.
                                ext.perform(&v, &expsb_current);
                                // OCCT cxx L521-542.
                                if !ext.is_done() || ext.nb_ext() == 0 {
                                    flag = 0;
                                    break;
                                } else {
                                    // OCCT cxx L528: SquareDistance(1).
                                    let mut dist2min = ext.square_distance(1);
                                    // OCCT cxx L529-532.
                                    for i in 2..=ext.nb_ext() {
                                        dist2min = dist2min.min(ext.square_distance(i));
                                    }
                                    // OCCT cxx L533-541.
                                    if dist2min
                                        >= brep_tool_tolerance(&v) * brep_tool_tolerance(&v)
                                    {
                                        flag = 0;
                                        break;
                                    } else {
                                        flag = 1;
                                    }
                                }
                            } else {
                                // OCCT cxx L546.
                                flag = 1;
                            }
                        }
                        // OCCT cxx L549-551: if (flag == 1) { } — the OCCT
                        // body is empty (the appended data is discarded).
                        if flag == 1 {
                            let _ = flag;
                        }
                    }
                }
            }
        }
    }
}
