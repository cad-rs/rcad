// OCCT LocOpe_Gluer.hxx L33-81 + LocOpe_Gluer.cxx L47-467 +
// LocOpe_Gluer.lxx L21-90 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Gluer.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Gluer.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Gluer.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// DEFERRED METHODS (the translation stops where the OCCT body calls
// not-yet-translated classes; close them in the next batch, in order):
// - Perform() (cxx L156-334) — needs LocOpe_WiresOnShape (1623 lines),
//   LocOpe_Spliter, LocOpe_Generator and LocOpe::TgtFaces (LocOpe.cxx).
// - AddEdges() (cxx L471-556) — needs BRepExtrema_ExtPF (its only call site
//   is Perform, which is itself deferred).
// Everything else (struct fields, Init, Bind x2, Perform-independent
// accessors, and the GetOrientation/Contains statics) is translated below.
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
//    deferred Perform follows the insertion order.
// 5. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate).
//
// first consumer: BRepFeat_Form family (3b) — BRepFeat_Gluer wraps
// LocOpe_Gluer for the glue operation of BRepFeat_Form.

use crate::feat::brep_feat_builder::explorer;
use indexmap::IndexMap;
use rcad_kernel::geom::{ Curve2d, Surface3 };
use rcad_kernel::SurfaceEval;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, TShape };
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
    // OCCT: myMapEE (NCollection_DataMap<Shape, Shape>).
    my_map_ee: HashMap<(u64, u32), Shape>,
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
        if let Some(bound) = self.my_map_ee.get(&key) {
            if shape_key(bound) != shape_key(the_ebase) {
                panic!("Standard_ConstructionError");
            }
        }
        // OCCT cxx L151.
        self.my_map_ee.insert(key, the_ebase.clone());
    }

    /// OCCT LocOpe_Gluer::OpeType() (lxx L73-76).
    pub fn ope_type(&self) -> LocOpeOperation {
        self.my_ope
    }

    /// OCCT LocOpe_Gluer::Perform() (cxx L156-334) — DEFERRED: the body
    /// needs LocOpe_WiresOnShape, LocOpe_Spliter, LocOpe_Generator and
    /// LocOpe::TgtFaces (see the header deferral note).

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

    // OCCT LocOpe_Gluer::AddEdges() (cxx L471-556) — DEFERRED: needs
    // BRepExtrema_ExtPF (see the header deferral note).
}
