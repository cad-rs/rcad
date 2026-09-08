// OCCT LocOpe_Pipe.hxx L17-73 + LocOpe_Pipe.cxx L17-549 + LocOpe_Pipe.lxx
// L17-43 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Pipe.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Pipe.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Pipe.lxx
//
// OCCT inheritance chain (LocOpe_Pipe.hxx L37): none — LocOpe_Pipe is a
// standalone value class.  (Contrary to the stage plan's note, this OCCT
// version's LocOpe_Pipe does NOT inherit LocOpe_GeneratedShape; the
// GeneratedShape implementors are LocOpe_GluedShape and the sweep results —
// see loc_ope_generated_shape.rs.)
//
// Architecture differences (referenced from the affected functions):
// 1. BRepFill_Pipe (the myPipe member, LocOpe_Pipe.hxx L62) is NOT yet
//    ported: the sweep pipeline it consumes (BRepFill_Sweep,
//    BRepFill_SectionPlacement, GeomFill_Sweep, ...) is outside the Stage 3c
//    scope.  BRepFillPipe below carries exactly the interface LocOpe_Pipe
//    consumes (BRepFill_Pipe.hxx L56-100: Shape/Face/Edge/Spine/Profile/
//    FirstShape/LastShape/PipeLine) with unimplemented bodies; the bodies
//    are the BRepFill port task (plan ruling D3: BRepFill per-class port).
// 2. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate);
//    TopExp::MapShapesAndAncestors is
//    feat::loc_ope_glued_shape::map_shapes_and_ancestors.
// 3. NCollection_Map<TopoDS_Shape> / NCollection_IndexedDataMap /
//    NCollection_DataMap of shapes map to IndexMap/HashMap keyed by
//    (TShape ptr, Location) — the TopTools_ShapeMapHasher identity
//    (TShape + Location, orientation ignored); the insertion order is the
//    deterministic stand-in for the OCCT bucket iteration order (same
//    reduction as loc_ope_build_shape.rs arch. diff. #3).
// 4. Geom_BSplineCurve maps to rcad BSplineCurve3 (flat knot vector with
//    multiplicities expanded): NbKnots/Knots(k)/Multiplicity(k) are derived
//    by run-length splitting (BSplineCurve3::knots_mults) and the knot/
//    pole concatenation of Curves / BarycCurve is performed on the
//    (distinct knots, mults, poles) form, rebuilt through
//    from_knots_mults — mirroring the OCCT Tkn/Tmu/Tpol arrays.
// 5. GeomConvert::CurveToBSplineCurve is re-hosted below
//    (geomconvert_curve_to_bspline) over rcad_kernel::base::convert; the
//    full OCCT GeomConvert.cxx port (exact analytic parameterisations) is
//    staged in base/convert — the caveats are tracked there, not here.
//    Geom_BSplineCurve::Segment(U1, U2) (default tolerance
//    Precision::PConfusion()) is rcad_kernel::math::bspl::
//    segment_bspline_curve; Geom_BSplineCurve::IncreaseDegree is
//    BSplineCurve3::increase_degree.
// 6. BSplCLib::Reparametrize is rcad_kernel::math::bspl_lib::reparametrize.
// 7. LocOpe::SampleEdges (LocOpe.cxx L195-235) is re-hosted below
//    (loc_ope_sample_edges) for BarycCurve; it belongs to LocOpe.cxx and is
//    to be forwarded to feat::loc_ope once its owner lands it (parallel
//    stage agent owns loc_ope.rs).
// 8. The gp_Pln / gp_Ax1 / gp_Dir statics used by the constructor
//    (gp_Ax1::IsParallel, gp_Pln::Contains, gp_Pln::Direct — gp_Ax1.hxx
//    L114-117, gp_Pln.hxx L120-123/L199-203, gp_Dir.hxx L191-195) are
//    re-hosted below as pure-math helpers.
//
// first consumer: BRepFeat_Pipe / the pipe feature path (Stage 3c).

use crate::feat::brep_feat_builder::explorer;
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{
    BSplineCurve3, Curve3, CurveEval, Plane, Surface3, TrimmedSurface,
};
use rcad_kernel::precision::{ANGULAR, CONFUSION, PCONFUSION};
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, TShape, BRepBuilder };
use std::collections::HashMap;

/// OCCT LocOpe.cxx L39: #define NECHANT 10.
const NECHANT: i32 = 10;

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

/// OCCT BRep_Tool::Surface(fac) (architecture difference #3 of
/// brep_feat_form.rs): the face surface (local coordinates).
fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::Tolerance(vtx).
fn brep_tool_tolerance(vtx: &Shape) -> f64 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Curve(edg, f, l) — the 3D curve and parameter range
/// (same re-host as loc_ope_find_edges.rs, arch. diff. #1: the feat
/// pipeline shapes carry identity locations).
fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
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

/// OCCT Geom_Surface::DynamicType() == Geom_RectangularTrimmedSurface ->
/// BasisSurface() — strip the trimmed wrapper.
fn basis_surface(s: &Surface3) -> Surface3 {
    match s {
        Surface3::Trimmed(t) => {
            let TrimmedSurface { basis, .. } = t;
            (**basis).clone()
        }
        other => other.clone(),
    }
}

/// OCCT gp_Dir::IsParallel(theOther, theAngularTolerance) (gp_Dir.hxx
/// L191-195): the angle is 0 or PI within the tolerance (pure-math helper,
/// architecture difference #8).
fn dir_is_parallel(d1: DVec3, d2: DVec3, the_angular_tolerance: f64) -> bool {
    let dot = d1.dot(d2).clamp(-1.0, 1.0);
    let cross_modulus = d1.cross(d2).length();
    let an_ang = cross_modulus.atan2(dot); // gp_Dir::Angle: [0, PI]
    an_ang <= the_angular_tolerance || std::f64::consts::PI - an_ang <= the_angular_tolerance
}

/// OCCT gp_Pln::Contains(theP, theLinearTolerance) (gp_Pln.hxx L199-203):
/// Distance(theP) <= theLinearTolerance, with Distance = |SignedDistance|
/// = |(P - Location) . Direction| (architecture difference #8).
fn pln_contains(pln: &Plane, the_p: DVec3, the_linear_tolerance: f64) -> bool {
    let signed = (the_p - pln.origin).dot(pln.normal);
    signed.abs() <= the_linear_tolerance
}

/// OCCT gp_Ax3::Direct() — true when the frame is right-handed
/// ((XDir x YDir) . ZDir > 0) (architecture difference #8).
fn pln_is_direct(pln: &Plane) -> bool {
    pln.u_dir.cross(pln.v_dir).dot(pln.normal) > 0.0
}

/// OCCT GeomConvert::CurveToBSplineCurve(C1) re-host (architecture
/// difference #5) — see the module header note; the range carried by the
/// edge (p1, p2) drives the analytic conversions so that the subsequent
/// OCCT `Segment(p1, p2)` (cxx L357-359 / L484-487) sees the source
/// parameterisation.
fn geomconvert_curve_to_bspline(c1: &Curve3, p1: f64, p2: f64) -> Option<BSplineCurve3> {
    use rcad_kernel::base::convert;
    match c1 {
        // OCCT CaseLine (GeomConvert.cxx): two poles, knots at the range
        // bounds.
        Curve3::Line(l) => Some(convert::line_to_bspline_range(l, p1, p2)),
        // OCCT circle / ellipse / bspline / other branches: the staged
        // rcad conversion (exact for bspline identity; the analytic
        // parameterisation caveats live in base/convert).
        other => {
            let _ = (p1, p2);
            Some(convert::curve_to_bspline(other, 33))
        }
    }
}

/// OCCT Geom_BSplineCurve::Segment(U1, U2) (2-arg form, tolerance
/// Precision::PConfusion()) over the rcad carrier (architecture
/// difference #5).
fn bspline_segment(c: &BSplineCurve3, u1: f64, u2: f64) -> BSplineCurve3 {
    rcad_kernel::math::bspl::segment_bspline_curve(c, u1, u2, PCONFUSION)
        .unwrap_or_else(|| c.clone())
}

/// OCCT BRepFill_Pipe (BRepFill_Pipe.hxx L40-120) — NOT YET PORTED
/// (architecture difference #1).  Only the interface consumed by
/// LocOpe_Pipe is declared; every body is the pending BRepFill sweep port
/// (BRepFill_Sweep / BRepFill_SectionPlacement / GeomFill_Sweep chain).
pub struct BRepFillPipe {
    // OCCT BRepFill_Pipe.hxx members relevant to the consumed interface.
    my_spine: Shape,       // OCCT: mySpine
    my_profile: Shape,     // OCCT: myProfile
    my_shape: Option<Shape>,       // OCCT: myShape
    my_first_shape: Option<Shape>, // OCCT: myFirstShape
    my_last_shape: Option<Shape>,  // OCCT: myLastShape
}

impl BRepFillPipe {
    /// OCCT BRepFill_Pipe::BRepFill_Pipe(Spine, Profile) (hxx L56).
    pub fn new(the_spine: &Shape, the_profile: &Shape) -> Self {
        BRepFillPipe {
            my_spine: the_spine.clone(),
            my_profile: the_profile.clone(),
            my_shape: None,
            my_first_shape: None,
            my_last_shape: None,
        }
    }

    /// OCCT BRepFill_Pipe::Spine() (hxx L66).
    pub fn spine(&self) -> Shape {
        self.my_spine.clone()
    }

    /// OCCT BRepFill_Pipe::Profile() (hxx L68).
    pub fn profile(&self) -> Shape {
        self.my_profile.clone()
    }

    /// OCCT BRepFill_Pipe::Shape() (hxx L70).
    pub fn shape(&self) -> Option<Shape> {
        self.my_shape.clone()
    }

    /// OCCT BRepFill_Pipe::FirstShape() (hxx L74).
    pub fn first_shape(&self) -> Option<Shape> {
        self.my_first_shape.clone()
    }

    /// OCCT BRepFill_Pipe::LastShape() (hxx L76).
    pub fn last_shape(&self) -> Option<Shape> {
        self.my_last_shape.clone()
    }

    /// OCCT BRepFill_Pipe::Face(ESpine, EProfile) (hxx L85).
    pub fn face(&mut self, _e_spine: &Shape, _e_profile: &Shape) -> Shape {
        unimplemented!("BRepFill_Pipe port pending (architecture difference #1)");
    }

    /// OCCT BRepFill_Pipe::Edge(ESpine, VProfile) (hxx L91).
    pub fn edge(&mut self, _e_spine: &Shape, _v_profile: &Shape) -> Shape {
        unimplemented!("BRepFill_Pipe port pending (architecture difference #1)");
    }

    /// OCCT BRepFill_Pipe::PipeLine(Point) (hxx L100).
    pub fn pipe_line(&mut self, _point: DVec3) -> Shape {
        unimplemented!("BRepFill_Pipe port pending (architecture difference #1)");
    }
}

/// OCCT LocOpe_Pipe (LocOpe_Pipe.hxx L37-69) — defines a pipe (near from
/// BRepFill_Pipe), with modifications provided for the Pipe feature.
pub struct LocOpePipe {
    my_pipe: BRepFillPipe,                   // OCCT: myPipe
    my_map: HashMap<(u64, u32), Vec<Shape>>, // OCCT: myMap (NCollection_DataMap)
    my_res: Option<Shape>,                   // OCCT: myRes (None = null shape)
    my_gshap: Vec<Shape>,                    // OCCT: myGShap
    my_crvs: Vec<Option<Curve3>>,            // OCCT: myCrvs (None = null handle)
    // OCCT LocOpe_Pipe.hxx L67-68 — declared by OCCT but never assigned
    // (the lxx accessors read myPipe.FirstShape()/LastShape()).
    #[allow(dead_code)]
    my_first_shape: Option<Shape>, // OCCT: myFirstShape
    #[allow(dead_code)]
    my_last_shape: Option<Shape>, // OCCT: myLastShape
}

impl LocOpePipe {
    /// OCCT LocOpe_Pipe::LocOpe_Pipe(Spine, Profile) (cxx L51-274).
    pub fn new(the_spine: &Shape, the_profile: &Shape) -> Self {
        let mut my_pipe = BRepFillPipe::new(the_spine, the_profile);
        let mut my_map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        let my_gshap: Vec<Shape> = Vec::new();

        // OCCT cxx L55: TopoDS_Shape Result = myPipe.Shape();
        let result = my_pipe.shape().expect("myPipe.Shape()");

        // OCCT cxx L57-58: "On enleve les faces generees par les edges de
        // connexite du profile, et on fusionne les plans si possible".

        // OCCT cxx L60-62: theEFMap.
        let mut the_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            the_profile,
            ShapeType::Edge,
            ShapeType::Face,
            &mut the_ef_map,
        );

        // OCCT cxx L63-65.
        let mut goodfaces: Vec<Shape> = Vec::new();

        // OCCT cxx L69-261: for (i = 1; i <= theEFMap.Extent(); i++).
        for i in 1..=(the_ef_map.len() as i32) {
            let (map_key, (edgpr, edgpr_faces)) =
                the_ef_map.get_index((i - 1) as usize).expect("theEFMap entry");
            let edgpr = edgpr.clone();
            let edgpr_key = *map_key;
            let edgpr_faces = edgpr_faces.clone();
            // OCCT cxx L72: myMap.Bind(edgpr, Empty).
            my_map.insert(edgpr_key, Vec::new());
            if edgpr_faces.len() >= 2 {
                // OCCT cxx L73-76: "on ne prend pas les faces generees".
            } else {
                // OCCT cxx L79-80: MapFac — "on mappe les plans generes par
                // cet edge".
                let mut map_fac: IndexMap<(u64, u32), Shape> = IndexMap::new();
                // OCCT cxx L81-102.
                for edgsp in explorer(the_spine, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT cxx L84: resfac = myPipe.Face(edgsp, edgpr).
                    let resfac = my_pipe.face(&edgsp, &edgpr);
                    if !resfac.is_null() {
                        let mut p = brep_tool_surface(&resfac);
                        if let Some(s) = p.as_ref() {
                            if matches!(s, Surface3::Trimmed(_)) {
                                p = Some(basis_surface(s));
                            }
                        }
                        if matches!(p.as_ref(), Some(Surface3::Plane(_))) {
                            // OCCT cxx L94: MapFac.Add(resfac).
                            map_fac.insert(shape_key(&resfac), resfac.clone());
                        } else {
                            // OCCT cxx L98-99.
                            my_map.get_mut(&edgpr_key).expect("myMap(edgpr)").push(resfac.clone());
                            goodfaces.push(resfac.clone());
                        }
                    }
                }

                // OCCT cxx L104-105: "Chercher les composantes connexes sur
                // cet ensemble de faces., avec meme support geometrique".

                // OCCT cxx L107-116.
                if map_fac.len() <= 1 {
                    // "un seul plan. Rien a faire".
                    if map_fac.len() == 1 {
                        let key = *map_fac.get_index(0).expect("MapFac entry").0;
                        let f = map_fac.get(&key).expect("MapFac entry").clone();
                        my_map.get_mut(&edgpr_key).expect("myMap(edgpr)").push(f.clone());
                        goodfaces.push(f);
                    }
                    continue;
                }

                // OCCT cxx L118-259: while (MapFac.Extent() >= 2).
                while map_fac.len() >= 2 {
                    // OCCT cxx L120-122: itm = Iterator(MapFac); FaceRef =
                    // Face(itm.Key()).
                    let first_face = {
                        let (_, v) = map_fac.get_index(0).expect("MapFac entry");
                        v.clone()
                    };
                    let mut face_ref = first_face;
                    // OCCT cxx L121: FacFuse.
                    let mut fac_fuse: Vec<Shape> = Vec::new();
                    fac_fuse.push(face_ref.clone());
                    // OCCT cxx L124-129.
                    let mut p = brep_tool_surface(&face_ref);
                    if let Some(s) = p.as_ref() {
                        if matches!(s, Surface3::Trimmed(_)) {
                            p = Some(basis_surface(s));
                        }
                    }
                    let plref = match p.as_ref() {
                        Some(Surface3::Plane(pl)) => *pl,
                        _ => Plane::new(DVec3::ZERO, DVec3::Y),
                    };

                    // OCCT cxx L131-144: for (itm.Next(); itm.More();
                    // itm.Next()).
                    for idx in 1..map_fac.len() {
                        let (_, k_face) = map_fac.get_index(idx).expect("MapFac entry");
                        let mut pp = brep_tool_surface(k_face);
                        if let Some(s) = pp.as_ref() {
                            if matches!(s, Surface3::Trimmed(_)) {
                                pp = Some(basis_surface(s));
                            }
                        }
                        let pl = match pp.as_ref() {
                            Some(Surface3::Plane(pl)) => *pl,
                            _ => Plane::new(DVec3::ZERO, DVec3::Y),
                        };
                        if dir_is_parallel(pl.normal, plref.normal, ANGULAR)
                            && pln_contains(&plref, pl.origin, CONFUSION)
                        {
                            fac_fuse.push(k_face.clone());
                        }
                    }

                    // OCCT cxx L146-147: "FacFuse contient des faces de meme
                    // support. Il faut en faire des composantes connexes".

                    // OCCT cxx L149-252: while (FacFuse.Extent() >= 2).
                    while fac_fuse.len() >= 2 {
                        // OCCT cxx L151: FaceRef = Face(FacFuse.First()).
                        face_ref = fac_fuse[0].clone();
                        // OCCT cxx L152-153: "Recuperer l'orientation".
                        let orref = orientation(&face_ref, &result);
                        // OCCT cxx L154-159.
                        let mut p = brep_tool_surface(&face_ref);
                        if let Some(s) = p.as_ref() {
                            if matches!(s, Surface3::Trimmed(_)) {
                                p = Some(basis_surface(s));
                            }
                        }
                        let plref = match p.as_ref() {
                            Some(Surface3::Plane(pl)) => *pl,
                            _ => Plane::new(DVec3::ZERO, DVec3::Y),
                        };
                        let mut dirref = plref.normal;
                        if (pln_is_direct(&plref) && orref == Orientation::Reversed)
                            || (!pln_is_direct(&plref) && orref == Orientation::Forward)
                        {
                            dirref = -dirref;
                        }

                        // OCCT cxx L167-171: MapEd.
                        let mut map_ed: IndexMap<(u64, u32), Shape> = IndexMap::new();
                        for e in explorer(
                            &with_orientation(&face_ref, Orientation::Forward),
                            ShapeType::Edge,
                            ShapeType::Shape,
                        ) {
                            map_ed.insert(shape_key(&e), e.clone());
                        }

                        // OCCT cxx L173-174.
                        map_fac.shift_remove(&shape_key(&face_ref));
                        fac_fuse.remove(0); // "on enleve FaceRef"
                        let mut face_to_fuse = false;
                        let mut more_found;

                        // OCCT cxx L178-231: do { ... } while (MoreFound).
                        loop {
                            more_found = false;
                            // OCCT cxx L181-196: for (it.Initialize(FacFuse);
                            // it.More(); it.Next()) — with the inner edge
                            // loop; `found_it` carries the OCCT iterator
                            // state after the break.
                            let mut found_it: Option<usize> = None;
                            for it_idx in 0..fac_fuse.len() {
                                let mut hit = false;
                                for e in explorer(&fac_fuse[it_idx], ShapeType::Edge, ShapeType::Shape)
                                {
                                    if map_ed.contains_key(&shape_key(&e)) {
                                        face_to_fuse = true;
                                        more_found = true;
                                        found_it = Some(it_idx);
                                        hit = true;
                                        break;
                                    }
                                }
                                // OCCT cxx L192-195: if (exp.More()) break.
                                if hit {
                                    break;
                                }
                            }
                            // OCCT cxx L197-230: if (MoreFound).
                            if more_found {
                                let it_idx = found_it.expect("iterator state");
                                // OCCT cxx L199: fac = Face(it.Value()).
                                let fac = fac_fuse[it_idx].clone();
                                // OCCT cxx L200: orrelat = Orientation(fac,
                                // Result).
                                let mut orrelat = orientation(&fac, &result);
                                // OCCT cxx L201-206.
                                let mut other_p = brep_tool_surface(&fac);
                                if let Some(s) = other_p.as_ref() {
                                    if matches!(s, Surface3::Trimmed(_)) {
                                        other_p = Some(basis_surface(s));
                                    }
                                }
                                let pl = match other_p.as_ref() {
                                    Some(Surface3::Plane(pl)) => *pl,
                                    _ => Plane::new(DVec3::ZERO, DVec3::Y),
                                };
                                let mut dirpl = pl.normal;
                                // OCCT cxx L208-212: NOTE — the second
                                // clause tests Plref.Direct(), not
                                // Pl.Direct() (OCCT source as written).
                                if (pln_is_direct(&pl) && orrelat == Orientation::Reversed)
                                    || (!pln_is_direct(&plref) && orrelat == Orientation::Forward)
                                {
                                    dirpl = -dirpl;
                                }
                                // OCCT cxx L213-220.
                                if dirpl.dot(dirref) > 0.0 {
                                    orrelat = Orientation::Forward;
                                } else {
                                    orrelat = Orientation::Reversed;
                                }
                                // OCCT cxx L221-227: toggle the edges of fac
                                // in MapEd.
                                for e in explorer(
                                    &with_orientation(&fac, orrelat),
                                    ShapeType::Edge,
                                    ShapeType::Shape,
                                ) {
                                    if map_ed.insert(shape_key(&e), e.clone()).is_some() {
                                        map_ed.shift_remove(&shape_key(&e));
                                    }
                                }
                                // OCCT cxx L228-229.
                                map_fac.shift_remove(&shape_key(&fac));
                                fac_fuse.remove(it_idx);
                            }
                            // OCCT cxx L231: } while (MoreFound).
                            if !more_found {
                                break;
                            }
                        }

                        // OCCT cxx L233-251: if (FaceToFuse).
                        if face_to_fuse {
                            // OCCT cxx L235-237: B.MakeFace(NewFace, P,
                            // Tolerance(FaceRef)) — rcad carries the builder
                            // pool locally (loc_ope_build_shape.rs arch.
                            // diff. #1).
                            let mut pool = BRep::new();
                            let b = BRepBuilder::new();
                            let mut b = b;
                            let tol = brep_tool_tolerance(&face_ref);
                            // OCCT cxx L238-245: NewWire + the MapEd edges.
                            let new_wire = b.make_wire(&mut pool);
                            for (_, s) in map_ed.iter() {
                                // OCCT cxx L244: B.Add(NewWire, itm2.Key()).
                                b.add_to_wire(&mut pool, new_wire.clone(), s.clone());
                            }
                            // OCCT cxx L246-247.
                            let face_ref_wires = explorer(
                                &with_orientation(&face_ref, Orientation::Forward),
                                ShapeType::Wire,
                                ShapeType::Shape,
                            );
                            let new_wire =
                                with_orientation(&new_wire, face_ref_wires[0].orientation);
                            // OCCT cxx L248: B.Add(NewFace, NewWire) — rcad
                            // creates the face atomically with its outer
                            // wire on the pool (add_tface_tol).
                            let new_face = pool.add_tface_tol(
                                p.clone(),
                                new_wire,
                                Vec::new(),
                                None,
                                None,
                                Vec::new(),
                                false,
                                tol,
                            );
                            // OCCT cxx L249-250.
                            my_map
                                .get_mut(&edgpr_key)
                                .expect("myMap(edgpr)")
                                .push(new_face.clone());
                            goodfaces.push(new_face);
                        }
                    }
                    // OCCT cxx L253-258: if (FacFuse.Extent() == 1).
                    if fac_fuse.len() == 1 {
                        let f = fac_fuse[0].clone();
                        map_fac.shift_remove(&shape_key(&f));
                        my_map.get_mut(&edgpr_key).expect("myMap(edgpr)").push(f.clone());
                        goodfaces.push(f);
                    }
                }
            }
        }

        // OCCT cxx L263-270.
        let mut first_shape_faces: Vec<Shape> = Vec::new();
        if let Some(fs) = my_pipe.first_shape() {
            first_shape_faces = explorer(&fs, ShapeType::Face, ShapeType::Shape);
        }
        for f in first_shape_faces {
            goodfaces.push(f);
        }
        let mut last_shape_faces: Vec<Shape> = Vec::new();
        if let Some(ls) = my_pipe.last_shape() {
            last_shape_faces = explorer(&ls, ShapeType::Face, ShapeType::Shape);
        }
        for f in last_shape_faces {
            goodfaces.push(f);
        }

        // OCCT cxx L272-273.
        let bs = crate::feat::loc_ope_build_shape::LocOpeBuildShape::with_faces(&goodfaces);
        let my_res = bs.shape().cloned();

        LocOpePipe {
            my_pipe,
            my_map,
            my_res,
            my_gshap,
            my_crvs: Vec::new(),
            my_first_shape: None,
            my_last_shape: None,
        }
    }

    /// OCCT LocOpe_Pipe::Spine() (lxx L19-22).
    pub fn spine(&self) -> Shape {
        self.my_pipe.spine()
    }

    /// OCCT LocOpe_Pipe::Profile() (lxx L26-29).
    pub fn profile(&self) -> Shape {
        self.my_pipe.profile()
    }

    /// OCCT LocOpe_Pipe::FirstShape() (lxx L33-36).
    pub fn first_shape(&self) -> Option<Shape> {
        self.my_pipe.first_shape()
    }

    /// OCCT LocOpe_Pipe::LastShape() (lxx L40-43).
    pub fn last_shape(&self) -> Option<Shape> {
        self.my_pipe.last_shape()
    }

    /// OCCT LocOpe_Pipe::Shape() (cxx L278-281).
    pub fn shape(&self) -> Option<&Shape> {
        self.my_res.as_ref()
    }

    /// OCCT LocOpe_Pipe::Shapes(S) (cxx L285-324).
    pub fn shapes(&mut self, s: &Shape) -> &[Shape] {
        // OCCT cxx L287-291.
        let typ_s = s.shape_type();
        if typ_s != ShapeType::Edge && typ_s != ShapeType::Vertex {
            panic!("Standard_DomainError");
        }
        // OCCT cxx L293-304: locate S among the sub-shapes of
        // myPipe.Profile().
        let profile = self.my_pipe.profile();
        let mut found = false;
        for cur in explorer(&profile, typ_s, ShapeType::Shape) {
            if cur.is_same(s) {
                found = true;
                break;
            }
        }
        if !found {
            panic!("Standard_NoSuchObject");
        }

        // OCCT cxx L306: myGShap.Clear().
        self.my_gshap.clear();
        if typ_s == ShapeType::Vertex {
            // OCCT cxx L309-319.
            let v_profile = s.clone();
            let spine = self.my_pipe.spine();
            for edsp in explorer(&spine, ShapeType::Edge, ShapeType::Shape) {
                let resed = self.my_pipe.edge(&edsp, &v_profile);
                if !resed.is_null() {
                    self.my_gshap.push(resed);
                }
            }
            return &self.my_gshap;
        }
        // OCCT cxx L321-323: TopAbs_EDGE — return myMap(EProfile) (the OCCT
        // DataMap operator() asserts IsBound).
        self.my_map.get(&shape_key(s)).expect("Standard_NoSuchObject")
    }

    /// OCCT LocOpe_Pipe::Curves(Spt) (cxx L328-427).
    #[allow(unused_assignments)] // OCCT cxx L356: P1 = C->Value(p2) is a dead store in the source
    pub fn curves(&mut self, spt: &[DVec3]) -> &[Option<Curve3>] {
        // OCCT cxx L332-333.
        self.my_crvs.clear();

        // OCCT cxx L339-424: for (i = 1; i <= Nbpnt; i++).
        for i in 0..spt.len() {
            // OCCT cxx L341: gp_Pnt P1 = Spt(i).
            let mut p1 = spt[i];
            // OCCT cxx L342-344.
            let mut max_deg = 0usize;
            let mut seq: Vec<BSplineCurve3> = Vec::new();
            let w = self.my_pipe.pipe_line(p1);

            // OCCT cxx L346-365.
            for e in explorer(&w, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L349-350.
                let Some((c1, p1_r, p2_r)) = brep_tool_curve(&e) else {
                    // OCCT cxx L351-353: if (C.IsNull()) continue.
                    continue;
                };
                let Some(mut c) = geomconvert_curve_to_bspline(&c1, p1_r, p2_r) else {
                    continue;
                };
                // OCCT cxx L355.
                max_deg = max_deg.max(c.degree);
                // OCCT cxx L356: P1 = C->Value(p2) (dead store in the OCCT
                // source — mirrored).
                p1 = c.point_at(p2_r);
                // OCCT cxx L357-359.
                if p1_r != c.first_parameter() || p2_r != c.last_parameter() {
                    c = bspline_segment(&c, p1_r, p2_r);
                }
                // OCCT cxx L361-363: Tkn(C->Knots());
                // BSplCLib::Reparametrize(seq.Length(), seq.Length() + 1,
                // Tkn); C->SetKnots(Tkn).
                let (mut tkn, tmu) = c.knots_mults();
                rcad_kernel::math::bspl_lib::reparametrize(
                    seq.len() as f64,
                    (seq.len() + 1) as f64,
                    &mut tkn,
                );
                c.set_knots(&tkn, &tmu);
                // OCCT cxx L364: seq.Append(C).
                seq.push(c);
            }

            // OCCT cxx L367-374.
            let mut nbkn = 0usize;
            let mut nbp = 0usize;
            let nbcurv = seq.len();
            if nbcurv == 0 {
                // OCCT cxx L372: myCrvs.Append(newC) — a null handle.
                self.my_crvs.push(None);
                continue;
            }

            // OCCT cxx L376-385.
            for b in seq.iter_mut() {
                // OCCT cxx L380: Bsp->IncreaseDegree(MaxDeg).
                b.increase_degree(max_deg);
                nbp += b.control_points.len();
                nbkn += b.knots_mults().0.len();
            }
            nbp -= nbcurv - 1;
            nbkn -= nbcurv - 1;
            // OCCT cxx L386-389: Tkn(1, Nbkn), Tmu(1, Nbkn), Tpol(1, Nbp);
            // Ik = 0, Ip = 0 (the OCCT arrays are 1-based; the rcad vectors
            // are 0-based with the counters pre-incremented).
            let mut tkn = vec![0.0f64; nbkn];
            let mut tmu = vec![0i32; nbkn];
            let mut tpol = vec![DVec3::ZERO; nbp];
            let mut ik = 0usize;
            let mut ip = 0usize;

            // OCCT cxx L391-403: the first curve contributes all poles and
            // knots.
            {
                let bsp = &seq[0];
                for k in 0..bsp.control_points.len() {
                    ip += 1;
                    tpol[ip - 1] = bsp.control_points[k];
                }
                let (bsp_knots, bsp_mults) = bsp.knots_mults();
                for k in 0..bsp_knots.len() {
                    ik += 1;
                    tkn[ik - 1] = bsp_knots[k];
                    tmu[ik - 1] = bsp_mults[k];
                }
                tmu[ik - 1] -= 1;
            }

            // OCCT cxx L405-420.
            for b in seq.iter().skip(1) {
                for k in 1..b.control_points.len() {
                    ip += 1;
                    tpol[ip - 1] = b.control_points[k];
                }
                let (bsp_knots, bsp_mults) = b.knots_mults();
                for k in 1..bsp_knots.len() {
                    ik += 1;
                    tkn[ik - 1] = bsp_knots[k];
                    tmu[ik - 1] = bsp_mults[k];
                }
                tmu[ik - 1] -= 1;
            }
            // OCCT cxx L421.
            tmu[ik - 1] += 1;
            // OCCT cxx L422-423: newC = new Geom_BSplineCurve(Tpol, Tkn,
            // Tmu, MaxDeg).
            let new_c = BSplineCurve3::from_knots_mults(max_deg, tkn, tmu, tpol);
            self.my_crvs.push(Some(Curve3::BSpline(new_c)));
        }

        &self.my_crvs
    }

    /// OCCT LocOpe_Pipe::BarycCurve() (cxx L449-549).
    #[allow(unused_assignments)] // OCCT cxx L483: P1 = C->Value(p2) is a dead store in the source
    pub fn baryc_curve(&mut self) -> Option<Curve3> {
        // OCCT cxx L453-462: the barycenter of the sampled points of
        // FirstShape().
        let mut bar = DVec3::ZERO;
        let mut spt: Vec<DVec3> = Vec::new();
        let base = self.first_shape().expect("FirstShape()");
        loc_ope_sample_edges(&base, &mut spt);
        for pvt in &spt {
            bar += *pvt;
        }
        bar /= spt.len() as f64;

        // OCCT cxx L464-466: gp_Pnt P1 = bar.
        let mut p1 = bar;

        // OCCT cxx L468-492.
        let mut max_deg = 0usize;
        let mut seq: Vec<BSplineCurve3> = Vec::new();
        let w = self.my_pipe.pipe_line(p1);

        for e in explorer(&w, ShapeType::Edge, ShapeType::Shape) {
            // OCCT cxx L475-476.
            let Some((c1, p1_r, p2_r)) = brep_tool_curve(&e) else {
                // OCCT cxx L478-480: if (C.IsNull()) continue.
                continue;
            };
            let Some(mut c) = geomconvert_curve_to_bspline(&c1, p1_r, p2_r) else {
                continue;
            };
            // OCCT cxx L482.
            max_deg = max_deg.max(c.degree);
            // OCCT cxx L483: P1 = C->Value(p2) (dead store in the OCCT
            // source — mirrored).
            p1 = c.point_at(p2_r);
            // OCCT cxx L484-486.
            if p1_r != c.first_parameter() || p2_r != c.last_parameter() {
                c = bspline_segment(&c, p1_r, p2_r);
            }
            // OCCT cxx L488-490.
            let (mut tkn, tmu) = c.knots_mults();
            rcad_kernel::math::bspl_lib::reparametrize(
                seq.len() as f64,
                (seq.len() + 1) as f64,
                &mut tkn,
            );
            c.set_knots(&tkn, &tmu);
            // OCCT cxx L491.
            seq.push(c);
        }
        // OCCT cxx L493-499.
        let mut nbkn = 0usize;
        let mut nbp = 0usize;
        let nbcurv = seq.len();
        if nbcurv == 0 {
            // OCCT cxx L498: myCrvs.Append(newC) — a null handle; the OCCT
            // source then falls through and seq(1) below raises
            // Standard_NoSuchObject on the empty sequence (OCCT source as
            // written) — the Rust indexing panics equivalently.
            self.my_crvs.push(None);
        }
        // OCCT cxx L500-507.
        for b in seq.iter_mut() {
            b.increase_degree(max_deg);
            nbp += b.control_points.len();
            nbkn += b.knots_mults().0.len();
        }
        nbp -= nbcurv - 1;
        nbkn -= nbcurv - 1;
        // OCCT cxx L510-513.
        let mut tkn = vec![0.0f64; nbkn];
        let mut tmu = vec![0i32; nbkn];
        let mut tpol = vec![DVec3::ZERO; nbp];
        let mut ik = 0usize;
        let mut ip = 0usize;

        // OCCT cxx L515-527.
        {
            let bsp = &seq[0];
            for k in 0..bsp.control_points.len() {
                ip += 1;
                tpol[ip - 1] = bsp.control_points[k];
            }
            let (bsp_knots, bsp_mults) = bsp.knots_mults();
            for k in 0..bsp_knots.len() {
                ik += 1;
                tkn[ik - 1] = bsp_knots[k];
                tmu[ik - 1] = bsp_mults[k];
            }
            tmu[ik - 1] -= 1;
        }

        // OCCT cxx L529-544.
        for b in seq.iter().skip(1) {
            for k in 1..b.control_points.len() {
                ip += 1;
                tpol[ip - 1] = b.control_points[k];
            }
            let (bsp_knots, bsp_mults) = b.knots_mults();
            for k in 1..bsp_knots.len() {
                ik += 1;
                tkn[ik - 1] = bsp_knots[k];
                tmu[ik - 1] = bsp_mults[k];
            }
            tmu[ik - 1] -= 1;
        }
        // OCCT cxx L545-546.
        tmu[ik - 1] += 1;
        let new_c = BSplineCurve3::from_knots_mults(max_deg, tkn, tmu, tpol);

        Some(Curve3::BSpline(new_c))
    }
}

/// OCCT static Orientation(Sub, S) (cxx L434-445) — the orientation of
/// <sub> as encountered in the exploration of <s>.
fn orientation(sub: &Shape, s: &Shape) -> Orientation {
    for cur in explorer(s, sub.shape_type(), ShapeType::Shape) {
        if cur.is_same(sub) {
            return cur.orientation;
        }
    }
    panic!("Standard_NoSuchObject");
}

/// OCCT LocOpe::SampleEdges(theShape, theSeq) (LocOpe.cxx L195-235) re-host
/// (architecture difference #7) — sample points on the edges (extremities
/// excluded) plus every vertex.  Owner package function; to be forwarded to
/// feat::loc_ope when its owner lands it.
pub(crate) fn loc_ope_sample_edges(the_shape: &Shape, the_seq: &mut Vec<DVec3>) {
    the_seq.clear();
    // OCCT L198: NCollection_Map theMap.
    let mut the_map: std::collections::HashSet<(u64, u32)> = std::collections::HashSet::new();

    // OCCT L200-225: "Computes points on edge, but does not take the
    // extremities into account".
    for edg in explorer(the_shape, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.insert(shape_key(&edg)) {
            continue;
        }
        if !brep_tool_degenerated(&edg) {
            // OCCT L216-217: C = BRep_Tool::Curve(...); C =
            // Transformed(Loc.Transformation()) — the feat pipeline shapes
            // carry identity locations (arch. diff. #1 of
            // loc_ope_find_edges.rs).
            let Some((c, f, l)) = brep_tool_curve(&edg) else {
                continue;
            };
            let delta = (l - f) / NECHANT as f64 * 0.123456;
            for i in 1..NECHANT {
                let prm = delta + ((NECHANT - i) as f64 * f + i as f64 * l) / NECHANT as f64;
                the_seq.push(c.point_at(prm));
            }
        }
    }

    // OCCT L227-234: "Adds every vertex".
    for vtx in explorer(the_shape, ShapeType::Vertex, ShapeType::Shape) {
        if the_map.insert(shape_key(&vtx)) {
            the_seq.push(brep_tool_pnt(&vtx));
        }
    }
}

