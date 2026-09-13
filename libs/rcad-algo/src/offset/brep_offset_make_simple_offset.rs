// OCCT BRepOffset_MakeSimpleOffset.cxx L1-703 + BRepOffset_MakeSimpleOffset.hxx
// L30-172 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_MakeSimpleOffset.cxx / .hxx
//
// OCCT inheritance chain (hxx L61): none — BRepOffset_MakeSimpleOffset is a
// standalone class.
//
// Architecture differences:
// 1. TopoDS_Shape -> rcad Shape (Arc<TShape> handle + location + orientation);
//    NCollection_DataMap<TopoDS_Vertex, TopoDS_Edge> (myMapVE) -> HashMap
//    keyed by (TShape ptr, Location) — the TopTools_ShapeMapHasher identity
//    (orientation ignored), the same map form as feat/loc_ope_wires_on_shape.
// 2. BRepTools_Modifier (TKTopAlgo/BRepTools) is the real class of
//    topalgo/brep_tools_modifier.rs; BRepOffset_SimpleOffset (TKOffset/BRepOffset)
//    — the rebuild engine's mapper — is the interface implementor below: the
//    six OCCT overrides are the BRepTools_Modification methods, translated 1:1
//    here, together with the FillOffsetData / FillFaceData / FillEdgeData /
//    FillVertexData construction (cxx L168-427).  Everything around the engine
//    calls is translated 1:1.
// 2.1 The OCCT NCollection_DataMap<TopoDS_Face|Edge|Vertex, *Data> with
//    TopTools_ShapeMapHasher maps to HashMap<(TShape ptr, Location), *Data>
//    (the hasher ignores the orientation, so the key is shape_key).
// 2.2 occ::handle<Geom_Curve|Geom_Surface|Geom2d_Curve> members map to
//    Option<...> (None = the OCCT null handle).
// 3. ShapeAnalysis_FreeBounds (TKShHealing), BRepTools_Quilt (TKTopAlgo),
//    ShapeFix_Edge::FixSameParameter (TKShHealing), GeomFill_Generator
//    (TKGeomAlgo) and the planar BRepLib_MakeFace(W, OnlyPlane) constructor
//    have no rcad translation yet — GAP carriers below.  (BRepLib::
//    BuildCurves3d is translated — topalgo/brep_lib/build_curves3d.rs;
//    the call sites call the real body over my_brep.)
// 4. OCCT mutates TShapes in place through a global arena; rcad carries the
//    arena as the BRep pool.  The class holds my_brep (the pool stand-in for
//    the arena; consumed by BRepBuilder mutations and the rcad
//    ShapeBuildReShape API, which takes &mut BRep).  The pool-local stand-in
//    for `BRep_Builder aBB;` locals is a BRepBuilder over that pool.
// 5. BRepLib_MakeEdge(V1, V2) — the rcad BRepBuilder::add_edge construction
//    cannot fail, so the OCCT IsDone() guards (cxx L549/L570) have no rcad
//    counterpart (annotated at the call sites).
// 6. BRepAdaptor_Surface(F, false)::D1 — the rcad stand-in evaluates the face
//    surface directly through Surface3::dn (the GeomAdaptor_Surface::DN
//    vehicle, geom/eval.rs); the BRepAdaptor_Surface wrapper needs the BRep
//    pool and carries no additional state for analytic surfaces.
// 7. The rcad ShapeBuildReShape API takes &mut self and &mut BRep
//    (shhealing/shape_build/reshape.rs); the OCCT const Generated/Modified
//    and the myReShape->Apply calls take the class by &mut self.

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, SurfaceEval, TrimmedCurve3,
};
use rcad_kernel::topo::topods::{tshape_flags, BRep, BRepBuilder, GeomAbsShape, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_wires_on_shape::{brep_tool_pnt, shape_key, top_exp_vertices};
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve, brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_parameter,
    brep_tool_range, brep_tool_surface, brep_tool_tolerance, ShapeKey,
};
use crate::offset::brep_offset_offset::BRepOffsetStatus;
use crate::offset::brep_offset_surface::{brep_offset_surface, collapse_singularities};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::topalgo::brep_lib_validate_edge::{
    Adaptor3dCurveOnSurface, BRepLibValidateEdge, Geom2dAdaptorCurve, GeomAdaptorCurve,
    GeomAdaptorSurface,
};
use crate::topalgo::brep_tools_modification::BRepToolsModification;
use crate::topalgo::brep_tools_modifier::BRepToolsModifier;

// ---------------------------------------------------------------------------
// OCCT BRepOffsetSimple_Status (hxx L30-38).
// ---------------------------------------------------------------------------

/// OCCT BRepOffsetSimple_Status (BRepOffset_MakeSimpleOffset.hxx L30-38).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepOffsetSimpleStatus {
    Ok,
    NullInputShape,
    ErrorOffsetComputation,
    ErrorWallFaceComputation,
    ErrorInvalidNbShells,
    ErrorNonClosedShell,
}

// ---------------------------------------------------------------------------
// OCCT BRepTools_Modifier + BRepOffset_SimpleOffset (architecture difference
// #2).
// ---------------------------------------------------------------------------

/// OCCT NewFaceData (BRepOffset_SimpleOffset.hxx L147-154).
struct NewFaceData {
    /// OCCT: myOffsetS.
    my_offset_s: Option<Surface3>,
    /// OCCT: myL.
    my_l: u32,
    /// OCCT: myTol.
    my_tol: f64,
    /// OCCT: myRevWires.
    my_rev_wires: bool,
    /// OCCT: myRevFace.
    my_rev_face: bool,
}

/// OCCT NewEdgeData (BRepOffset_SimpleOffset.hxx L156-161).
struct NewEdgeData {
    /// OCCT: myOffsetC (Resulting curve).
    my_offset_c: Option<Curve3>,
    /// OCCT: myL.
    my_l: u32,
    /// OCCT: myTol.
    my_tol: f64,
}

/// OCCT NewVertexData (BRepOffset_SimpleOffset.hxx L163-167).
struct NewVertexData {
    /// OCCT: myP.
    my_p: DVec3,
    /// OCCT: myTol.
    my_tol: f64,
}

/// OCCT BRepOffset_SimpleOffset (BRepOffset_SimpleOffset.hxx L44-190;
/// BRepOffset_SimpleOffset.cxx L1-427) — the BRepTools_Modification mapper of
/// the simple offset algorithm (architecture difference #2).
pub struct BRepOffsetSimpleOffset {
    /// OCCT: myFaceInfo (hxx L177) — Map of faces to new faces information.
    my_face_info: HashMap<ShapeKey, NewFaceData>,
    /// OCCT: myEdgeInfo (hxx L180) — Map of edges to new edges information.
    my_edge_info: HashMap<ShapeKey, NewEdgeData>,
    /// OCCT: myVertexInfo (hxx L183) — Map of vertices to new vertices
    /// information.
    my_vertex_info: HashMap<ShapeKey, NewVertexData>,
    /// OCCT: myOffsetValue (hxx L186) — Offset value.
    my_offset_value: f64,
    /// OCCT: myTolerance (hxx L189) — Tolerance.
    my_tolerance: f64,
    /// The rcad arena stand-in for the BRep_Builder locals of FillEdgeData
    /// (architecture difference #4).
    my_brep: BRep,
}

impl BRepOffsetSimpleOffset {
    /// OCCT BRepOffset_SimpleOffset::BRepOffset_SimpleOffset(theInputShape,
    /// theOffsetValue, theTolerance) (BRepOffset_SimpleOffset.cxx L40-47).
    pub fn new(the_input_shape: &Shape, the_offset_value: f64, the_tolerance: f64) -> Self {
        // OCCT L43-44: myOffsetValue(theOffsetValue), myTolerance(theTolerance).
        let mut a_mapper = BRepOffsetSimpleOffset {
            my_face_info: HashMap::new(),
            my_edge_info: HashMap::new(),
            my_vertex_info: HashMap::new(),
            my_offset_value: the_offset_value,
            my_tolerance: the_tolerance,
            my_brep: BRep::new(),
        };
        // OCCT L46: FillOffsetData(theInputShape);
        a_mapper.fill_offset_data(the_input_shape);
        a_mapper
    }

    /// OCCT BRepOffset_SimpleOffset::FillOffsetData (cxx L168-202) — Fills
    /// offset data.
    fn fill_offset_data(&mut self, the_shape: &Shape) {
        // Clears old data.
        // OCCT L171-173: myFaceInfo.Clear(); myEdgeInfo.Clear();
        // myVertexInfo.Clear();
        self.my_face_info.clear();
        self.my_edge_info.clear();
        self.my_vertex_info.clear();

        // Faces loop. Compute offset surface for each face.
        // OCCT L176-181: TopExp_Explorer anExpSF(theShape, TopAbs_FACE);
        for a_curr_face in explorer(the_shape, ShapeType::Face, ShapeType::Shape) {
            // OCCT L179-180: const TopoDS_Face& aCurrFace =
            // TopoDS::Face(anExpSF.Current()); FillFaceData(aCurrFace);
            self.fill_face_data(&a_curr_face);
        }

        // Iterate over edges to compute 3d curve.
        // OCCT L184-186: NCollection_IndexedDataMap<...> aEdgeFaceMap;
        // TopExp::MapShapesAndAncestors(theShape, TopAbs_EDGE, TopAbs_FACE,
        // aEdgeFaceMap);
        let mut a_edge_face_map: indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)> =
            indexmap::IndexMap::new();
        map_shapes_and_ancestors(
            the_shape,
            ShapeType::Edge,
            ShapeType::Face,
            &mut a_edge_face_map,
        );
        // OCCT L187-191: for (int anIdx = 1; anIdx <= aEdgeFaceMap.Length();
        // ++anIdx).
        for an_idx in 1..=a_edge_face_map.len() {
            // OCCT L189: const TopoDS_Edge& aCurrEdge =
            // TopoDS::Edge(aEdgeFaceMap.FindKey(anIdx));
            let a_curr_edge = a_edge_face_map
                .get_index(an_idx - 1)
                .expect("indexed data map index")
                .1
                 .0
                .clone();
            self.fill_edge_data(&a_curr_edge, &a_edge_face_map, an_idx);
        }

        // Iterate over vertices to compute new vertex.
        // OCCT L194-196.
        let mut a_vertex_edge_map: indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)> =
            indexmap::IndexMap::new();
        map_shapes_and_ancestors(
            the_shape,
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut a_vertex_edge_map,
        );
        // OCCT L197-201.
        for an_idx in 1..=a_vertex_edge_map.len() {
            // OCCT L199: const TopoDS_Vertex& aCurrVertex =
            // TopoDS::Vertex(aVertexEdgeMap.FindKey(anIdx));
            let a_curr_vertex = a_vertex_edge_map
                .get_index(an_idx - 1)
                .expect("indexed data map index")
                .1
                 .0
                .clone();
            self.fill_vertex_data(&a_curr_vertex, &a_vertex_edge_map, an_idx);
        }
    }

    /// OCCT BRepOffset_SimpleOffset::FillFaceData (cxx L206-233) — Method to
    /// fill new face data for single face.
    fn fill_face_data(&mut self, the_face: &Shape) {
        // OCCT L208-211: NewFaceData aNFD; aNFD.myRevWires = false;
        // aNFD.myRevFace = false; aNFD.myTol = BRep_Tool::Tolerance(theFace);
        let mut a_nfd = NewFaceData {
            my_offset_s: None,
            my_l: 0,
            my_tol: brep_tool_tolerance(the_face),
            my_rev_wires: false,
            my_rev_face: false,
        };

        // Create offset surface.

        // Any existing transformation is applied to the surface.
        // New face will have null transformation.
        // OCCT L217: occ::handle<Geom_Surface> aS = BRep_Tool::Surface(theFace);
        let a_s = brep_tool_surface(the_face)
            .expect("BRep_Tool::Surface: the face carries no surface");
        // OCCT L218: aS = BRepOffset::CollapseSingularities(aS, theFace,
        // myTolerance);
        let a_s = collapse_singularities(&a_s, the_face, self.my_tolerance);

        // Take into account face orientation.
        // OCCT L221-225.
        let a_mult = if the_face.orientation == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };

        // OCCT L227-228: BRepOffset_Status aStatus; aNFD.myOffsetS =
        // BRepOffset::Surface(aS, aMult * myOffsetValue, aStatus, true);
        let mut a_status = BRepOffsetStatus::Good;
        a_nfd.my_offset_s = Some(brep_offset_surface(
            &a_s,
            a_mult * self.my_offset_value,
            &mut a_status,
            true,
        ));
        // OCCT L229: aNFD.myL = TopLoc_Location(); // Null transformation.
        a_nfd.my_l = 0;

        // Save offset surface in map.
        // OCCT L232: myFaceInfo.Bind(theFace, aNFD);
        self.my_face_info.insert(shape_key(the_face), a_nfd);
    }

    /// OCCT BRepOffset_SimpleOffset::FillEdgeData (cxx L237-314) — Method to
    /// fill new edge data for single edge.
    ///
    /// The `unused_assignments` allowance covers the OCCT aF/aL out-parameters:
    /// `BRep_Tool::Curve(aNewEdge, aNED.myL, aF, aL)` stores the 3d-curve
    /// range, which the very first `BRep_Tool::CurveOnSurface(theEdge,
    /// aCurFace, aF, aL)` of the loop below overwrites (the OCCT dead store is
    /// kept for form).
    #[allow(unused_assignments)]
    fn fill_edge_data(
        &mut self,
        the_edge: &Shape,
        the_edge_face_map: &indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)>,
        the_idx: usize,
    ) {
        // OCCT L244: const NCollection_List<TopoDS_Shape>& aFacesList =
        // theEdgeFaceMap(theIdx);
        let a_faces_list = indexed_data_map_value(the_edge_face_map, the_idx);

        if a_faces_list.is_empty() {
            return; // Free edges are skipped.
        }

        // Get offset surface.
        // OCCT L252: const TopoDS_Face& aCurrFace = TopoDS::Face(aFacesList.First());
        let a_curr_face = a_faces_list[0].clone();

        // OCCT L254-257: if (!myFaceInfo.IsBound(aCurrFace)) return;
        let Some(a_nfd) = self.my_face_info.get(&shape_key(&a_curr_face)) else {
            return;
        };

        // No need to deal with transformation - it is applied in fill faces
        // data method.
        // OCCT L260-261.
        let an_offset_surf = a_nfd.my_offset_s.clone();

        // Compute offset 3d curve.
        // OCCT L264-265: double aF, aL; occ::handle<Geom2d_Curve> aC2d =
        // BRep_Tool::CurveOnSurface(theEdge, aCurrFace, aF, aL);
        let a_c2d = brep_tool_curve_on_surface(the_edge, &a_curr_face);
        let (mut a_f, mut a_l) = match &a_c2d {
            Some((_, f, l)) => (*f, *l),
            None => (0.0, 0.0),
        };

        // OCCT L267-268: BRepBuilderAPI_MakeEdge anEdgeMaker(aC2d,
        // anOffsetSurf, aF, aL); TopoDS_Edge aNewEdge = anEdgeMaker.Edge();
        // The rcad carrier of that maker is the ShapeBuild_Edge::MakeEdge(E,
        // pcurve, S, L, p1, p2) wrapper (ShapeBuild_Edge.cxx L851-881, whose
        // body is exactly `BRepBuilderAPI_MakeEdge ME(pcurve, S, p1, p2); E =
        // ME.Edge();`); the not-done maker yields the null edge, which the
        // rcad form carries as Shape::null() (architecture difference #5).
        let mut a_new_edge = Shape::null();
        if let (Some((a_c2d, _a_f, _a_l)), Some(an_offset_surf)) = (&a_c2d, &an_offset_surf) {
            let a_shape_build_edge = ShapeBuildEdge;
            // OCCT: the maker's TopLoc_Location default is identity (0).
            a_shape_build_edge.make_edge_pcurve_surface_loc_params(
                &mut self.my_brep,
                &mut a_new_edge,
                a_c2d,
                an_offset_surf,
                0,
                a_f,
                a_l,
            );
        }

        // Compute max tolerance. Vertex tolerance usage is taken from existing
        // offset computation algorithm. This piece of code significantly
        // influences resulting performance.
        // OCCT L272: double aTol = BRep_Tool::MaxTolerance(theEdge, TopAbs_VERTEX);
        let a_tol = brep_tool_max_tolerance(the_edge, ShapeType::Vertex);
        // OCCT L273: BRepLib::BuildCurves3d(aNewEdge, aTol);
        crate::topalgo::brep_lib::brep_lib::BRepLib::build_curves3d_tol(
            &mut self.my_brep,
            &a_new_edge,
            a_tol,
        );

        // OCCT L275-276: NewEdgeData aNED; aNED.myOffsetC =
        // BRep_Tool::Curve(aNewEdge, aNED.myL, aF, aL);
        let mut a_ned = NewEdgeData {
            my_offset_c: None,
            my_l: 0,
            my_tol: 0.0,
        };
        match brep_tool_curve(&a_new_edge) {
            Some((a_offset_c, a_curve_f, a_curve_l)) => {
                a_ned.my_offset_c = Some(a_offset_c);
                // OCCT: L = E.Location() * GC->Location().
                a_ned.my_l = a_new_edge.location;
                a_f = a_curve_f;
                a_l = a_curve_l;
            }
            None => {
                // OCCT: L.Identity(); First = Last = 0.;
                a_ned.my_l = 0;
                a_f = 0.0;
                a_l = 0.0;
            }
        }

        // Iterate over adjacent faces for the current edge and compute max
        // deviation.
        // OCCT L279-281: double anEdgeTol = 0.0; NCollection_List<...>::Iterator
        // anIter(aFacesList); for (; !aNED.myOffsetC.IsNull() && anIter.More();
        // anIter.Next())
        let mut an_edge_tol = 0.0f64;
        let mut an_iter = 0usize;
        while a_ned.my_offset_c.is_some() && an_iter < a_faces_list.len() {
            // OCCT L283: const TopoDS_Face& aCurFace = TopoDS::Face(anIter.Value());
            let a_cur_face = a_faces_list[an_iter].clone();
            // OCCT anIter.Next() — placed right after the value read so that
            // the `continue` paths below match the OCCT for-loop increment.
            an_iter += 1;

            // OCCT L285-288: if (!myFaceInfo.IsBound(aCurFace)) continue;
            let Some(a_cur_nfd) = self.my_face_info.get(&shape_key(&a_cur_face)) else {
                continue;
            };

            // Create offset curve on surface.
            // OCCT L291: const occ::handle<Geom2d_Curve> aC2dNew =
            // BRep_Tool::CurveOnSurface(theEdge, aCurFace, aF, aL);
            // The OCCT out-parameters aF/aL are re-assigned here and feed the
            // adaptors below (the Curve(aNewEdge, ...) range from above is
            // overwritten on every iteration).
            let a_c2d_new = match brep_tool_curve_on_surface(the_edge, &a_cur_face) {
                Some((a_c, a_pcurve_f, a_pcurve_l)) => {
                    a_f = a_pcurve_f;
                    a_l = a_pcurve_l;
                    a_c
                }
                // OCCT reaches this only through the BRep_Tool::CurveOnPlane
                // fallback; the rcad re-host carries no fallback (arch. diff.
                // #1 of loc_ope_wires_on_shape_b), so the iteration is skipped.
                None => continue,
            };
            // OCCT L292: const occ::handle<Adaptor2d_Curve2d> aHCurve2d =
            // new Geom2dAdaptor_Curve(aC2dNew, aF, aL);
            let a_h_curve2d = Geom2dAdaptorCurve::new(a_c2d_new, a_f, a_l);
            // OCCT L293-294: const occ::handle<Adaptor3d_Surface> aHSurface =
            // new GeomAdaptor_Surface(myFaceInfo.Find(aCurFace).myOffsetS);
            let a_h_surface = GeomAdaptorSurface::new(
                a_cur_nfd
                    .my_offset_s
                    .clone()
                    .expect("BRepOffset_SimpleOffset: the face carries an offset surface"),
            );
            // OCCT L295-296: const occ::handle<Adaptor3d_CurveOnSurface>
            // aCurveOnSurf = new Adaptor3d_CurveOnSurface(aHCurve2d, aHSurface);
            let a_curve_on_surf = Adaptor3dCurveOnSurface::new(a_h_curve2d, a_h_surface);

            // Extract 3d-curve (it is not null).
            // OCCT L299: const occ::handle<Adaptor3d_Curve> aCurve3d =
            // new GeomAdaptor_Curve(aNED.myOffsetC, aF, aL);
            let a_curve3d = match &a_ned.my_offset_c {
                Some(a_offset_c) => GeomAdaptorCurve::new(a_offset_c.clone(), a_f, a_l),
                // The loop condition guarantees a non-null myOffsetC.
                None => break,
            };

            // It is necessary to compute maximal deviation (tolerance).
            // OCCT L302-308.
            let mut a_validate_edge = BRepLibValidateEdge::new(a_curve3d, a_curve_on_surf, true);
            a_validate_edge.process();
            if a_validate_edge.is_done() {
                let a_max_tol1 = a_validate_edge.get_max_distance();
                an_edge_tol = an_edge_tol.max(a_max_tol1);
            }
        }
        // OCCT L310: aNED.myTol = std::max(BRep_Tool::Tolerance(aNewEdge),
        // anEdgeTol);
        a_ned.my_tol = brep_tool_tolerance(&a_new_edge).max(an_edge_tol);

        // Save computed 3d curve in map.
        // OCCT L313: myEdgeInfo.Bind(theEdge, aNED);
        self.my_edge_info.insert(shape_key(the_edge), a_ned);
    }

    /// OCCT BRepOffset_SimpleOffset::FillVertexData (cxx L318-427) — Method to
    /// fill new vertex data for single vertex.
    fn fill_vertex_data(
        &mut self,
        the_vertex: &Shape,
        the_vertex_edge_map: &indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)>,
        the_idx: usize,
    ) {
        // Algorithm:
        // Find adjacent edges for the given vertex.
        // Find corresponding end on the each adjacent edge.
        // Get offset points for founded end.
        // Set result vertex position as barycenter of founded points.

        // OCCT L331: gp_Pnt aCurrPnt = BRep_Tool::Pnt(theVertex);
        let a_curr_pnt = brep_tool_pnt(the_vertex).unwrap_or(DVec3::ZERO);

        // OCCT L333: const NCollection_List<TopoDS_Shape>& aEdgesList =
        // theVertexEdgeMap(theIdx);
        let a_edges_list = indexed_data_map_value(the_vertex_edge_map, the_idx);

        if a_edges_list.is_empty() {
            return; // Free verices are skipped.
        }

        // Array to store offset points.
        // OCCT L341: NCollection_DynamicArray<gp_Pnt> anOffsetPointVec;
        let mut an_offset_point_vec: Vec<DVec3> = Vec::new();

        // OCCT L343: double aMaxEdgeTol = 0.0;
        let mut a_max_edge_tol = 0.0f64;

        // Iterate over adjacent edges.
        // OCCT L346-347: NCollection_List<...>::Iterator anIterEdges(aEdgesList);
        // for (; anIterEdges.More(); anIterEdges.Next())
        for a_curr_edge in a_edges_list.iter() {
            // OCCT L351-354: if (!myEdgeInfo.IsBound(aCurrEdge)) continue;
            if !self.my_edge_info.contains_key(&shape_key(a_curr_edge)) {
                continue; // Skip shared edges with wrong orientation.
            }

            // Find the closest bound.
            // OCCT L357-358: double aF, aL; occ::handle<Geom_Curve> aC3d =
            // BRep_Tool::Curve(aCurrEdge, aF, aL);
            // Protection from degenerated edges.
            // OCCT L361-364: if (aC3d.IsNull()) continue;
            let (a_c3d, a_f, a_l) = match brep_tool_curve(a_curr_edge) {
                Some((a_c3d, a_f, a_l)) => (a_c3d, a_f, a_l),
                None => continue,
            };

            // OCCT L366-367: const gp_Pnt aPntF = aC3d->Value(aF); aPntL =
            // aC3d->Value(aL);
            let a_pnt_f = CurveEval::point_at(&a_c3d, a_f);
            let a_pnt_l = CurveEval::point_at(&a_c3d, a_l);

            // OCCT L369-370.
            let a_sq_dist_f = a_pnt_f.distance_squared(a_curr_pnt);
            let a_sq_dist_l = a_pnt_l.distance_squared(a_curr_pnt);

            // OCCT L372-378: double aMinParam = aF, aMaxParam = aL; if
            // (aSqDistL < aSqDistF) { aMinParam = aL; aMaxParam = aF; }
            let (a_min_param, a_max_param) = if a_sq_dist_l < a_sq_dist_f {
                // Square distance to last point is closer.
                (a_l, a_f)
            } else {
                (a_f, a_l)
            };

            // Compute point on offset edge.
            // OCCT L381-384.
            let a_ned = self
                .my_edge_info
                .get(&shape_key(a_curr_edge))
                .expect("bound above");
            let an_offset_curve = a_ned.my_offset_c.clone();
            let an_offset_point = match &an_offset_curve {
                Some(a_offset_c) => CurveEval::point_at(a_offset_c, a_min_param),
                // OCCT dereferences the handle unconditionally; the rcad null
                // handle path skips the edge (architecture difference #2.2).
                None => continue,
            };
            an_offset_point_vec.push(an_offset_point);

            // Handle situation when edge is closed.
            // OCCT L386-393.
            let (a_v1, a_v2) = top_exp_vertices_shape(a_curr_edge);
            if a_v1.is_same(&a_v2) {
                let an_offset_point_last = match &an_offset_curve {
                    Some(a_offset_c) => CurveEval::point_at(a_offset_c, a_max_param),
                    None => continue,
                };
                an_offset_point_vec.push(an_offset_point_last);
            }

            // OCCT L395: aMaxEdgeTol = std::max(aMaxEdgeTol, aNED.myTol);
            a_max_edge_tol = a_max_edge_tol.max(a_ned.my_tol);
        }

        // NCollection_DynamicArray starts from 0 by default.
        // It's better to use lower() and upper() in this case instead of direct
        // indexes range.
        // OCCT L400-405.
        let mut a_center = DVec3::ZERO;
        for an_offset_point in &an_offset_point_vec {
            a_center += *an_offset_point;
        }
        a_center /= an_offset_point_vec.len() as f64;

        // Compute max distance.
        // OCCT L408-416.
        let mut a_sq_max_dist = 0.0f64;
        for an_offset_point in &an_offset_point_vec {
            let a_sq_dist = a_center.distance_squared(*an_offset_point);
            if a_sq_dist > a_sq_max_dist {
                a_sq_max_dist = a_sq_dist;
            }
        }

        // OCCT L418: const double aResTol = std::max(aMaxEdgeTol,
        // std::sqrt(aSqMaxDist));
        let a_res_tol = a_max_edge_tol.max(a_sq_max_dist.sqrt());

        // OCCT L420: const double aMultCoeff = 1.001; // Avoid tolernace problems.
        let a_mult_coeff = 1.001;
        let a_nvd = NewVertexData {
            my_p: a_center,
            my_tol: a_res_tol * a_mult_coeff,
        };

        // Save computed vertex info.
        // OCCT L426: myVertexInfo.Bind(theVertex, aNVD);
        self.my_vertex_info.insert(shape_key(the_vertex), a_nvd);
    }
}

/// OCCT NCollection_IndexedDataMap::operator()(theIdx) — the 1-based value
/// read (the OCCT map indexes run 1..Length()).
fn indexed_data_map_value(
    the_map: &indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)>,
    the_idx: usize,
) -> &[Shape] {
    the_map
        .get_index(the_idx - 1)
        .expect("indexed data map index")
        .1
         .1
        .as_slice()
}

/// OCCT BRepOffset_SimpleOffset -> BRepTools_Modification (hxx L44:
/// `class BRepOffset_SimpleOffset : public BRepTools_Modification`).  The six
/// OCCT overrides (hxx L68-127) are the BRepTools_Modification methods.
impl BRepToolsModification for BRepOffsetSimpleOffset {
    /// OCCT BRepOffset_SimpleOffset::NewSurface (cxx L51-72).
    fn new_surface(
        &mut self,
        the_f: &Shape,
        the_s: &mut Option<Surface3>,
        the_l: &mut u32,
        the_tol: &mut f64,
        the_rev_wires: &mut bool,
        the_rev_face: &mut bool,
    ) -> bool {
        // OCCT L58-61: if (!myFaceInfo.IsBound(F)) return false;
        let Some(a_nfd) = self.my_face_info.get(&shape_key(the_f)) else {
            return false;
        };

        // OCCT L63-69: const NewFaceData& aNFD = myFaceInfo.Find(F); S =
        // aNFD.myOffsetS; L = aNFD.myL; Tol = aNFD.myTol; RevWires =
        // aNFD.myRevWires; RevFace = aNFD.myRevFace;
        *the_s = a_nfd.my_offset_s.clone();
        *the_l = a_nfd.my_l;
        *the_tol = a_nfd.my_tol;
        *the_rev_wires = a_nfd.my_rev_wires;
        *the_rev_face = a_nfd.my_rev_face;

        // OCCT L71: return true;
        true
    }

    /// OCCT BRepOffset_SimpleOffset::NewCurve (cxx L76-93).
    fn new_curve(
        &mut self,
        the_e: &Shape,
        the_c: &mut Option<Curve3>,
        the_l: &mut u32,
        the_tol: &mut f64,
    ) -> bool {
        // OCCT L81-84: if (!myEdgeInfo.IsBound(E)) return false;
        let Some(a_ned) = self.my_edge_info.get(&shape_key(the_e)) else {
            return false;
        };

        // OCCT L86-90: C = aNED.myOffsetC; L = aNED.myL; Tol = aNED.myTol;
        *the_c = a_ned.my_offset_c.clone();
        *the_l = a_ned.my_l;
        *the_tol = a_ned.my_tol;

        // OCCT L92: return true;
        true
    }

    /// OCCT BRepOffset_SimpleOffset::NewPoint (cxx L97-110).
    fn new_point(&mut self, the_v: &Shape, the_p: &mut DVec3, the_tol: &mut f64) -> bool {
        // OCCT L99-102: if (!myVertexInfo.IsBound(V)) return false;
        let Some(a_nvd) = self.my_vertex_info.get(&shape_key(the_v)) else {
            return false;
        };

        // OCCT L104-107: const NewVertexData& aNVD = myVertexInfo.Find(V); P =
        // aNVD.myP; Tol = aNVD.myTol;
        *the_p = a_nvd.my_p;
        *the_tol = a_nvd.my_tol;

        // OCCT L109: return true;
        true
    }

    /// OCCT BRepOffset_SimpleOffset::NewCurve2d (cxx L114-132).
    fn new_curve2d(
        &mut self,
        the_e: &Shape,
        the_f: &Shape,
        _the_new_e: &mut Shape,
        _the_new_f: &Shape,
        the_c: &mut Option<Curve2d>,
        the_tol: &mut f64,
    ) -> bool {
        // Use original pcurve.
        // OCCT L122-124: double aF, aL; C = BRep_Tool::CurveOnSurface(E, F,
        // aF, aL); Tol = BRep_Tool::Tolerance(E);
        *the_c = brep_tool_curve_on_surface(the_e, the_f).map(|(a_c2d, _a_f, _a_l)| a_c2d);
        *the_tol = brep_tool_tolerance(the_e);

        // OCCT L126-129: if (myEdgeInfo.IsBound(E)) Tol =
        // myEdgeInfo.Find(E).myTol;
        if let Some(a_ned) = self.my_edge_info.get(&shape_key(the_e)) {
            *the_tol = a_ned.my_tol;
        }

        // OCCT L131: return true;
        true
    }

    /// OCCT BRepOffset_SimpleOffset::NewParameter (cxx L136-151).
    fn new_parameter(
        &mut self,
        the_v: &Shape,
        the_e: &Shape,
        the_p: &mut f64,
        the_tol: &mut f64,
    ) -> bool {
        // Use original parameter.
        // OCCT L142-143: P = BRep_Tool::Parameter(V, E); Tol =
        // BRep_Tool::Tolerance(V);
        *the_p = brep_tool_parameter(the_v, the_e);
        *the_tol = brep_tool_tolerance(the_v);

        // OCCT L145-148: if (myVertexInfo.IsBound(V)) Tol =
        // myVertexInfo.Find(V).myTol;
        if let Some(a_nvd) = self.my_vertex_info.get(&shape_key(the_v)) {
            *the_tol = a_nvd.my_tol;
        }

        // OCCT L150: return true;
        true
    }

    /// OCCT BRepOffset_SimpleOffset::Continuity (cxx L155-164).
    fn continuity(
        &mut self,
        the_e: &Shape,
        the_f1: &Shape,
        the_f2: &Shape,
        _the_new_e: &Shape,
        _the_new_f1: &Shape,
        _the_new_f2: &Shape,
    ) -> GeomAbsShape {
        // Compute result using original continuity.
        // OCCT L163: return BRep_Tool::Continuity(E, F1, F2);
        brep_tool_continuity(the_e, the_f1, the_f2)
    }
}

/// OCCT ShapeAnalysis_FreeBounds (TKShHealing) — the free-bounds explorer of
/// BuildMissingWalls (architecture difference #3; GAP: no rcad translation
/// yet — the GAP panic is the §0.6 annotation; GetClosedWires keeps the OCCT
/// accessor surface).
pub struct ShapeAnalysisFreeBounds {
    my_closed_wires: Shape, // OCCT: myClosedWires
}

impl ShapeAnalysisFreeBounds {
    /// OCCT ShapeAnalysis_FreeBounds::ShapeAnalysis_FreeBounds(theShape)
    /// (defaults: theSewConnected = false, theShared = false, theSetProjPCur
    /// = false).  GAP: the free-bounds computation is not translated.
    pub fn new(_the_shape: &Shape) -> Self {
        panic!("GAP: ShapeAnalysis_FreeBounds (TKShHealing not translated)");
    }

    /// OCCT ShapeAnalysis_FreeBounds::GetClosedWires().
    pub fn get_closed_wires(&self) -> Shape {
        self.my_closed_wires.clone()
    }
}

/// OCCT BRepTools_Quilt (TKBRep/BRepTools/BRepTools_Quilt.hxx / .cxx) — the
/// real body lives in crate::topalgo::brep_tools_quilt; the sewing of
/// BuildMissingWalls keeps the OCCT import path through the re-export.
pub use crate::topalgo::brep_tools_quilt::BRepToolsQuilt;

/// OCCT ShapeFix_Edge (TKShHealing) — the FixSameParameter of
/// BuildMissingWalls (architecture difference #3; GAP: no rcad translation
/// yet; the context argument keeps the OCCT SetContext form).
fn shape_fix_edge_fix_same_parameter(the_context: &mut ShapeBuildReShape, the_e: &Shape) {
    let _ = (the_context, the_e);
    panic!("GAP: ShapeFix_Edge::FixSameParameter (TKShHealing not translated)");
}

/// OCCT GeomFill_Generator (TKGeomAlgo/GeomFill) — the thrusection generator
/// of BuildWallFace (architecture difference #3; GAP: no rcad translation yet
/// — the GAP panics are the §0.6 annotation).
pub struct GeomFillGenerator;

impl GeomFillGenerator {
    /// OCCT GeomFill_Generator::AddCurve(Curve).
    pub fn add_curve(&mut self, the_curve: &TrimmedCurve3) {
        let _ = the_curve;
        panic!("GAP: GeomFill_Generator::AddCurve (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Generator::Perform(Pres3d).
    pub fn perform(&mut self, the_pres3d: f64) {
        let _ = the_pres3d;
        panic!("GAP: GeomFill_Generator::Perform (TKGeomAlgo/GeomFill not translated)");
    }

    /// OCCT GeomFill_Generator::Surface().
    pub fn surface(&self) -> Surface3 {
        panic!("GAP: GeomFill_Generator::Surface (TKGeomAlgo/GeomFill not translated)");
    }
}

/// OCCT BRepLib_MakeFace(W, OnlyPlane) (TKTopAlgo/BRepLib_MakeFace) — the
/// planar face maker of BuildWallFace (architecture difference #3; GAP: the
/// planar-surface fitting is not translated yet; the carrier keeps the OCCT
/// constructor/IsDone/Face surface with IsDone() = false so the caller takes
/// the OCCT failure branch, the loc_ope_wires_on_shape_b.rs #9 precedent).
pub struct BRepLibMakeFace {
    my_face: Shape, // OCCT: myFace
}

impl BRepLibMakeFace {
    /// OCCT BRepLib_MakeFace::BRepLib_MakeFace(const TopoDS_Wire& W, const
    /// bool OnlyPlane).
    pub fn from_wire(_the_w: &Shape, _the_only_plane: bool) -> Self {
        // GAP: the planar face maker (BRepLib_FindSurface vehicle) is not
        // translated; IsDone() = false reproduces the OCCT failure path.
        BRepLibMakeFace {
            my_face: Shape::null(),
        }
    }

    /// OCCT BRepLib_MakeFace::IsDone().
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT BRepLib_MakeFace::Face().
    pub fn face(&self) -> Shape {
        self.my_face.clone()
    }
}

// ---------------------------------------------------------------------------
// BRep_Tool / TopExp re-hosts (module-private; the loc_ope_* precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::MaxTolerance(theShape, theSubShape)
/// (BRep_Tool.cxx L1792-1830) — the max tolerance over the explored
/// sub-shapes (Face/Edge/Vertex branches).
fn brep_tool_max_tolerance(the_shape: &Shape, the_sub_shape: ShapeType) -> f64 {
    let mut a_tol: f64 = 0.0;

    // Explorer Shape-Subshape.
    let an_exp_ss = explorer(the_shape, the_sub_shape, ShapeType::Shape);
    if the_sub_shape == ShapeType::Face
        || the_sub_shape == ShapeType::Edge
        || the_sub_shape == ShapeType::Vertex
    {
        for a_current_sub_shape in &an_exp_ss {
            a_tol = a_tol.max(brep_tool_tolerance(a_current_sub_shape));
        }
    }

    a_tol
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) (BRep_Tool.cxx L1180-1188 -> the
/// representation scan) — the regularity stored on the (E, F1/F2) curve
/// representations; the default is GeomAbs_C0 when nothing matches
/// (hlr/brep/shape_to_hlr.rs precedent, reduced to the identity-location rcad
/// form).
fn brep_tool_continuity(the_e: &Shape, the_f1: &Shape, the_f2: &Shape) -> GeomAbsShape {
    use rcad_kernel::topo::topods::CurveRepresentation;

    // OCCT: S1 = Surface(F1); S2 = Surface(F2) — the local surfaces.
    let (Some(s1), Some(s2)) = (brep_tool_surface(the_f1), brep_tool_surface(the_f2)) else {
        // The OCCT null-surface case has no matching representation.
        return GeomAbsShape::C0;
    };
    let Some(ed) = the_e.as_edge() else {
        return GeomAbsShape::C0;
    };
    for cr in &ed.representations {
        // OCCT: cr->IsRegularity(S1, S2, l1, l2) — the rcad face locations
        // are identity in this pipeline (arch. diff.: compose_pcurve_location
        // reduced to 0).
        if cr.is_regularity_on(&s1, &s2, 0, 0) {
            // OCCT: return cr->Continuity();
            if let CurveRepresentation::CurveOn2Surfaces { continuity, .. } = cr {
                return *continuity;
            }
        }
    }
    // OCCT: return GeomAbs_C0;
    GeomAbsShape::C0
}

/// OCCT BRep_Tool::IsClosed(theShape) (BRep_Tool.cxx L1707-1771) — the
/// shell/wire/edge branches (the default branch returns theShape.Closed(),
/// read here from the TShape flags).
fn brep_tool_is_closed(the_shape: &Shape) -> bool {
    if the_shape.shape_type() == ShapeType::Shell {
        let mut a_map: HashMap<(u64, u32), ()> = HashMap::new();
        let mut has_bound = false;
        let mut oriented = the_shape.clone();
        oriented.orientation = Orientation::Forward;
        for e in explorer(&oriented, ShapeType::Edge, ShapeType::Shape) {
            if brep_tool_degenerated(&e)
                || e.orientation == Orientation::Internal
                || e.orientation == Orientation::External
            {
                continue;
            }
            has_bound = true;
            // OCCT: if (!aMap.Add(E)) { aMap.Remove(E); } — the parity toggle.
            if a_map.remove(&shape_key(&e)).is_none() {
                a_map.insert(shape_key(&e), ());
            }
        }
        has_bound && a_map.is_empty()
    } else if the_shape.shape_type() == ShapeType::Wire {
        let mut a_map: HashMap<(u64, u32), ()> = HashMap::new();
        let mut has_bound = false;
        let mut oriented = the_shape.clone();
        oriented.orientation = Orientation::Forward;
        for v in explorer(&oriented, ShapeType::Vertex, ShapeType::Shape) {
            if v.orientation == Orientation::Internal || v.orientation == Orientation::External {
                continue;
            }
            has_bound = true;
            if a_map.remove(&shape_key(&v)).is_none() {
                a_map.insert(shape_key(&v), ());
            }
        }
        has_bound && a_map.is_empty()
    } else if the_shape.shape_type() == ShapeType::Edge {
        let (a_v_first, a_v_last) = top_exp_vertices_shape(the_shape);
        !a_v_first.is_null() && a_v_first.is_same(&a_v_last)
    } else {
        // OCCT: return theShape.Closed() — the TShape CLOSED flag.
        tshape_is_closed(the_shape.data.as_ref())
    }
}

/// OCCT TopoDS_Shape::Closed() — the TShape CLOSED flag read across the rcad
/// TShape variants (the flag is a member of every OCCT TopoDS_TShape).
fn tshape_is_closed(the_tshape: &TShape) -> bool {
    match the_tshape {
        TShape::Vertex(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Edge(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Wire(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Face(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Shell(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::Solid(d) => d.flags & tshape_flags::CLOSED != 0,
        TShape::CompSolid(_) => false,
        TShape::Compound(_) => false,
    }
}

/// OCCT TopExp::Vertices(E, VFirst, VLast) — the null Shape form of the
/// loc_ope_wires_on_shape re-host.
pub(crate) fn top_exp_vertices_shape(edg: &Shape) -> (Shape, Shape) {
    let (v1, v2) = top_exp_vertices(edg);
    (
        v1.unwrap_or_else(Shape::null),
        v2.unwrap_or_else(Shape::null),
    )
}

/// OCCT BRepAdaptor_Surface(F, false)::D1(U, V, P, D1U, D1V) — the rcad
/// stand-in evaluates the face surface directly (architecture difference #6):
/// the point through SurfaceEval::point_at, the derivatives through the
/// GeomAdaptor_Surface::DN vehicle.
fn brep_adaptor_surface_d1(the_f: &Shape, the_u: f64, the_v: f64) -> (DVec3, DVec3, DVec3) {
    let s = brep_tool_surface(the_f).expect("BRepAdaptor_Surface: the face carries no surface");
    let p = s.point_at(the_u, the_v);
    let d1u = s.dn(the_u, the_v, 1, 0);
    let d1v = s.dn(the_u, the_v, 0, 1);
    (p, d1u, d1v)
}

/// OCCT Geom_Surface::Bounds(U1, U2, V1, V2) — the rcad stand-in through the
/// SurfaceEval::default_domain vehicle.
fn surface_bounds(the_s: &Surface3) -> (f64, f64, f64, f64) {
    let [u1, u2, v1, v2] = the_s.default_domain();
    (u1, u2, v1, v2)
}

/// OCCT Geom_Surface::UIso(U) — GAP: the iso-curve construction has no rcad
/// translation yet (the kernel iso sampling is private to base/convert).
fn surface_u_iso(the_s: &Surface3, the_u: f64) -> Curve3 {
    let _ = (the_s, the_u);
    panic!("GAP: Geom_Surface::UIso (iso-curve construction not translated)");
}

// ---------------------------------------------------------------------------
// OCCT statics (BRepOffset_MakeSimpleOffset.cxx).
// ---------------------------------------------------------------------------

//=============================================================================
// function : tgtfaces
// purpose  : check the angle at the border between two squares.
//           Two shares should have a shared front edge.
//=============================================================================
// OCCT BRepOffset_MakeSimpleOffset.cxx L215-322.
fn tgtfaces(ed: &Shape, f1: &Shape, f2: &Shape, couture: bool, the_res_angle: &mut f64) {
    // Check that pcurves exist on both faces of edge.
    // OCCT L224-233: aCurve = BRep_Tool::CurveOnSurface(Ed, F1/F2, aFirst,
    // aLast); the handles are consumed by the null checks only.
    if brep_tool_curve_on_surface(ed, f1).is_none() {
        return;
    }
    if brep_tool_curve_on_surface(ed, f2).is_none() {
        return;
    }

    let mut e = ed.clone();
    // OCCT L237-249: BRepAdaptor_Surface aBAS1/aBAS2, HS1/HS2 (reduced to the
    // face handles — architecture difference #6).
    let hs1 = f1.clone();
    let hs2 = if couture { hs1.clone() } else { f2.clone() };
    // case when edge lies on the one face

    e.orientation = Orientation::Forward;
    // OCCT L253: BRepAdaptor_Curve2d C2d1(E, F1).
    let c2d1 = brep_tool_curve_on_surface(&e, f1).map(|(c, _, _)| c);
    if couture {
        e.orientation = Orientation::Reversed;
    }
    // OCCT L258: BRepAdaptor_Curve2d C2d2(E, F2).
    let c2d2 = brep_tool_curve_on_surface(&e, f2).map(|(c, _, _)| c);
    let (Some(c2d1), Some(c2d2)) = (c2d1, c2d2) else {
        // The OCCT adaptor cannot fail here (the pcurves were checked above);
        // the rcad Option form is unpacked defensively.
        return;
    };

    let rev1 = f1.orientation == Orientation::Reversed;
    let rev2 = f2.orientation == Orientation::Reversed;
    let (mut f, mut l) = brep_tool_range(&e);
    // OCCT L264: Extrema_LocateExtPC ext; — declared but never used in the
    // OCCT body; kept as the unused binding.
    let _ext: Option<u8> = None;

    let eps = (l - f) / 100.0;
    f += eps; // to avoid calculations on
    l -= eps; // points of pointed squares.

    const NBPNT: i32 = 23;
    for i in 0..=NBPNT {
        // First suppose that this is sameParameter
        let u = f + (l - f) * i as f64 / NBPNT as f64;

        // take derivatives of surfaces at the same u, and compute normals
        let p = c2d1.point_at(u);
        let (_pp1, du1, dv1) = brep_adaptor_surface_d1(&hs1, p.x, p.y);
        let mut d1 = du1.cross(dv1);
        let norm = d1.length();
        if norm > 1.0e-12 {
            d1 /= norm;
        } else {
            continue; // skip degenerated point
        }
        if rev1 {
            d1 = -d1;
        }

        let p = c2d2.point_at(u);
        let (_pp2, du2, dv2) = brep_adaptor_surface_d1(&hs2, p.x, p.y);
        let mut d2 = du2.cross(dv2);
        let norm = d2.length();
        if norm > 1.0e-12 {
            d2 /= norm;
        } else {
            continue; // skip degenerated point
        }
        if rev2 {
            d2 = -d2;
        }

        // Compute angle.
        let a_current_ang = d1.angle_between(d2);

        *the_res_angle = the_res_angle.max(a_current_ang);
    }
}

//=============================================================================
// function : ComputeMaxAngleOnShape
// purpose  : Code the regularities on all edges of the shape, boundary of
//            two faces that do not have it.
//=============================================================================
// OCCT BRepOffset_MakeSimpleOffset.cxx L329-388.
fn compute_max_angle_on_shape(s: &Shape, the_res_angle: &mut f64) {
    // OCCT L331-333: NCollection_IndexedDataMap<TopoDS_Shape,
    // NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> M;
    // TopExp::MapShapesAndAncestors(S, TopAbs_EDGE, TopAbs_FACE, M);
    let mut m: indexmap::IndexMap<ShapeKey, (Shape, Vec<Shape>)> = indexmap::IndexMap::new();
    map_shapes_and_ancestors(s, ShapeType::Edge, ShapeType::Face, &mut m);
    for (_key, (e, ancestors)) in m.iter() {
        let e = e.clone();
        let mut found = false;
        let mut couture = false;
        let mut f1 = Shape::null();
        let mut f2 = Shape::null();
        // OCCT L344-358: for (It.Initialize(M.FindFromIndex(i));
        // It.More() && !found; It.Next()) { ... }
        for it in ancestors {
            if found {
                break;
            }
            if f1.is_null() {
                f1 = it.clone();
            } else if !f1.is_same(it) {
                found = true;
                f2 = it.clone();
            }
        }
        if !found && !f1.is_null() {
            // is it a sewing edge?
            let or_e = e.orientation;
            // OCCT L363-373: for (Ex.Init(F1, TopAbs_EDGE); Ex.More() &&
            // !found; Ex.Next()) { ... }
            for cur_e in explorer(&f1, ShapeType::Edge, ShapeType::Shape) {
                if found {
                    break;
                }
                if e.is_same(&cur_e) && or_e != cur_e.orientation {
                    found = true;
                    couture = true;
                    f2 = f1.clone();
                }
            }
        }
        if found {
            // OCCT L376: if (BRep_Tool::Continuity(E, F1, F2) <= GeomAbs_C0).
            if brep_tool_continuity(&e, &f1, &f2) <= GeomAbsShape::C0 {
                // OCCT L378-384: the Standard_Failure catch is a no-op; the
                // rcad tgtfaces does not panic on the same inputs.
                tgtfaces(&e, &f1, &f2, couture, the_res_angle);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT class (BRepOffset_MakeSimpleOffset.hxx L61-172, cxx L46-703).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_MakeSimpleOffset (BRepOffset_MakeSimpleOffset.hxx L61-172).
pub struct BRepOffsetMakeSimpleOffset {
    // Input data.
    my_input_shape: Shape,           // OCCT: myInputShape (hxx L134)
    my_offset_value: f64,            // OCCT: myOffsetValue (hxx L137)
    my_tolerance: f64,               // OCCT: myTolerance (hxx L140)
    my_is_build_solid: bool,         // OCCT: myIsBuildSolid (hxx L143)

    // Internal data.
    my_max_angle: f64,               // OCCT: myMaxAngle (hxx L148)
    my_error: BRepOffsetSimpleStatus, // OCCT: myError (hxx L151)
    my_is_done: bool,                // OCCT: myIsDone (hxx L154)
    my_map_ve: HashMap<ShapeKey, Shape>, // OCCT: myMapVE (hxx L158)
    my_builder: BRepToolsModifier,   // OCCT: myBuilder (hxx L161)
    my_re_shape: ShapeBuildReShape,  // OCCT: myReShape (hxx L164)
    my_brep: BRep,                   // rcad arena stand-in (arch. diff. #4)

    // Output data.
    my_res_shape: Shape,             // OCCT: myResShape (hxx L169)
}

impl Default for BRepOffsetMakeSimpleOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetMakeSimpleOffset {
    /// OCCT BRepOffset_MakeSimpleOffset::BRepOffset_MakeSimpleOffset()
    /// (cxx L48-57).
    pub fn new() -> Self {
        BRepOffsetMakeSimpleOffset {
            my_input_shape: Shape::null(),
            my_offset_value: 0.0,
            my_tolerance: rcad_kernel::core::precision::CONFUSION,
            my_is_build_solid: false,
            my_max_angle: 0.0,
            my_error: BRepOffsetSimpleStatus::Ok,
            my_is_done: false,
            my_map_ve: HashMap::new(),
            my_builder: BRepToolsModifier::new(false),
            my_re_shape: ShapeBuildReShape::new(),
            my_brep: BRep::new(),
            my_res_shape: Shape::null(),
        }
    }

    /// OCCT BRepOffset_MakeSimpleOffset::BRepOffset_MakeSimpleOffset(
    /// theInputShape, theOffsetValue) (cxx L61-72) — Rust has no overloading.
    pub fn with_input(the_input_shape: &Shape, the_offset_value: f64) -> Self {
        let mut res = Self::new();
        res.my_input_shape = the_input_shape.clone();
        res.my_offset_value = the_offset_value;
        res
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Initialize (cxx L76-82) — Initialise
    /// shape for modifications.
    pub fn initialize(&mut self, the_input_shape: &Shape, the_offset_value: f64) {
        self.my_input_shape = the_input_shape.clone();
        self.my_offset_value = the_offset_value;
        self.clear();
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetErrorMessage (cxx L86-117) —
    /// Gets error message.
    pub fn get_error_message(&self) -> &'static str {
        if self.my_error == BRepOffsetSimpleStatus::NullInputShape {
            return "Null input shape";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorOffsetComputation {
            return "Error during offset construction";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorWallFaceComputation {
            return "Error during building wall face";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorInvalidNbShells {
            return "Result contains two or more shells";
        } else if self.my_error == BRepOffsetSimpleStatus::ErrorNonClosedShell {
            return "Result shell is not closed";
        }
        ""
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetError() (hxx L81).
    pub fn get_error(&self) -> BRepOffsetSimpleStatus {
        self.my_error
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetBuildSolidFlag() (hxx L85).
    pub fn get_build_solid_flag(&self) -> bool {
        self.my_is_build_solid
    }

    /// OCCT BRepOffset_MakeSimpleOffset::SetBuildSolidFlag(theBuildFlag)
    /// (hxx L88).
    pub fn set_build_solid_flag(&mut self, the_build_flag: bool) {
        self.my_is_build_solid = the_build_flag;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetOffsetValue() (hxx L91).
    pub fn get_offset_value(&self) -> f64 {
        self.my_offset_value
    }

    /// OCCT BRepOffset_MakeSimpleOffset::SetOffsetValue(theOffsetValue)
    /// (hxx L94).
    pub fn set_offset_value(&mut self, the_offset_value: f64) {
        self.my_offset_value = the_offset_value;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetTolerance() (hxx L97).
    pub fn get_tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT BRepOffset_MakeSimpleOffset::SetTolerance(theValue) (hxx L100).
    pub fn set_tolerance(&mut self, the_value: f64) {
        self.my_tolerance = the_value;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::IsDone() (hxx L103).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetResultShape() (hxx L106).
    pub fn get_result_shape(&self) -> Shape {
        self.my_res_shape.clone()
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Clear (cxx L121-128) — Clears
    /// previous result.
    fn clear(&mut self) {
        self.my_is_done = false;
        self.my_error = BRepOffsetSimpleStatus::Ok;
        self.my_max_angle = 0.0;
        self.my_map_ve.clear();
        self.my_re_shape.clear(); // Clear possible stored modifications.
    }

    /// OCCT BRepOffset_MakeSimpleOffset::GetSafeOffset (cxx L132-151) —
    /// Computes max safe offset value for the given tolerance.
    pub fn get_safe_offset(&mut self, the_expected_toler: f64) -> f64 {
        if self.my_input_shape.is_null() {
            return 0.0; // Input shape is null.
        }

        // Compute max angle in faces junctions.
        if self.my_max_angle == 0.0 {
            // Non-initialized.
            self.compute_max_angle();
        }

        // OCCT L146: aMaxTol = BRep_Tool::MaxTolerance(myInputShape,
        // TopAbs_VERTEX).
        let a_max_tol = brep_tool_max_tolerance(&self.my_input_shape, ShapeType::Vertex);

        // OCCT L148-149: std::max((theExpectedToler - aMaxTol) /
        // (2.0 * myMaxAngle), 0.0) — Minimal distance can't be lower than 0.0.
        let an_exp_offset = ((the_expected_toler - a_max_tol) / (2.0 * self.my_max_angle)).max(0.0);
        an_exp_offset
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Perform (cxx L155-208) — Computes
    /// offset shape.
    pub fn perform(&mut self) {
        // Clear result of previous computations.
        self.clear();

        // Check shape existence.
        if self.my_input_shape.is_null() {
            self.my_error = BRepOffsetSimpleStatus::NullInputShape;
            return;
        }

        if self.my_max_angle == 0.0 {
            // Non-initialized.
            self.compute_max_angle();
        }

        // OCCT L172: myBuilder.Init(myInputShape).
        self.my_builder.init(&self.my_input_shape);
        // OCCT L173-174: occ::handle<BRepOffset_SimpleOffset> aMapper = new
        // BRepOffset_SimpleOffset(myInputShape, myOffsetValue, myTolerance).
        let mut a_mapper = BRepOffsetSimpleOffset::new(
            &self.my_input_shape,
            self.my_offset_value,
            self.my_tolerance,
        );
        // OCCT L175: myBuilder.Perform(aMapper).
        self.my_builder.perform(&mut a_mapper);

        if !self.my_builder.is_done() {
            self.my_error = BRepOffsetSimpleStatus::ErrorOffsetComputation;
            return;
        }

        self.my_res_shape = self.my_builder.modified_shape(&self.my_input_shape);

        // Fix degeneracy. Degenerated edge should be mapped to the degenerated.
        // OCCT L186-199: BRep_Builder aBB; — the rcad carrier of the
        // BRep_Builder::Degenerated call is
        // brep_algo::tool::builder_set_degenerated (the rebuilt shapes of
        // BRepTools_Modifier are pool-free, so the pool-addressed BRepBuilder
        // API cannot reach them — arch. diff. #4).
        let an_exp_se = explorer(&self.my_input_shape, ShapeType::Edge, ShapeType::Shape);
        for a_curr_edge in &an_exp_se {
            if !brep_tool_degenerated(a_curr_edge) {
                continue;
            }

            // OCCT L197: const TopoDS_Edge& anEdge = TopoDS::Edge(myBuilder.ModifiedShape(aCurrEdge));
            let mut an_edge = self.my_builder.modified_shape(a_curr_edge);
            // OCCT L198: aBB.Degenerated(anEdge, true).
            crate::brep_algo::tool::builder_set_degenerated(&mut an_edge, true);
        }

        // Restore walls for solid.
        if self.my_is_build_solid && !self.build_missing_walls() {
            return;
        }

        self.my_is_done = true;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::ComputeMaxAngle (cxx L394-397) —
    /// Computes max angle in faces junction.
    fn compute_max_angle(&mut self) {
        let mut res_angle = self.my_max_angle;
        compute_max_angle_on_shape(&self.my_input_shape, &mut res_angle);
        self.my_max_angle = res_angle;
    }

    /// OCCT BRepOffset_MakeSimpleOffset::BuildMissingWalls (cxx L403-509) —
    /// Builds walls to the result solid.
    fn build_missing_walls(&mut self) -> bool {
        // Internal list of new faces.
        let mut a_bb = BRepBuilder::new();
        let a_new_faces = a_bb.make_compound(&mut self.my_brep, Vec::new());

        // Compute outer bounds of original shape.
        // OCCT L411-412: ShapeAnalysis_FreeBounds aFB(myInputShape);
        // GetClosedWires.
        let a_fb = ShapeAnalysisFreeBounds::new(&self.my_input_shape);
        let a_free_wires = a_fb.get_closed_wires();

        // Build linear faces on each edge and its image.
        let an_exp_cw = explorer(&a_free_wires, ShapeType::Wire, ShapeType::Shape);
        for a_cur_wire in &an_exp_cw {
            // Iterate over outer edges in outer wires.
            let an_exp_we = explorer(a_cur_wire, ShapeType::Edge, ShapeType::Shape);
            for a_cur_edge in &an_exp_we {
                let a_new_face = self.build_wall_face(a_cur_edge);

                if a_new_face.is_null() {
                    self.my_error = BRepOffsetSimpleStatus::ErrorWallFaceComputation;
                    return false;
                }

                // OCCT L434: aBB.Add(aNewFaces, aNewFace).
                a_bb.add_to_compound(&mut self.my_brep, a_new_faces.clone(), a_new_face);
            }
        }

        // Update edges from wall faces.
        // OCCT L439-447: ShapeFix_Edge aSFE; aSFE.SetContext(myReShape);
        // aSFE.FixSameParameter(aCurrEdge).
        for a_curr_edge in explorer(&a_new_faces, ShapeType::Edge, ShapeType::Shape) {
            // Fix same parameter and same range flags.
            shape_fix_edge_fix_same_parameter(&mut self.my_re_shape, &a_curr_edge);
        }

        // Update result to be compound.
        // OCCT L450-451: TopoDS_Compound aResCompound; aBB.MakeCompound.
        let a_res_compound = a_bb.make_compound(&mut self.my_brep, Vec::new());

        // Add old faces the result.
        for a_f in explorer(&self.my_input_shape, ShapeType::Face, ShapeType::Shape) {
            a_bb.add_to_compound(&mut self.my_brep, a_res_compound.clone(), a_f);
        }

        // Add new faces the result.
        for a_f in explorer(&self.my_res_shape, ShapeType::Face, ShapeType::Shape) {
            a_bb.add_to_compound(&mut self.my_brep, a_res_compound.clone(), a_f);
        }

        // Add wall faces to the result.
        for a_f in explorer(&a_new_faces, ShapeType::Face, ShapeType::Shape) {
            a_bb.add_to_compound(&mut self.my_brep, a_res_compound.clone(), a_f);
        }

        // Apply stored modifications.
        // OCCT L476: aResCompound = TopoDS::Compound(myReShape->Apply(...)).
        let a_res_compound = self
            .my_re_shape
            .apply(&mut self.my_brep, &a_res_compound, ShapeType::Shape);

        // Create result shell.
        // OCCT L479-481: BRepTools_Quilt aQuilt; aQuilt.Add(aResCompound);
        // aShells = aQuilt.Shells().
        let mut a_quilt = BRepToolsQuilt::new();
        a_quilt.add(&a_res_compound);
        let a_shells = a_quilt.shells();

        let mut a_res_shell = Shape::null();
        for a_shell in explorer(&a_shells, ShapeType::Shell, ShapeType::Shape) {
            if !a_res_shell.is_null() {
                // Shell is not null -> explorer contains two or more shells.
                self.my_error = BRepOffsetSimpleStatus::ErrorInvalidNbShells;
                return false;
            }
            a_res_shell = a_shell;
        }

        // OCCT L496: if (!BRep_Tool::IsClosed(aResShell)).
        if !brep_tool_is_closed(&a_res_shell) {
            self.my_error = BRepOffsetSimpleStatus::ErrorNonClosedShell;
            return false;
        }

        // Create result solid.
        // OCCT L503-506: aBB.MakeSolid(aResSolid); aBB.Add(aResSolid,
        // aResShell); myResShape = aResSolid.
        let a_res_solid = a_bb.make_solid(&mut self.my_brep, vec![a_res_shell]);
        self.my_res_shape = a_res_solid;

        true
    }

    /// OCCT BRepOffset_MakeSimpleOffset::BuildWallFace (cxx L513-667) —
    /// Builds face on specified wall.
    fn build_wall_face(&mut self, the_orig_edge: &Shape) -> Shape {
        let a_res_face = Shape::null();

        // Get offset edge. offset edge is reversed to create correct wire.
        // OCCT L518: aNewEdge = TopoDS::Edge(myBuilder.ModifiedShape(...)).
        let mut a_new_edge = self.my_builder.modified_shape(the_orig_edge);
        a_new_edge.orientation = Orientation::Reversed;

        // OCCT L522: TopExp::Vertices(aNewEdge, aNewV1, aNewV2).
        let (a_new_v1, a_new_v2) = top_exp_vertices_shape(&a_new_edge);

        // Wire contour is:
        // theOrigEdge (forcible forward) -> wall1 -> aNewEdge (forcible reversed) -> wall2
        // Firstly it is necessary to create copy of original shape with forward direction.
        // This simplifies walls creation.
        let mut an_orig_copy = the_orig_edge.clone();
        an_orig_copy.orientation = Orientation::Forward;
        // OCCT L531: TopExp::Vertices(anOrigCopy, aV1, aV2).
        let (a_v1, a_v2) = top_exp_vertices_shape(&an_orig_copy);

        // To simplify work with map.
        // OCCT L534-535: TopoDS::Vertex(aV1.Oriented(TopAbs_FORWARD)).
        let a_forward_v1 = oriented_vertex(&a_v1, Orientation::Forward);
        let a_forward_v2 = oriented_vertex(&a_v2, Orientation::Forward);

        // Check existence of edges in stored map: Edge1
        let a_wall1 = if self.my_map_ve.contains_key(&shape_key(&a_forward_v2)) {
            // Edge exists - get it from map.
            self.my_map_ve[&shape_key(&a_forward_v2)].clone()
        } else {
            // Edge does not exist - create it and add to the map.
            // OCCT L547-552: BRepLib_MakeEdge aME1(...); if (!aME1.IsDone())
            // return aResFace; — the rcad BRepBuilder::add_edge cannot fail
            // (architecture difference #5).
            let mut a_bb = BRepBuilder::new();
            let a_me1 = a_bb.add_edge(
                &mut self.my_brep,
                None,
                oriented_vertex(&a_v2, Orientation::Forward),
                oriented_vertex(&a_new_v2, Orientation::Reversed),
                [0.0, 0.0],
            );
            let a_wall1 = a_me1;

            self.my_map_ve
                .insert(shape_key(&a_forward_v2), a_wall1.clone());
            a_wall1
        };

        // Check existence of edges in stored map: Edge2
        let a_wall2 = if self.my_map_ve.contains_key(&shape_key(&a_forward_v1)) {
            // Edge exists - get it from map.
            // OCCT L563: TopoDS::Edge(myMapVE(aForwardV1).Oriented(
            // TopAbs_REVERSED)).
            let mut w = self.my_map_ve[&shape_key(&a_forward_v1)].clone();
            w.orientation = Orientation::Reversed;
            w
        } else {
            // Edge does not exist - create it and add to the map.
            // OCCT L568-580: BRepLib_MakeEdge aME2(...) — the IsDone guard has
            // no rcad counterpart (architecture difference #5).
            let mut a_bb = BRepBuilder::new();
            let a_me2 = a_bb.add_edge(
                &mut self.my_brep,
                None,
                oriented_vertex(&a_v1, Orientation::Forward),
                oriented_vertex(&a_new_v1, Orientation::Reversed),
                [0.0, 0.0],
            );
            let mut a_wall2 = a_me2;

            self.my_map_ve
                .insert(shape_key(&a_forward_v1), a_wall2.clone());

            // Orient it in reversed direction.
            a_wall2.orientation = Orientation::Reversed;
            a_wall2
        };

        let mut a_bb = BRepBuilder::new();

        // OCCT L584-589: aBB.MakeWire(aWire); aBB.Add(aWire, ...) x4.
        let a_wire = a_bb.make_wire(&mut self.my_brep);
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), an_orig_copy.clone());
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), a_wall1.clone());
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), a_new_edge.clone());
        a_bb.add_to_wire(&mut self.my_brep, a_wire.clone(), a_wall2.clone());

        // Build 3d curves on wire
        // OCCT L592: BRepLib::BuildCurves3d(aWire) (BRepLib.cxx L460-464).
        crate::topalgo::brep_lib::brep_lib::BRepLib::build_curves3d(&mut self.my_brep, &a_wire);

        // Try to build using simple planar approach.
        // OCCT L595-607: the face maker is wrapped by try/catch since it
        // generates exceptions sometimes; the rcad BRepLibMakeFace carrier
        // keeps IsDone() = false (architecture difference #3).
        let mut a_f = Shape::null();
        {
            // Call of face maker is wrapped by try/catch since it generates exceptions sometimes.
            let a_fm = BRepLibMakeFace::from_wire(&a_wire, true);
            if a_fm.is_done() {
                a_f = a_fm.face();
            }
        }

        if a_f.is_null() {
            // Exception in face maker or result is not computed.
            // Build using thrusections.
            // OCCT L612: bool ToReverse = false.
            let to_reverse = false;
            // OCCT L614-615: EdgeCurve = BRep_Tool::Curve(theOrigEdge, fpar,
            // lpar); TrEdgeCurve = new Geom_TrimmedCurve(EdgeCurve, fpar,
            // lpar).
            let (edge_curve, fpar, lpar) = match edge_curve_of(the_orig_edge) {
                Some(v) => v,
                None => return a_res_face,
            };
            let tr_edge_curve = TrimmedCurve3::new(edge_curve.clone(), fpar, lpar);
            // OCCT L616-618: OffsetCurve = BRep_Tool::Curve(aNewEdge, fparOE,
            // lparOE); TrOffsetCurve = new Geom_TrimmedCurve(...).
            let (offset_curve, fpar_oe, lpar_oe) = match edge_curve_of(&a_new_edge) {
                Some(v) => v,
                None => return a_res_face,
            };
            let tr_offset_curve = TrimmedCurve3::new(offset_curve.clone(), fpar_oe, lpar_oe);

            // OCCT L620-624: GeomFill_Generator ThrusecGenerator;
            // AddCurve x2; Perform(Precision::PConfusion()); theSurf =
            // Surface().
            let mut thrusec_generator = GeomFillGenerator;
            thrusec_generator.add_curve(&tr_edge_curve);
            thrusec_generator.add_curve(&tr_offset_curve);
            thrusec_generator.perform(rcad_kernel::core::precision::PCONFUSION);
            let the_surf = thrusec_generator.surface();
            // OCCT L626-627: theSurf->Bounds(Uf, Ul, Vf, Vl).
            let (uf, ul, vf, vl) = surface_bounds(&the_surf);
            // OCCT L628: TopLoc_Location Loc; — the rcad location index
            // (0 = identity).
            let loc: u32 = 0;
            // OCCT L629-633: Geom2d_Line pcurves bound at (0, Vf)/(0, Vl) in
            // the X direction.  GAP: the BRep_Builder::UpdateEdge(E, C2d, S,
            // L, Tol) surface-keyed pcurve storage has no rcad carrier — the
            // calls are kept as no-op re-hosts (annotated).
            let edge_line2d = line2d_x(0.0, vf);
            update_edge_pcurve_on_surface(&the_orig_edge, &edge_line2d, &the_surf, loc);
            let oe_line2d = line2d_x(0.0, vl);
            update_edge_pcurve_on_surface(&a_new_edge, &oe_line2d, &the_surf, loc);
            // OCCT L634-637.
            let u_on_v1 = if to_reverse { ul } else { uf };
            let u_on_v2 = if to_reverse { uf } else { ul };
            let a_line2d = line2d_y(u_on_v2, 0.0);
            let a_line2d2 = line2d_y(u_on_v1, 0.0);
            if a_wall1.is_same(&a_wall2) {
                // OCCT L640: aBB.UpdateEdge(aWall1, aLine2d, aLine2d2,
                // theSurf, Loc, Precision::Confusion()).
                update_edge_pcurves_on_surface(&a_wall1, &a_line2d, &a_line2d2, &the_surf, loc);
                // OCCT L641-642: BSplC34 = theSurf->UIso(Uf);
                // aBB.UpdateEdge(aWall1, BSplC34, Precision::Confusion()).
                let bspl_c34 = surface_u_iso(&the_surf, uf);
                update_edge_curve3d_gap(&a_wall1, &bspl_c34);
                // OCCT L643: aBB.Range(aWall1, Vf, Vl).
                a_bb.set_edge_range(&mut self.my_brep, a_wall1.clone(), vf, vl);
            } else {
                // OCCT L647-650: aBB.SameParameter/SameRange(aWall1/aWall2,
                // false).
                a_bb.set_edge_same_parameter(&mut self.my_brep, a_wall1.clone(), false);
                a_bb.set_edge_same_range(&mut self.my_brep, a_wall1.clone(), false);
                a_bb.set_edge_same_parameter(&mut self.my_brep, a_wall2.clone(), false);
                a_bb.set_edge_same_range(&mut self.my_brep, a_wall2.clone(), false);
                // OCCT L651-654.
                update_edge_pcurve_on_surface(&a_wall1, &a_line2d, &the_surf, loc);
                set_edge_range_on_surface_gap(&a_wall1, &the_surf, loc, vf, vl);
                update_edge_pcurve_on_surface(&a_wall2, &a_line2d2, &the_surf, loc);
                set_edge_range_on_surface_gap(&a_wall2, &the_surf, loc, vf, vl);
                // OCCT L655-660: BSplC3 = theSurf->UIso(UonV2);
                // aBB.UpdateEdge(aWall1, BSplC3, ...); aBB.Range(aWall1, Vf,
                // Vl, true); BSplC4 = theSurf->UIso(UonV1);
                // aBB.UpdateEdge(aWall2, BSplC4, ...); aBB.Range(aWall2, Vf,
                // Vl, true).
                let bspl_c3 = surface_u_iso(&the_surf, u_on_v2);
                update_edge_curve3d_gap(&a_wall1, &bspl_c3);
                set_edge_range3d_gap(&a_wall1, vf, vl); // only for 3d curve
                let bspl_c4 = surface_u_iso(&the_surf, u_on_v1);
                update_edge_curve3d_gap(&a_wall2, &bspl_c4);
                set_edge_range3d_gap(&a_wall2, vf, vl); // only for 3d curve
            }

            // OCCT L663: aF = BRepLib_MakeFace(theSurf, aWire) — the
            // rcad BRepBuilder::make_face is the same storage-only stand-in
            // for the BRepLib_MakeFace(S, W) constructor.
            a_f = a_bb.make_face(&mut self.my_brep, Some(the_surf), a_wire);
        }

        a_f
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Generated (cxx L671-686) — Returns
    /// result shape for the given one (if exists).  The OCCT const method
    /// takes &mut self (architecture difference #7: the rcad ReShape Apply).
    pub fn generated(&mut self, the_shape: &Shape) -> Shape {
        // Shape generated by modification.
        let mut a_res = self.my_builder.modified_shape(the_shape);

        if a_res.is_null() {
            return a_res;
        }

        // Shape modifications obtained in scope of shape healing.
        // OCCT L683: aRes = myReShape->Apply(aRes) (the until/until shape
        // argument defaults to TopAbs_SHAPE).
        a_res = self
            .my_re_shape
            .apply(&mut self.my_brep, &a_res, ShapeType::Shape);

        a_res
    }

    /// OCCT BRepOffset_MakeSimpleOffset::Modified (cxx L690-703) — Returns
    /// modified shape for the given one (if exists).  The OCCT const method
    /// takes &mut self (architecture difference #7: the rcad ReShape Status).
    pub fn modified(&mut self, the_shape: &Shape) -> Shape {
        let an_empty_shape = Shape::null();

        // Get modification status and new shape.
        // OCCT L695: int aModStatus = myReShape->Status(theShape, aRes).
        let (a_mod_status, a_res) = self
            .my_re_shape
            .status(&mut self.my_brep, the_shape, false);

        if a_mod_status == 0 {
            return an_empty_shape; // No modifications are applied to the shape or its sub-shapes.
        }

        a_res
    }
}

// ---------------------------------------------------------------------------
// Small local helpers (OCCT expression stand-ins).
// ---------------------------------------------------------------------------

/// OCCT aV.Oriented(TopAbs_FORWARD) — TopoDS::Vertex cast included.
pub(crate) fn oriented_vertex(the_v: &Shape, the_orient: Orientation) -> Shape {
    let mut s = the_v.clone();
    s.orientation = the_orient;
    s
}

/// OCCT BRep_Tool::Curve(E, f, l) — the (curve, fpar, lpar) triple (the
/// loc_ope_wires_on_shape_b re-host in the Option form).
pub(crate) fn edge_curve_of(the_e: &Shape) -> Option<(Curve3, f64, f64)> {
    match the_e.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT new Geom2d_Line(gp_Pnt2d(x, y), gp_Dir2d(gp_Dir2d::D::X)) — the
/// X-direction 2d line stand-in (Geom2d_Line -> Curve2d::Line(Line2d)).
fn line2d_x(x: f64, y: f64) -> Curve2d {
    Curve2d::Line(Line2d::new(DVec2::new(x, y), DVec2::new(1.0, 0.0)))
}

/// OCCT new Geom2d_Line(gp_Pnt2d(x, y), gp_Dir2d(gp_Dir2d::D::Y)) — the
/// Y-direction 2d line stand-in.
fn line2d_y(x: f64, y: f64) -> Curve2d {
    Curve2d::Line(Line2d::new(DVec2::new(x, y), DVec2::new(0.0, 1.0)))
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, S, L, Tol) — GAP: the rcad
/// BRepBuilder carries only face-keyed pcurve storage; the surface-keyed
/// form is a no-op re-host (annotated at the call site).
fn update_edge_pcurve_on_surface(the_e: &Shape, the_c2d: &Curve2d, the_s: &Surface3, the_l: u32) {
    let _ = (the_e, the_c2d, the_s, the_l);
}

/// OCCT BRep_Builder::UpdateEdge(E, C1, C2, S, L, Tol) — the two-pcurve
/// surface-keyed no-op re-host.
fn update_edge_pcurves_on_surface(
    the_e: &Shape,
    the_c1: &Curve2d,
    the_c2: &Curve2d,
    the_s: &Surface3,
    the_l: u32,
) {
    let _ = (the_e, the_c1, the_c2, the_s, the_l);
}

/// OCCT BRep_Builder::UpdateEdge(E, C3d, Tol) — GAP no-op re-host (the rcad
/// 3d-curve representation is index-based, not value-based).
fn update_edge_curve3d_gap(the_e: &Shape, the_c: &Curve3) {
    let _ = (the_e, the_c);
}

/// OCCT BRep_Builder::Range(E, S, L, f, l) — GAP no-op re-host
/// (surface-keyed range storage).
fn set_edge_range_on_surface_gap(
    the_e: &Shape,
    the_s: &Surface3,
    the_l: u32,
    the_f: f64,
    the_last: f64,
) {
    let _ = (the_e, the_s, the_l, the_f, the_last);
}

/// OCCT BRep_Builder::Range(E, f, l, only3d = true) — GAP no-op re-host (the
/// rcad BRepBuilder has no only-3d range form).
fn set_edge_range3d_gap(the_e: &Shape, the_f: f64, the_l: f64) {
    let _ = (the_e, the_f, the_l);
}
