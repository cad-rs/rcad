//! OCCT TopOpeBRepBuild_HBuilder facade for the ChFi3d call surface.
//!
//! User ruling 2026-09-07 (D6): TKBool code (TopOpeBRepBuild / TopOpeBRepDS /
//! TopOpeBRepTool) is NOT translated 1:1.  Wherever TKFillet / TKOffset /
//! TKFeat depend on it, the rcad implementation goes through the aligned
//! TKBO pipeline (`bop::algo` / `bop::ds` / `bop::brep_algo_api`).
//!
//! The HBuilder surface consumed by ChFi3d (Builder.hxx) is the 7-method
//! set Perform / MergeSolid / IsSplit / Splits / Merged / NewEdges /
//! NewFaces.  Stage 1g wires those methods over the rcad TKBO pipeline:
//!
//!   - Perform(HDS) rebuilds the reconstruction tables from the ChFi3d DS
//!     facade (fillet/chfi3d_ds.rs): BuildVertices / BuildEdges / BuildFaces
//!     materialize the DS geometry payloads (side tables + interference
//!     records) into concrete topods shapes driven through the rcad BRep
//!     builders, and the split pass mirrors TopOpeBRepBuild_Builder::
//!     SplitEdge over the facade per-shape interference records.
//!   - MergeSolid(S, TB) runs the TKBO-equivalent solid reconstruction:
//!     BOPAlgo_BuilderSolid over the object solid's faces plus the new
//!     faces (the BuildFaces products, attached through the DS
//!     SolidSurface interferences) — the MakeSolids semantics (faces ->
//!     closed shells -> solids with classification) of OCCT
//!     MergeShapes/SplitSolid/MakeSolids.
//!   - IsSplit / Splits / Merged read the split / merged result tables.
//!   - NewEdges / NewFaces read the reconstruction products of Perform.
//!
//! Result-table data structures mirror OCCT TopOpeBRepBuild_Builder.hxx:
//! `mySplitOUT/IN/ON` and `myMergedOUT/IN/ON` are per-shape maps of
//! TopOpeBRepDS_ListOfShapeOn1State (keyed by TShape + Location, the
//! TopTools_ShapeMapHasher identity, carried as `(ptr_id, location)`), and
//! `myNewEdges` / `myNewFaces` / `myNewVertices` are the DS
//! curve / surface / point indexed lists.

use std::collections::HashMap;

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::Curve2dEval as _;
use rcad_kernel::geom::CurveEval as _;
use rcad_kernel::topo::topods::{self, BRep, BRepBuilder, BRepTool as _, Orientation, Shape};
use rcad_kernel::topods::TShape;

use super::chfi3d_builder_2::TopAbsState;
use super::chfi3d_ds::{TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference, TopOpeBRepDSKind};
use crate::bop::algo::builder_solid::BuilderSolid;

// =========================================================================
// OCCT TopOpeBRepDS_ListOfShapeOn1State — the per-shape value type of the
// Builder split / merged maps (TopOpeBRepBuild_Builder.hxx L100-108):
// the IsSplit / merge flag plus the list of shapes lying on the state.
// =========================================================================
#[derive(Debug, Clone, Default)]
pub struct ListOfShapeOn1State {
    /// OCCT: bool myIsSplit (MarkSplit(S, ToBuild, Bval) sets it).
    pub is_split: bool,
    /// OCCT: NCollection_List<TopoDS_Shape> myListOnState.
    pub list_on_state: Vec<Shape>,
}

/// OCCT NCollection_DataMap<TopoDS_Shape, TopOpeBRepDS_ListOfShapeOn1State,
/// TopTools_ShapeMapHasher>.  TopTools_ShapeMapHasher identity is
/// TShape + Location (orientation ignored) — the rcad key is
/// `(ptr_id, location)`, the same identity the bop `myImages` map uses.
pub type ShapeStateMap = HashMap<(u64, u32), ListOfShapeOn1State>;

/// OCCT TopOpeBRepBuild_HBuilder (TopOpeBRepBuild_Builder handle wrapper).
/// D6 ruling: implemented over the rcad TKBO pipeline; see the module doc.
#[derive(Debug, Clone, Default)]
pub struct TopOpeBRepBuildHBuilder {
    /// OCCT Builder.hxx: mySplitOUT / mySplitIN / mySplitON.
    pub my_split_out: ShapeStateMap,
    pub my_split_in: ShapeStateMap,
    pub my_split_on: ShapeStateMap,
    /// OCCT Builder.hxx: myMergedOUT / myMergedIN / myMergedON.
    pub my_merged_out: ShapeStateMap,
    pub my_merged_in: ShapeStateMap,
    pub my_merged_on: ShapeStateMap,
    /// OCCT Builder.hxx: myNewEdges (IndexedDataMap<int, list>) — the edges
    /// built on DS curve I (1-based DS curve index).
    pub my_new_edges: HashMap<i32, Vec<Shape>>,
    /// OCCT Builder.hxx: myNewFaces (HArray1(0, NbSurfaces)) — the faces
    /// built on DS surface I (1-based DS surface index).
    pub my_new_faces: HashMap<i32, Vec<Shape>>,
    /// OCCT Builder.hxx: myNewVertices (HArray1(0, NbPoints)) — the vertex
    /// built on DS point I (1-based DS point index).
    pub my_new_vertices: HashMap<i32, Shape>,
    /// OCCT Builder.hxx L122: myDataStructure = HDS (stored at Perform).
    pub my_data_structure: Option<TopOpeBRepDSHDataStructure>,
}

// D6 architecture (pool identity, converged): OCCT shapes are pool-free
// handles and rcad TopoDS_Shape resolution is pool-local
// (`tshapes[index]`).  All shapes materialized by Perform / MergeSolid are
// built into the CALLER's pool (the ChFi3d builder's my_brep, passed to
// perform / merge_solid) — the same pool the DS shapes and the result
// compound live in — so every handle the HBuilder hands back resolves in
// one pool.  The former private my_build_brep pool produced handles that
// resolved only in that private pool and broke every downstream
// index-based consumer (the same_parameter_pass "edge_mut: Shape N is not
// an Edge" defect).

impl TopOpeBRepBuildHBuilder {
    // OCCT TopOpeBRepBuild_HBuilder.hxx L52 — Perform(HDS).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // TopOpeBRepBuild_Builder::Perform (Builder.cxx L116-134):
    //   L121  Clear();
    //   L122  myDataStructure = HDS;
    //   L123  BuildVertices(HDS);
    //   L124  SplitEvisoONperiodicF();   (periodic-face ON split — no
    //         ChFi3d carrier, D6 no-op)
    //   L125  BuildEdges(HDS);
    //   L126  BuildFaces(HDS);
    //   L127  myIsKPart = 0;             (KPart state is not consulted over
    //         the TKBO path)
    //   L128-129 InitSection / SplitSectionEdges; L130-133 Filter / Reducer.
    //         (old-TopOpeBRepDS section-edge and interference-filter
    //         machinery — the TKBO pipeline carries its own section /
    //         filter passes, D6 no-op on this facade)
    //
    // The ChFi3d DS facade (chfi3d_ds.rs) carries the shape table in
    // `bopds`, the geometry payloads (surface / curve / point tables) in
    // `side`, and the interference quadruplets in the facade maps (D6
    // architecture difference: BOPDS has no per-shape quadruplet lists).
    pub fn perform(&mut self, brep: &mut BRep, ds: &mut TopOpeBRepDSHDataStructure) {
        // Builder.cxx L121: Clear();
        self.clear();
        // Builder.cxx L122: myDataStructure = HDS;
        self.my_data_structure = Some(ds.clone());

        // Builder.cxx L123: BuildVertices(HDS) — TopOpeBRepBuild_
        // BuildVertices.cxx L25-36: for every DS point, MakeVertex.
        self.build_vertices(brep, ds);

        // Builder.cxx L125: BuildEdges(HDS) — TopOpeBRepBuild_BuildEdges.cxx
        // L38-100 (per DS curve) and L104-146 (curve walk with the
        // mother-curve filter).
        self.build_edges(brep, ds);

        // Builder.cxx L126: BuildFaces(HDS) — TopOpeBRepBuild_BuildFaces.cxx
        // L99-107: myNewFaces = new HArray1(0, NbSurfaces); for each DS
        // surface BuildFaces(iS, HDS) (L40-95).
        self.build_faces(brep, ds);

        // Builder.cxx L466-475 (the SplitEdge pass over the DS shapes with
        // point interferences, reached from MergeShapes' SplitShapes walk —
        // performed once here so the tables are ready before MergeSolid,
        // mirroring the OCCT fill order Build* -> split tables).
        self.split_ds_edges(brep, ds);
    }

    /// OCCT TopOpeBRepBuild_Builder::Clear (Builder.cxx L182-245) — drops
    /// the split / merge state.  The OCCT body keeps split entries of
    /// section edges (BDS.IsSectionEdge); over the D6 facade the section
    /// edges are ChFi3d-private records, so the clear is unconditional.
    fn clear(&mut self) {
        self.my_split_out.clear();
        self.my_split_in.clear();
        self.my_split_on.clear();
        self.my_merged_out.clear();
        self.my_merged_in.clear();
        self.my_merged_on.clear();
        self.my_new_edges.clear();
        self.my_new_faces.clear();
        self.my_new_vertices.clear();
    }

    /// OCCT TopOpeBRepBuild_Builder::BuildVertices (BuildVertices.cxx
    /// L25-36): for iP = 1..NbPoints, MakeVertex(ChangeNewVertex(iP),
    /// HDS->Point(iP)).  The DS point table routes to the facade `side`
    /// table (D6).
    fn build_vertices(&mut self, brep: &mut BRep, ds: &TopOpeBRepDSHDataStructure) {
        let n = ds.side.points.len() as i32;
        let mut b1 = BRepBuilder::new();
        for ip in 1..=n {
            // BuildTool::MakeVertex(V, DSP): point + tolerance of the DS point.
            let dsp = ds.point(ip);
            let v = brep.add_tvertex_unique(dsp.point());
            b1.update_vertex_tolerance(brep, v.clone(), dsp.tolerance());
            self.my_new_vertices.insert(ip, v);
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::BuildEdges (BuildEdges.cxx L38-100 per
    /// curve, L104-146 walk): for each DS curve without a mother, build the
    /// edge on the curve and split it at the curve's point interferences;
    /// the pieces append to ChangeNewEdges(iC).
    ///
    /// D6 routing: the curve payload + its range live in the facade `side`
    /// table, the point interferences in the facade curve-interference map
    /// (BOPDS has no curve/point payload equivalents).  The OCCT tail
    /// RecomputeCurves / UpdateEdge / RemoveCurve(iC) (BuildEdges.cxx
    /// L77-99) is the replaced-curve remap of the old DS — the facade side
    /// table carries `mother` for the filter and has no removal, so only
    /// the L133-140 mother gate is translated.
    fn build_edges(&mut self, brep: &mut BRep, ds: &TopOpeBRepDSHDataStructure) {
        let nb_curves = ds.side.curves.len() as i32;
        for ic in 1..=nb_curves {
            let c = ds.curve(ic);
            // BuildEdges.cxx L45-53: skip when the 3d curve and both SCIs
            // are null.
            if c.curve.is_none() && c.sci1.is_none() && c.sci2.is_none() {
                continue;
            }
            // BuildEdges.cxx L133-140: curves with a mother are rebuilt by
            // their mother's pass.
            if c.mother != 0 {
                continue;
            }
            let curve = match c.curve.as_ref() {
                Some(cv) => cv.clone(),
                None => continue,
            };
            let (first, last) = (c.first, c.last);

            // BuildEdges.cxx L59-60: the paves = the curve's point
            // interferences (CurvePoint records: DS point index + parameter).
            let mut paves: Vec<(f64, i32)> = ds
                .curve_interferences(ic)
                .iter()
                .filter_map(|i| match i {
                    TopOpeBRepDSInterference::CurvePoint(ci) => Some((ci.parameter, ci.index_g)),
                    _ => None,
                })
                .collect();
            paves.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            // BuildEdges.cxx L69-73: no paves -> no new edges.
            if paves.is_empty() {
                continue;
            }

            // BuildEdges.cxx L52-63: equalpar = PVS.HasEqualParameters()
            // (PaveSet.cxx L363-407: two distinct paves share one parameter
            // within PConfusion); closvert = PVS.ClosedVertices()
            // (PaveSet.cxx L462-487).  On a closed base curve the
            // EdgeBuilder closed loop yields the single closed edge over
            // the base edge's natural domain — BuildTool::MakeEdge
            // (BuildTool.cxx L157-171) builds the base edge on the curve's
            // natural bounds (its setrange branch is a constant-false dead
            // path), so the degenerate DS range plays no part in the
            // product edge.
            let p_confusion = rcad_kernel::core::precision::p_confusion();
            let equal_par = paves[0].0;
            let equalpar = paves
                .iter()
                .all(|(p, _)| (*p - equal_par).abs() < p_confusion);
            if equalpar {
                let [nat_first, nat_last] = curve.default_domain();
                let closed_curve = (curve.point_at(nat_first).distance(curve.point_at(nat_last)))
                    <= CONFUSION;
                if closed_curve {
                    let v = self.pave_vertex(
                        brep,
                        ds,
                        paves[0].1,
                        TopOpeBRepDSKind::Point,
                        &Some(curve.clone()),
                        equal_par,
                    );
                    let piece = brep.add_tedge(
                        Some(curve.clone()),
                        v.clone(),
                        v.clone(),
                        [nat_first, nat_last],
                    );
                    self.my_new_edges.insert(ic, vec![piece]);
                    continue;
                }
            }

            // MakeEdges (Merge.cxx L440-622): consecutive paves bound the
            // new edges; paves sitting on the curve bounds identify the
            // bound vertices, interior paves get the DS-point vertices.
            let mut first_pave: i32 = -1;
            let mut last_pave: i32 = -1;
            let mut interior: Vec<(f64, i32)> = Vec::new();
            for (p, ip) in &paves {
                if (*p - first).abs() <= CONFUSION {
                    first_pave = *ip;
                } else if (*p - last).abs() <= CONFUSION {
                    last_pave = *ip;
                } else if *p > first + CONFUSION && *p < last - CONFUSION {
                    interior.push((*p, *ip));
                }
            }
            let mut bounds: Vec<(f64, i32)> = Vec::with_capacity(interior.len() + 2);
            bounds.push((first, first_pave));
            bounds.extend(interior);
            bounds.push((last, last_pave));

            let mut pieces: Vec<Shape> = Vec::new();
            for w in bounds.windows(2) {
                let (pa, ipa) = w[0];
                let (pb, ipb) = w[1];
                if pb - pa <= CONFUSION {
                    continue;
                }
                let va = self.bound_vertex(brep, ds, ipa, &curve, pa);
                let vb = self.bound_vertex(brep, ds, ipb, &curve, pb);
                // BuildTool::MakeEdge on the DS curve with the piece range.
                let e = brep.add_tedge(Some(curve.clone()), va, vb, [pa, pb]);
                pieces.push(e);
            }
            // BuildEdges.cxx L75: NCollection_List& EL = ChangeNewEdges(iC).
            self.my_new_edges.insert(ic, pieces);
        }
    }

    /// The vertex of a piece boundary: the DS-point vertex for paves
    /// (MakeEdges L563-598), a fresh vertex on the curve for the plain
    /// bounds (BuildTool::MakeEdge bound vertices).  `ds` is kept on the
    /// signature for the OCCT form (MakeEdges reads the HDS point state);
    /// the D6 facade carries no extra bound-vertex payload beyond the
    /// myNewVertices table filled by build_vertices.
    fn bound_vertex(
        &mut self,
        brep: &mut BRep,
        _ds: &TopOpeBRepDSHDataStructure,
        ipoint: i32,
        curve: &rcad_kernel::geom::Curve3,
        param: f64,
    ) -> Shape {
        if ipoint > 0 {
            if let Some(v) = self.my_new_vertices.get(&ipoint) {
                return v.clone();
            }
        }
        let p = curve.point_at(param);
        brep.add_tvertex_unique(p)
    }

    /// OCCT TopOpeBRepBuild_Builder::BuildFaces (BuildFaces.cxx L99-107 per
    /// surface, L40-95 per BuildFaces(iS)): make the face on the DS surface,
    /// feed it the new edges of the surface's curves (SurfaceCurves(iS)),
    /// build the wires and store into ChangeNewFaces(iS).
    ///
    /// D6 routing: SurfaceCurves(iS) = the facade surface-interference
    /// records of iS (SurfaceCurve quadruplets carry curve index + pcurve +
    /// transition), with the SCI index_s scan as the fallback binding.  The
    /// closing curves of the fillet patches reach this surface-curve
    /// binding through the ChFi3d_FilDS emission
    /// (chfi3d_builder_0_filds.rs — OCCT ChFi3d_Builder_0.cxx L2855-2889 /
    /// L2918-2950 / L3171-3240), so each patch boundary curve carries its
    /// pcurve + end paves into the L72-89 walk below.
    ///
    /// L78-86 per start element: the edge tolerance is raised to the DS
    /// surface tolerance (UpdateEdge), the element takes the SCI
    /// orientation (Orientation(TopAbs_IN)) and the SCI pcurve is put on
    /// the edge for the built face (myBuildTool.PCurve).  Architecture
    /// note (FaceBuilder): OCCT TopOpeBRepBuild_FaceBuilder builds closed
    /// areas from the WireEdgeSet.  The rcad equivalent chains the new
    /// edges by shared end vertices; only closed chains become faces —
    /// open chains produce no face, as in OCCT.
    fn build_faces(&mut self, brep: &mut BRep, ds: &TopOpeBRepDSHDataStructure) {
        let nb_surfaces = ds.side.surfaces.len() as i32;
        // BuildFaces.cxx L102: the array spans 0..NbSurfaces — every index
        // is bound (possibly to an empty list).
        for is in 1..=nb_surfaces {
            self.my_new_faces.entry(is).or_default();
        }
        let mut b1 = BRepBuilder::new();
        for is in 1..=nb_surfaces {
            // BuildFaces.cxx L62: the curves of the surface.
            let mut bound_curves: Vec<(i32, Option<rcad_kernel::geom::Curve2d>, Orientation)> =
                Vec::new();
            for i in ds.surface_interferences(is) {
                if let TopOpeBRepDSInterference::SurfaceCurve(sc) = i {
                    if sc.index_g > 0
                        && !bound_curves.iter().any(|(ic, _, _)| *ic == sc.index_g)
                    {
                        // BuildFaces.cxx L84: SCurves.Orientation(TopAbs_IN).
                        bound_curves.push((
                            sc.index_g,
                            sc.pcurve.clone(),
                            sc.transition.orientation_in(),
                        ));
                    }
                }
            }
            if bound_curves.is_empty() {
                // Fallback binding: the SCI records of the side-table curves
                // (no transition carrier — the neutral orientation).
                for ic in 1..=ds.side.curves.len() as i32 {
                    let c = ds.curve(ic);
                    let hit = c
                        .sci1
                        .as_ref()
                        .map(|s| s.index_s == is)
                        .unwrap_or(false)
                        || c.sci2
                            .as_ref()
                            .map(|s| s.index_s == is)
                            .unwrap_or(false);
                    if hit {
                        let pc = c.sci1.as_ref().and_then(|s| s.pcurve.clone());
                        bound_curves.push((ic, pc, Orientation::Forward));
                    }
                }
            }

            // BuildFaces.cxx L72-89: the start elements = NewEdges(iC), each
            // raised to the surface tolerance and oriented by its SCI.
            let a_tbs_tol = ds.surface(is).tolerance();
            let mut start_elements: Vec<Shape> = Vec::new();
            let mut piece_pcurves: Vec<(u64, Option<rcad_kernel::geom::Curve2d>)> = Vec::new();
            for (ic, pc, ori) in &bound_curves {
                if let Some(pieces) = self.my_new_edges.get(ic) {
                    for piece in pieces {
                        // BuildFaces.cxx L78-82: aTBCTol = Tolerance(aE);
                        // if (aTBCTol < aTBSTol) aBB.UpdateEdge(aE, aTBSTol).
                        // In-place edit (the OCCT BRep_Builder::UpdateEdge
                        // semantics): the piece handle recorded in
                        // myNewEdges observes the change.
                        let piece_tol = brep.tolerance(piece);
                        if piece_tol < a_tbs_tol {
                            b1.update_edge_tolerance(brep, piece.clone(), a_tbs_tol);
                        }
                        // BuildFaces.cxx L84: the SCI orientation.
                        let mut oriented = piece.clone();
                        oriented.orientation = *ori;
                        start_elements.push(oriented);
                        piece_pcurves.push((piece.ptr_id(), pc.clone()));
                    }
                }
            }
            if start_elements.is_empty() {
                continue;
            }

            // FaceBuilder equivalent: chain by shared end vertices; closed
            // chains become faces on the DS surface.
            let loops = chain_closed_loops(brep, &start_elements);
            let surface = ds.surface(is).surface().clone();
            let mut faces: Vec<Shape> = Vec::new();
            for loop_edges in loops {
                let w = b1.build_wire(brep, loop_edges.clone());
                let f = b1.make_face(brep, Some(surface.clone()), w);
                // BuildFaces.cxx L85-86: myBuildTool.PCurve(aFace, anEdge,
                // CDS, PC) — the SCI pcurve of the piece's curve is put on
                // the edge for the built face (BRep_Tool::CurveOnSurface
                // identity: (face TShape, location)).  In-place edit: the
                // pieces are already wired above.
                for e in &loop_edges {
                    let pc = piece_pcurves
                        .iter()
                        .find(|(pid, _)| *pid == e.ptr_id())
                        .and_then(|(_, pc)| pc.clone());
                    if let Some(pc) = pc {
                        let [t1, t2] = pc.default_domain();
                        brep.edge_mut_inplace(e.clone())
                            .pcurves
                            .insert((f.ptr_id(), f.location), (pc, t1, t2));
                    }
                }
                faces.push(f);
            }
            self.my_new_faces.insert(is, faces);
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::SplitEdge pass (Builder.cxx L925-1152,
    /// reached from SplitFace1 over the DS faces' edges — the MergeShapes /
    /// SplitShapes walk): a DS edge carrying point interferences is marked
    /// split and its pieces enter mySplitIN (and, per SplitEdge1 L1033, the
    /// ChangeMerged table of the edge).  Performed once at Perform so the
    /// tables are ready for the split-edge tolerance pass of the ChFi3d
    /// compute tail (ChFi3d_Builder.cxx L480-509).
    ///
    /// Per-edge body = SplitEdge1 (Builder.cxx L941-1081): FORWARD edge,
    /// ToSplit gate, PaveSet over EdgePoints (the edge's DS point
    /// interferences) + the edge's own bound vertices (PaveSet::Prepare),
    /// MarkSplit, then MakeEdges (Merge.cxx L516-622) builds one new edge
    /// per consecutive pave pair: CopyEdge (BuildTool.cxx L293-318 — the
    /// copy keeps the source curve, tolerance and pcurve representations)
    /// + AddEdgeVertex/Parameter per pave vertex.
    ///
    /// D6 architecture note (pool identity, converged): OCCT pieces are free
    /// TShape handles; rcad materializes the piece TShapes into the caller's
    /// pool (brep — the same pool the DS shapes and the Build* products
    /// live in), and the recorded piece handles are the real restricted
    /// pieces (range + carried pcurves + the CopyEdge tolerance), so every
    /// index-based consumer resolving through the caller's BRep
    /// (chfi3d_perform.rs L187-200) reads the OCCT piece form directly.
    fn split_ds_edges(&mut self, brep: &mut BRep, ds: &TopOpeBRepDSHDataStructure) {
        let entries: Vec<i32> = ds.shape_interferences.keys().copied().collect();
        for ishape in entries {
            let list = ds.shape_interferences(ishape);
            // Builder.cxx L999-1000: the paves = EdgePoints(Eforward) — the
            // edge's DS point interferences (CurvePoint records: DS
            // point/vertex index + parameter), ordered like the
            // TopOpeBRepDS_PointIterator walk.
            let mut paves: Vec<(f64, i32, TopOpeBRepDSKind)> = list
                .iter()
                .filter_map(|i| match i {
                    TopOpeBRepDSInterference::CurvePoint(ci) => {
                        Some((ci.parameter, ci.index_g, ci.kind_g))
                    }
                    _ => None,
                })
                .collect();
            if paves.is_empty() {
                continue;
            }
            let s = ds.shape(ishape);
            if s.shape_type() != topods::ShapeType::Edge {
                continue;
            }
            // Builder.cxx L948 + L308-328: tosplit = !IsSplit(E, ToBuild1)
            // && (HasGeometry(E) || HasSameDomain(E)).  HasGeometry = the
            // non-empty shape interference list (already true here);
            // HasSameDomain has no D6 facade equivalent (ChFi3d fills no
            // SameDomain tables).
            if self.is_split(s, TopAbsState::In) {
                continue;
            }
            // Builder.cxx L944-947: work on a FORWARD edge.
            let e_forward = {
                let mut e = s.clone();
                e.orientation = Orientation::Forward;
                e
            };
            let ed = e_forward.as_edge().expect("DS edge");
            let (first, last) = (ed.range[0], ed.range[1]);

            // Builder.cxx L999-1000 FillVertexSet (L1991-2053): kind POINT
            // paves resolve to the BuildVertices products (myNewVertices),
            // kind VERTEX paves to the DS vertex shapes.
            paves.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            // Builder.cxx L1010: MarkSplit(Eforward, ToBuild1) — bound even
            // before the pave walk (the no-vertex split marks, L1011-1025).
            let key = (s.ptr_id(), s.location);
            let entry = self.my_split_in.entry(key).or_default();
            entry.is_split = true;

            // Builder.cxx L1028-1034: EdgeBuilder over the pave set;
            // EdgeList = ChangeMerged(Eforward, ToBuild1); MakeEdges
            // (Merge.cxx L516-622) appends one new edge per consecutive
            // pave pair.  The edge's own bound vertices enter the pave set
            // first (PaveSet::Prepare, PaveSet.cxx L113+): a bound-identical
            // interference pave replaces the fresh bound vertex
            // (MakeEdges L553-577 equafound betonnage).
            let curve = ed.curve.clone();
            let (v_first, v_last) = (ed.first.clone(), ed.last.clone());
            let src_tol = ed.tolerance;
            let src_pcurves = ed.pcurves.clone();
            let src_degenerated = ed.degenerated;

            let mut bounds: Vec<(f64, Shape)> = Vec::with_capacity(paves.len() + 2);
            let mut vstart = v_first.clone();
            let mut vend = v_last.clone();
            for (par, ip, kind) in &paves {
                if (*par - first).abs() <= CONFUSION {
                    vstart = self.pave_vertex(brep, ds, *ip, *kind, &curve, *par);
                } else if (*par - last).abs() <= CONFUSION {
                    vend = self.pave_vertex(brep, ds, *ip, *kind, &curve, *par);
                }
            }
            bounds.push((first, vstart));
            for (par, ip, kind) in &paves {
                if *par > first + CONFUSION && *par < last - CONFUSION {
                    bounds.push((*par, self.pave_vertex(brep, ds, *ip, *kind, &curve, *par)));
                }
            }
            bounds.push((last, vend));

            let mut edge_list: Vec<Shape> = Vec::new();
            for w in bounds.windows(2) {
                let (pa, va) = (&w[0].0, &w[0].1);
                let (pb, vb) = (&w[1].0, &w[1].1);
                if pb - pa <= CONFUSION {
                    // Merge.cxx L525-537: a single-vertex loop is dropped.
                    continue;
                }
                // BuildTool.cxx L293-318 CopyEdge: the piece keeps the
                // source curve, tolerance and pcurve representations;
                // Merge.cxx L590-594 AddEdgeVertex/Parameter records the
                // pave vertex parameters.  In-place edit (the OCCT
                // BRep_Builder::UpdateEdge semantics): the recorded piece
                // handle must observe the edits, so the pool slot Arc is
                // mutated in place (edge_mut would copy-on-write and leave
                // the recorded handle on the pre-edit data).
                let piece = brep.add_tedge(curve.clone(), va.clone(), vb.clone(), [*pa, *pb]);
                {
                    let pd = brep.edge_mut_inplace(piece.clone());
                    pd.tolerance = src_tol;
                    pd.degenerated = src_degenerated;
                    for (fk, (pcv, t1, t2)) in &src_pcurves {
                        pd.pcurves.insert(*fk, (pcv.clone(), *t1, *t2));
                    }
                    pd.vertex_params.insert(va.ptr_id(), *pa);
                    pd.vertex_params.insert(vb.ptr_id(), *pb);
                }
                edge_list.push(piece);
            }

            // Builder.cxx L1033: the pieces are the ChangeMerged list of the
            // edge; L1037-1053 (ConnectTo1): ChangeSplit(Ecur, ToBuild1)
            // receives the same list; LE2 stays empty (no SameDomain on the
            // D6 facade, ConnectTo2 = false).
            self.my_merged_in.entry(key).or_default().list_on_state = edge_list.clone();
            self.my_split_in.entry(key).or_default().list_on_state = edge_list;
        }
    }

    /// The pave vertex (Builder.cxx FillVertexSetOnValue L2003-2053): kind
    /// POINT resolves to the BuildVertices product (NewVertex), kind VERTEX
    /// to the DS vertex shape; a pave that resolved to neither (its table
    /// slot is absent) falls back to a fresh vertex on the curve, the
    /// BuildTool::MakeEdge bound-vertex form.
    fn pave_vertex(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        ip: i32,
        kind: TopOpeBRepDSKind,
        curve: &Option<rcad_kernel::geom::Curve3>,
        param: f64,
    ) -> Shape {
        if kind == TopOpeBRepDSKind::Point {
            if let Some(v) = self.my_new_vertices.get(&ip) {
                return v.clone();
            }
        } else if kind == TopOpeBRepDSKind::Vertex && ip >= 1 && ip <= ds.nb_shapes() {
            return ds.shape(ip).clone();
        }
        if let Some(c) = curve {
            let p = c.point_at(param);
            return brep.add_tvertex_unique(p);
        }
        Shape::null()
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L85 — MergeSolid(S, TB).
    //
    // OCCT TopOpeBRepBuild_Builder::MergeSolid (Merge.cxx L366-370) is
    // MergeShapes(S, TB, Snull, TB): the ShellFaceSet fill (SplitShapes
    // over S L260 / SplitSolid L1535-1713 / FillSolid L1946-1950), then
    // MakeSolids (L374-405) appends the rebuilt classified solids to
    // ChangeMerged(S, ToBuild); SplitSolid L1653-1669 additionally marks
    // the split and connects the same SolidList to ChangeSplit(S, ToBuild).
    //
    // D6 ruling: the reconstruction goes through the rcad TKBO pipeline.
    // BOPAlgo_BuilderSolid (bop/algo/builder_solid.rs) carries the
    // MakeSolids semantics — faces to closed shells to solids with the
    // growth/hole classification (PerformShapesToAvoid / PerformLoops /
    // PerformAreas, the BOPAlgo equivalent of the SolidBuilder area walk):
    //   - the ShellFaceSet face fill = the object solid's own faces
    //     (FillShape walk, orientations as-found: Reverse(TB, TB) = false,
    //     Builder.cxx L621-634) plus the intersection surfaces — SplitSolid
    //     L1631-1662 walks the DS SolidSurface interferences of the solid
    //     (the ChFi3d_FilDS SolidSurface records, ChFi3d_Builder_0.cxx
    //     L2687) and adds NewFaces(iS) oriented by the interference
    //     transition (the fillet patches);
    //   - MakeSolids = BuilderSolid::Perform, Areas() = the SolidList.
    pub fn merge_solid(&mut self, brep: &mut BRep, s: &Shape, tb: TopAbsState) {
        let _ = brep;
        // MergeShapes L192-195: myState1/myState2 = ToBuild; the rcad
        // tables key by state below.
        // SplitSolid L1602-1607 (FillSolid over LS1): the object solid's
        // own faces enter the face set with their as-found orientations.
        // Pending boundary (recorded 2026-09-09): OCCT SplitShapes /
        // SplitFace1 (Builder.cxx L1171-1289) rebuilds every face carrying
        // split edges or intersection edges (FillFace + AddIntersectionEdges
        // + FaceBuilder) so the support faces share edge TShapes with the
        // patches; the rcad merge feeds the as-found faces, the patch faces
        // cannot chain by edge identity and BuilderSolid drops them
        // (perform_shapes_to_avoid) — the blend-2 closed-patch front (E3-I
        // item 3 / q4).  A first SplitFace1 translation attempt regressed
        // a3/q1 (loop classification / pcurve-transfer semantics of
        // FillFace+FaceBuilder need the full 1:1 study) and was reverted;
        // re-land it together with FillFace/FaceBuilder proper.
        let mut faces = solid_faces(s);
        // SplitSolid L1631-1662: the DS intersection surfaces of the solid.
        let ds = self.my_data_structure.as_ref().expect("Perform first");
        let isolid = ds.bopds.index(s);
        if isolid >= 0 {
            for i in ds.shape_interferences(isolid as i32 + 1) {
                if let TopOpeBRepDSInterference::SolidSurface(ssi) = i {
                    // SSurfaces.Orientation(ToBuild1) — the SolidSurface
                    // interference transition orientation.
                    let ori = ssi.transition.orientation_in();
                    for f in self.new_faces(ssi.index_s) {
                        let mut face = f;
                        face.orientation = ori;
                        faces.push(face);
                    }
                }
            }
        }

        // SplitSolid L1666-1669: SolidBuilder SOBU(SFS); MakeSolids
        // (Merge.cxx L374-405) builds the classified solids.  The D6
        // carrier is BOPAlgo_BuilderSolid over the facade's BOPDS (the
        // shape/locations registration of Perform; its location table
        // stays identity — the ChFi3d DS registers unlocated shapes).
        let mut bs = BuilderSolid::new(&ds.bopds);
        bs.my_shapes = faces;
        bs.perform();
        // Pool identity (see the struct doc): BuilderSolid produces
        // pool-free solid handles (index 0, the bop pipeline convention);
        // their trees are consumed data-side (solid walkers, STEP export,
        // area), so the plain handles are recorded.
        let pieces: Vec<Shape> = bs.my_solids.clone();

        // SplitSolid L1647-1650: ChangeMerged(S1oriented, ToBuild1) =
        // SolidList; L1653-1669: MarkSplit(Scur, ToBuild1) and
        // ChangeSplit(Scur, ToBuild1) = SolidList (ConnectTo1; LS2 is
        // empty — no SameDomain on the D6 facade).
        let key = (s.ptr_id(), s.location);
        match tb {
            TopAbsState::Out => {
                self.my_merged_out.entry(key).or_default().list_on_state = pieces.clone();
                let entry = self.my_split_out.entry(key).or_default();
                entry.is_split = true;
                entry.list_on_state = pieces;
            }
            TopAbsState::In => {
                self.my_merged_in.entry(key).or_default().list_on_state = pieces.clone();
                let entry = self.my_split_in.entry(key).or_default();
                entry.is_split = true;
                entry.list_on_state = pieces;
            }
            TopAbsState::On => {
                self.my_merged_on.entry(key).or_default().list_on_state = pieces.clone();
                let entry = self.my_split_on.entry(key).or_default();
                entry.is_split = true;
                entry.list_on_state = pieces;
            }
            TopAbsState::Unknown => return,
        }
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L88 — IsSplit(S, ToBuild).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // Builder.cxx L378-410: res = bound(S) && losos.IsSplit() on the
    // mySplit(ToBuild) map; Unknown selects no map (L396-399).
    pub fn is_split(&self, s: &Shape, to_build: TopAbsState) -> bool {
        let map = match to_build {
            TopAbsState::Out => &self.my_split_out,
            TopAbsState::In => &self.my_split_in,
            TopAbsState::On => &self.my_split_on,
            TopAbsState::Unknown => return false,
        };
        // OCCT L401-409.
        match map.get(&(s.ptr_id(), s.location)) {
            Some(losos) => losos.is_split,
            None => false,
        }
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L91 — Splits(S, ToBuild).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // Builder.cxx L414-445: the ListOnState of the mySplit(ToBuild) entry;
    // unbound or Unknown -> the empty list.
    pub fn splits(&self, s: &Shape, to_build: TopAbsState) -> Vec<Shape> {
        let map = match to_build {
            TopAbsState::Out => &self.my_split_out,
            TopAbsState::In => &self.my_split_in,
            TopAbsState::On => &self.my_split_on,
            TopAbsState::Unknown => return Vec::new(),
        };
        // OCCT L438-444.
        map.get(&(s.ptr_id(), s.location))
            .map(|losos| losos.list_on_state.clone())
            .unwrap_or_default()
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L98 — Merged(S, ToBuild).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // Merge.cxx L663-696: the ListOnState of the myMerged(ToBuild) entry;
    // unbound or Unknown -> the empty list.
    pub fn merged(&self, s: &Shape, to_build: TopAbsState) -> Vec<Shape> {
        let map = match to_build {
            TopAbsState::Out => &self.my_merged_out,
            TopAbsState::In => &self.my_merged_in,
            TopAbsState::On => &self.my_merged_on,
            TopAbsState::Unknown => return Vec::new(),
        };
        // Merge.cxx L687-695.
        map.get(&(s.ptr_id(), s.location))
            .map(|losos| losos.list_on_state.clone())
            .unwrap_or_default()
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L105 — NewEdges(I).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // Builder.cxx L265-275: the edges built on DS curve I; unbound ->
    // the empty list (myEmptyShapeList).
    pub fn new_edges(&self, i: i32) -> Vec<Shape> {
        self.my_new_edges.get(&i).cloned().unwrap_or_default()
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L111 — NewFaces(I).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // Builder.cxx L249-253: the faces built on DS surface I (array slot).
    pub fn new_faces(&self, i: i32) -> Vec<Shape> {
        self.my_new_faces.get(&i).cloned().unwrap_or_default()
    }
}

/// FaceBuilder-equivalent chaining: order the start elements into closed
/// loops by shared end-vertex identity (BRepTool::FirstVertex /
/// LastVertex).  Edges matching only at one end are flipped to a
/// reversed-orientation handle (OCCT gives the wire edges their
/// orientation the same way).  Returns the closed loops only.
fn chain_closed_loops(brep: &BRep, edges: &[Shape]) -> Vec<Vec<Shape>> {
    let ends = |e: &Shape| -> ((u64, u32), (u64, u32)) {
        let ed = brep.edge(e.clone());
        let kf = (ed.first.ptr_id(), ed.first.location);
        let kl = (ed.last.ptr_id(), ed.last.location);
        // A REVERSED handle runs from LastVertex to FirstVertex — the
        // flip applied during chaining must be honored here, otherwise
        // the walk never turns a corner (the loop would stay open
        // whenever the wire direction opposes the stored edge range).
        if e.orientation == Orientation::Reversed {
            (kl, kf)
        } else {
            (kf, kl)
        }
    };
    let flip = |e: &Shape| -> Shape {
        let orientation = match e.orientation {
            Orientation::Reversed => Orientation::Forward,
            _ => Orientation::Reversed,
        };
        Shape {
            data: e.data.clone(),
            index: e.index,
            location: e.location,
            orientation,
        }
    };

    let mut pool: Vec<Shape> = edges.to_vec();
    let mut loops: Vec<Vec<Shape>> = Vec::new();
    'seed: while let Some(seed) = pool.first().cloned() {
        let (seed_first, _) = ends(&seed);
        let mut chain: Vec<Shape> = vec![seed];
        pool.remove(0);
        loop {
            let cur = chain.last().expect("non-empty chain").clone();
            let (_, cur_last) = ends(&cur);
            // A closed edge (first vertex == last vertex) closes on itself.
            if cur_last == seed_first {
                loops.push(std::mem::take(&mut chain));
                continue 'seed;
            }
            let mut advanced = false;
            for i in 0..pool.len() {
                let cand = pool[i].clone();
                let (cf, cl) = ends(&cand);
                if cf == cur_last {
                    // continues in the same direction
                    chain.push(cand);
                    pool.remove(i);
                    advanced = true;
                    break;
                }
                if cl == cur_last {
                    // continues reversed: flip to a reversed-orientation
                    // handle, the wire keeps the chain direction
                    let flipped = flip(&cand);
                    chain.push(flipped);
                    pool.remove(i);
                    advanced = true;
                    break;
                }
            }
            if !advanced {
                // Open chain: not a face boundary (see build_faces).
                break;
            }
        }
    }
    loops
}

/// The face fill of the object solid (SplitSolid L1602-1607 FillShape over
/// the solid: shells -> faces; OCCT TopOpeBRepTool_ShapeExplorer over
/// TopAbs_SHELL / TopAbs_FACE).  Data-side walk — the pool-free handle
/// semantics of OCCT exploration.
fn solid_faces(s: &Shape) -> Vec<Shape> {
    let mut faces: Vec<Shape> = Vec::new();
    fn walk(sh: &Shape, faces: &mut Vec<Shape>) {
        match &*sh.data {
            TShape::Solid(sd) => {
                for shell in &sd.shells {
                    walk(shell, faces);
                }
            }
            TShape::Shell(sd) => {
                for f in &sd.faces {
                    faces.push(f.clone());
                }
            }
            TShape::Face(_) => faces.push(sh.clone()),
            TShape::Compound(cs) | TShape::CompSolid(cs) => {
                for c in cs {
                    walk(c, faces);
                }
            }
            _ => {}
        }
    }
    walk(s, &mut faces);
    faces
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fillet::brep_fillet_api::BRepFilletAPIMakeFillet;
    use glam::DVec3;
    use rcad_kernel::topods::TShape;
    use rcad_modeling::make_box_brep;

    fn edge_shapes(brep: &topods::BRep) -> Vec<Shape> {
        brep.tshapes
            .iter()
            .enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Edge(_)))
            .map(|(idx, ts)| Shape::from_parts(ts.clone(), idx, 0, Orientation::Forward))
            .collect()
    }

    fn solid_shapes(brep: &topods::BRep) -> Vec<Shape> {
        brep.tshapes
            .iter()
            .enumerate()
            .filter(|(_, ts)| matches!(ts.as_ref(), TShape::Solid(_)))
            .map(|(idx, ts)| Shape::from_parts(ts.clone(), idx, 0, Orientation::Forward))
            .collect()
    }

    /// Regression anchor for the HBuilder TKBO facade (the box-fillet
    /// smoke): the ChFi3d compute chain over the public surface —
    /// perform_set_of_surf, FilDS, perform_hbuilder_reconstruction — runs
    /// to completion, the SplitEdge pass marks the DS edges carrying point
    /// interferences, and the MergeSolid reconstruction (BOPAlgo_
    /// BuilderSolid) yields a Solid in the merged table.
    ///
    /// The corner stage (perform_fillet_on_vertex) is intentionally not
    /// part of this smoke: the box single-edge spine ends on vertices whose
    /// corner tails are under separate translation, so the stripe's end
    /// curves stay unrecorded (the OCCT IsBound sink semantics of the DS
    /// facade absorb the unset-index accesses).
    #[test]
    fn box_fillet_hbuilder_reconstruction_smoke() {
        let brep =
            make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 100.0, 100.0, 100.0).expect("box");
        let solids = solid_shapes(&brep);
        assert_eq!(solids.len(), 1, "box has one solid");
        let mut mk = BRepFilletAPIMakeFillet::new(&brep, &solids[0]);
        let edges = edge_shapes(&brep);
        mk.add_radius(10.0, &edges[4]);
        // Stage-by-stage replication of ChFi3d_Builder::Compute (chfi3d.rs
        // L424+) over the public call surface (update_tolesp /
        // extent_analyse are private and are skipped).
        let b = &mut mk.my_builder.base;
        b.reset();
        b.my_ds = Some(super::super::chfi3d_ds::TopOpeBRepDSHDataStructure::default());
        b.done = true;
        b.hasresult = false;
        for itel in b.my_list_stripe.clone() {
            assert!(b.perform_set_of_surf(&itel, false), "perform_set_of_surf");
        }
        for itel in b.my_list_stripe.clone() {
            let st = itel.read().expect("stripe lock");
            b.filds_stripe(&st);
            drop(st);
        }
        let mut map_ind_so: Vec<i32> = Vec::new();
        for cursol in crate::fillet::brep_fillet_api::explore_solids(&b.my_brep) {
            let indcursol = b.my_ds.as_mut().expect("DS").add_shape(&cursol);
            if !map_ind_so.contains(&indcursol) {
                map_ind_so.push(indcursol);
            }
        }
        assert_eq!(map_ind_so, vec![1], "the box solid is DS shape 1");
        b.perform_hbuilder_reconstruction(&map_ind_so);

        let coup = b.my_coup.as_ref().expect("coup");
        // The split tables: the spine/support edges carrying point
        // interferences are marked split (SplitEdge pass).
        assert!(
            !coup.my_split_in.is_empty(),
            "SplitEdge pass marks the DS edges with point interferences"
        );
        for e in coup.my_split_in.values() {
            assert!(e.is_split, "MarkSplit(E, TopAbs_IN)");
        }
        // The merge tables: the solid reconstruction yields Solid entries.
        let solid_merges: Vec<&ListOfShapeOn1State> = coup
            .my_merged_in
            .values()
            .filter(|e| e.list_on_state.iter().any(|s| s.shape_type() == topods::ShapeType::Solid))
            .collect();
        assert!(
            !solid_merges.is_empty(),
            "MergeSolid yields a Solid in myMerged(IN)"
        );
        assert!(b.done, "the reconstruction keeps done for the smoke flow");
    }
}
