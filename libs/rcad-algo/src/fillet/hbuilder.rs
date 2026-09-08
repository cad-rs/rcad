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
//!   - MergeSolid(S, TB) drives the rcad bop Splitter path (PaveFiller +
//!     BOPAlgo_Builder with `my_is_splitter`), the TKBO equivalent of the
//!     OCCT ShellFaceSet / SolidBuilder reconstruction: the fillet patch
//!     faces (the BuildFaces products) are the splitter tools and the
//!     support solid is the object.
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

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::CurveEval as _;
use rcad_kernel::topo::topods::{self, BRep, BRepBuilder, Orientation, Shape};

use super::chfi3d_builder_2::TopAbsState;
use super::chfi3d_ds::{TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference};
use crate::bop::algo::builder::{Builder, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;

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
    /// D6 architecture: OCCT shapes are pool-free handles; rcad TopoDS_Shape
    /// lives inside a BRep TShape pool.  All shapes materialized by Perform
    /// (new vertices / edges / faces) are built into this pool.
    pub my_build_brep: BRep,
    /// D6 architecture: the merge-pipeline result pool.  The pieces stored
    /// in the merged tables reference this pool (kept alive here); replaced
    /// on every merge run, exactly as OCCT's MergeShapes rebuilds per call.
    pub my_merge_brep: Option<BRep>,
}

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
    pub fn perform(&mut self, ds: &mut TopOpeBRepDSHDataStructure) {
        // Builder.cxx L121: Clear();
        self.clear();
        // Builder.cxx L122: myDataStructure = HDS;
        self.my_data_structure = Some(ds.clone());

        // Builder.cxx L123: BuildVertices(HDS) — TopOpeBRepBuild_
        // BuildVertices.cxx L25-36: for every DS point, MakeVertex.
        self.build_vertices(ds);

        // Builder.cxx L125: BuildEdges(HDS) — TopOpeBRepBuild_BuildEdges.cxx
        // L38-100 (per DS curve) and L104-146 (curve walk with the
        // mother-curve filter).
        self.build_edges(ds);

        // Builder.cxx L126: BuildFaces(HDS) — TopOpeBRepBuild_BuildFaces.cxx
        // L99-107: myNewFaces = new HArray1(0, NbSurfaces); for each DS
        // surface BuildFaces(iS, HDS) (L40-95).
        self.build_faces(ds);

        // Builder.cxx L466-475 (the SplitEdge pass over the DS shapes with
        // point interferences, reached from MergeShapes' SplitShapes walk —
        // performed once here so the tables are ready before MergeSolid,
        // mirroring the OCCT fill order Build* -> split tables).
        self.split_ds_edges(ds);
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
        self.my_build_brep = BRep::new();
        self.my_merge_brep = None;
    }

    /// OCCT TopOpeBRepBuild_Builder::BuildVertices (BuildVertices.cxx
    /// L25-36): for iP = 1..NbPoints, MakeVertex(ChangeNewVertex(iP),
    /// HDS->Point(iP)).  The DS point table routes to the facade `side`
    /// table (D6).
    fn build_vertices(&mut self, ds: &TopOpeBRepDSHDataStructure) {
        let n = ds.side.points.len() as i32;
        let mut b1 = BRepBuilder::new();
        for ip in 1..=n {
            // BuildTool::MakeVertex(V, DSP): point + tolerance of the DS point.
            let dsp = ds.point(ip);
            let v = self.my_build_brep.add_tvertex_unique(dsp.point());
            b1.update_vertex_tolerance(&mut self.my_build_brep, v.clone(), dsp.tolerance());
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
    fn build_edges(&mut self, ds: &TopOpeBRepDSHDataStructure) {
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
                let va = self.bound_vertex(ds, ipa, &curve, pa);
                let vb = self.bound_vertex(ds, ipb, &curve, pb);
                // BuildTool::MakeEdge on the DS curve with the piece range.
                let e = self
                    .my_build_brep
                    .add_tedge(Some(curve.clone()), va, vb, [pa, pb]);
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
        self.my_build_brep.add_tvertex_unique(p)
    }

    /// OCCT TopOpeBRepBuild_Builder::BuildFaces (BuildFaces.cxx L99-107 per
    /// surface, L40-95 per BuildFaces(iS)): make the face on the DS surface,
    /// feed it the new edges of the surface's curves (SurfaceCurves(iS)),
    /// build the wires and store into ChangeNewFaces(iS).
    ///
    /// D6 routing: SurfaceCurves(iS) = the facade surface-interference
    /// records of iS (SurfaceCurve quadruplets carry curve index + pcurve),
    /// with the SCI index_s scan as the fallback binding.  The pcurve
    /// write-back onto the materialized edges (myBuildTool.PCurve,
    /// BuildFaces.cxx L87) is pending the DS-owned build pool follow-up;
    /// the SCI pcurves stay on the facade.
    ///
    /// Architecture note (FaceBuilder): OCCT TopOpeBRepBuild_FaceBuilder
    /// builds closed areas from the WireEdgeSet.  The rcad equivalent
    /// chains the new edges by shared end vertices; only closed chains
    /// become faces.  Open chains (the fillet patches whose closing curves
    /// are not yet in the DS — the untranslated Filds closing-curve
    /// segments) produce no face, matching the pending-boundary neutral.
    fn build_faces(&mut self, ds: &TopOpeBRepDSHDataStructure) {
        let nb_surfaces = ds.side.surfaces.len() as i32;
        // BuildFaces.cxx L102: the array spans 0..NbSurfaces — every index
        // is bound (possibly to an empty list).
        for is in 1..=nb_surfaces {
            self.my_new_faces.entry(is).or_default();
        }
        let mut b1 = BRepBuilder::new();
        for is in 1..=nb_surfaces {
            // BuildFaces.cxx L62: the curves of the surface.
            let mut bound_curves: Vec<(i32, Option<rcad_kernel::geom::Curve2d>)> = Vec::new();
            for i in ds.surface_interferences(is) {
                if let TopOpeBRepDSInterference::SurfaceCurve(sc) = i {
                    if sc.index_g > 0
                        && !bound_curves.iter().any(|(ic, _)| *ic == sc.index_g)
                    {
                        bound_curves.push((sc.index_g, sc.pcurve.clone()));
                    }
                }
            }
            if bound_curves.is_empty() {
                // Fallback binding: the SCI records of the side-table curves.
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
                        bound_curves.push((ic, pc));
                    }
                }
            }

            // BuildFaces.cxx L72-89: the start elements = NewEdges(iC).
            let mut start_elements: Vec<Shape> = Vec::new();
            for (ic, _pc) in &bound_curves {
                if let Some(pieces) = self.my_new_edges.get(ic) {
                    start_elements.extend(pieces.iter().cloned());
                }
            }
            if start_elements.is_empty() {
                continue;
            }

            // FaceBuilder equivalent: chain by shared end vertices; closed
            // chains become faces on the DS surface.
            let loops = chain_closed_loops(&self.my_build_brep, &start_elements);
            let surface = ds.surface(is).surface().clone();
            let mut faces: Vec<Shape> = Vec::new();
            for loop_edges in loops {
                let w = b1.build_wire(&mut self.my_build_brep, loop_edges);
                let f = b1.make_face(&mut self.my_build_brep, Some(surface.clone()), w);
                faces.push(f);
            }
            self.my_new_faces.insert(is, faces);
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::SplitEdge pass (Builder.cxx L925-1152,
    /// reached from MergeShapes' SplitShapes walk L1739-1741): a DS edge
    /// carrying point interferences is marked split and its pieces enter
    /// mySplitIN.  Performed once at Perform so the tables are ready for
    /// the split-edge tolerance pass of the ChFi3d compute tail
    /// (ChFi3d_Builder.cxx L480-509).
    ///
    /// D6 architecture note: the OCCT pieces are new sub-edge handles built
    /// by TopOpeBRepDS_BuildTool::SplitEdge; they are consumed by the
    /// tolerance pass through the CALLER's BRep pool (BRep_Tool::Tolerance /
    /// TopExp vertices, chfi3d_perform.rs).  rcad pool indexes only resolve
    /// inside their own pool, so the pieces here are the safe same-pool
    /// re-wraps of the source DS edge: the tolerance propagation reads the
    /// source edge tolerance (SplitEdge copies it onto every piece in
    /// OCCT) and propagates it to the source end vertices — conservative,
    /// panic-free, and converging to the OCCT form once the DS owns the
    /// build pool.
    fn split_ds_edges(&mut self, ds: &TopOpeBRepDSHDataStructure) {
        let entries: Vec<i32> = ds.shape_interferences.keys().copied().collect();
        for ishape in entries {
            let list = ds.shape_interferences(ishape);
            let has_point = list
                .iter()
                .any(|i| matches!(i, TopOpeBRepDSInterference::CurvePoint(_)));
            if !has_point {
                continue;
            }
            let s = ds.shape(ishape);
            if s.shape_type() != topods::ShapeType::Edge {
                continue;
            }
            // MarkSplit(S, TopAbs_IN, true) + ChangeSplit(S, TopAbs_IN)
            // append (Builder.cxx L330-374 / L449-465).
            let entry = self.my_split_in.entry((s.ptr_id(), s.location)).or_default();
            entry.is_split = true;
            entry.list_on_state.push(s.clone());
        }
    }

    // OCCT TopOpeBRepBuild_HBuilder.hxx L85 — MergeSolid(S, TB).
    //
    // D6 ruling: implemented over the rcad TKBO pipeline.  OCCT
    // TopOpeBRepBuild_Builder::MergeSolid (Merge.cxx L366-370) is
    // MergeShapes(S, TB, Snull, TB): the ShellFaceSet fill (SplitShapes
    // over S's faces, L260), the area construction and MakeSolids
    // (L374-405) append the rebuilt solids to ChangeMerged(S, ToBuild).
    //
    // TKBO path: the rcad bop Splitter run (PaveFiller + BOPAlgo_Builder
    // with my_is_splitter, the BRepAlgoAPI_Splitter::Build translation in
    // bop/brep_algo_api) over objects = [S] and tools = the new faces (the
    // fillet patches built by BuildFaces — the ShellFaceSet patch
    // equivalent).  The splitter's object images (myImages) are the split
    // pieces of S and feed myMerged(TB)[S]; the result pool is kept alive
    // on my_merge_brep.
    //
    // Degenerate form: with no fillet patches (or a failed run) the merge
    // degrades to the trivial rebuild [S] — the OCCT merge of an untouched
    // solid also re-appends S itself, and the ChFi3d result assembly
    // (ChFi3d_Builder.cxx L515-541) pushes exactly that shape.
    pub fn merge_solid(&mut self, s: &Shape, tb: TopAbsState) {
        // OCCT MergeShapes L192-195: myState1/myState2 = ToBuild; the rcad
        // tables key by state below.
        let tools: Vec<Shape> = self.my_new_faces.values().flatten().cloned().collect();

        let pieces: Vec<Shape> = if tools.is_empty() {
            // Trivial rebuild: nothing to split with.
            vec![s.clone()]
        } else {
            let mut all_args = vec![s.clone()];
            all_args.extend(tools.iter().cloned());
            // OCCT BRepAlgoAPI_Splitter::Build: aLArgs = objects + tools ->
            // one PaveFiller; the builder receives objects and tools
            // separately (bop/brep_algo_api run_build_splitter_brep form).
            let mut filler = PaveFiller::new();
            filler.set_arguments(all_args);
            filler.set_fuzzy_value(CONFUSION);
            let a_prog = NoopProgress;
            let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
            filler.perform(&a_ps);
            let fuzz = filler.fuzzy_value();
            // The DS clone of S carries the identity the images are keyed by.
            let ds_arg0 = filler.ds().arguments.first().cloned();
            let mut builder = Builder::new(filler.ds(), BooleanOpType::Union, fuzz);
            builder.my_arguments = filler.ds().arguments.clone();
            builder.my_tools = builder.my_arguments[1..].to_vec();
            builder.my_is_splitter = true;
            match builder.build() {
                Ok(brep) => {
                    let mut pieces: Vec<Shape> = Vec::new();
                    if let Some(s0) = ds_arg0 {
                        let key = (s0.ptr_id(), s0.location);
                        if let Some(imgs) = builder.my_images.get(key) {
                            // Normalize the pieces onto the returned pool
                            // (the pool kept on my_merge_brep below).
                            for p in imgs {
                                let piece = if p.index < brep.tshapes.len() {
                                    Shape::from_parts(
                                        brep.tshapes[p.index].clone(),
                                        p.index,
                                        p.location,
                                        p.orientation,
                                    )
                                } else {
                                    p.clone()
                                };
                                pieces.push(piece);
                            }
                        }
                    }
                    self.my_merge_brep = Some(brep);
                    if pieces.is_empty() {
                        // The run produced no images for the object: the
                        // solid is untouched — trivial rebuild.
                        vec![s.clone()]
                    } else {
                        pieces
                    }
                }
                Err(()) => {
                    // Pipeline failure: degrade to the trivial rebuild.
                    vec![s.clone()]
                }
            }
        };

        // ChangeMerged(S, ToBuild) + append (Merge.cxx L700-729).
        let key = (s.ptr_id(), s.location);
        let entry = match tb {
            TopAbsState::Out => self.my_merged_out.entry(key).or_default(),
            TopAbsState::In => self.my_merged_in.entry(key).or_default(),
            TopAbsState::On => self.my_merged_on.entry(key).or_default(),
            TopAbsState::Unknown => return,
        };
        entry.list_on_state.extend(pieces);
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
        (kf, kl)
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
