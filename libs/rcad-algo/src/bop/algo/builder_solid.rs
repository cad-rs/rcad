// OCCT BOPAlgo_BuilderSolid — solid building from shells.
//
// OCCT BOPAlgo_BuilderSolid.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::{Alert, Report};
use crate::bop::algo::shell_splitter::{make_connexity_blocks, make_shells};
use crate::bop::closest_point_on_surface;
use crate::bop::ds::DS;
use crate::bop::int_tools::bean_face_intersector::{
    BeanFaceIntersector, BRepAdaptorCurve, BRepAdaptorSurface,
};
use crate::bop::int_tools::context::IntToolsContext;
use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;
use indexmap::IndexMap;
use rcad_kernel::geom::{Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, Orientation, ShapeType, TShape, TSolidData, TShellData, TFaceData, TWireData, tshape_flags};
use std::collections::{HashSet, HashMap, VecDeque};
use std::sync::Arc;
use glam::DVec3;

/// OCCT BOPAlgo_BuilderSolid — builds solids from a set of faces.
pub struct BuilderSolid<'a> {
    ds: &'a DS,
    my_report: Report,
    // BOPAlgo_Algo (inherited)
    my_run_parallel: bool,
    my_context: IntToolsContext,         // OCCT: myContext
    // BOPAlgo_BuilderSolid
    pub my_shapes: Vec<Shape>,          // OCCT: myShapes
    pub my_solids: Vec<Shape>,          // OCCT: myAreas (Areas())
    // OCCT BOPAlgo_SplitSolid::mySolid (Builder_3.cxx L404) — the solid to
    // split. BOPAlgo_BuilderSolid::Perform does not use it; the SplitSolid
    // subclass stores it for Solid() (used as the aSolidsIm key in
    // BuildSplitSolids L542). rcad stores the source solid here as well.
    pub my_solid: Option<Shape>,        // OCCT: BOPAlgo_SplitSolid::mySolid
    my_avoid_internal_shapes: bool,     // OCCT: BOPAlgo_BuilderSolid::myAvoidInternalShapes
    my_shapes_to_avoid: IndexMap<(u64, u32, Orientation), Shape>, // OCCT: NCollection_IndexedMap<TopoDS_Shape> myShapesToAvoid — default hasher, orientation-sensitive
    my_loops: Vec<Vec<Shape>>,          // OCCT: myLoops
    my_loops_internal: Vec<Shape>,      // OCCT: myLoopsInternal (shells with Closed flag)
}

impl<'a> BuilderSolid<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderSolid {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_context: IntToolsContext::new(),
            my_shapes: Vec::new(),
            my_solids: Vec::new(),
            my_solid: None,
            my_avoid_internal_shapes: false,
            my_shapes_to_avoid: IndexMap::new(),
            my_loops: Vec::new(),
            my_loops_internal: Vec::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::Perform (BOPAlgo_BuilderSolid.cxx L76-125).
    pub fn perform(&mut self) {
        // OCCT L80: GetReport()->Clear().
        self.my_report.clear();
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT L106: PerformShapesToAvoid
        self.perform_shapes_to_avoid();
        if self.has_errors() { return; }
        // OCCT L112: PerformLoops — group faces into shells
        self.perform_loops();
        if self.has_errors() { return; }
        // OCCT L118: PerformAreas — classify shells, build solids
        self.perform_areas();
        if self.has_errors() { return; }
        // OCCT L124: PerformInternalShapes
        self.perform_internal_shapes();
    }


    /// OCCT BOPAlgo_BuilderSolid::SetAvoidInternalShapes (BuilderSolid.hxx
    /// L50-55): avoid creating internal shells in the resulting solids.
    pub fn set_avoid_internal_shapes(&mut self, the_flag: bool) {
        self.my_avoid_internal_shapes = the_flag;
    }

    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }

    /// Access to the report (OCCT BOPAlgo_Algo::GetReport) for merging
    /// sub-split alerts into the parent Builder report.
    pub fn report(&self) -> &Report {
        &self.my_report
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformShapesToAvoid (BuilderSolid.cxx L129-220).
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT L138: myShapesToAvoid.Clear()
        self.my_shapes_to_avoid.clear();
        // OCCT L142-218: iterative — mark faces with free boundary edges,
        // repeat until no new faces are marked.
        loop {
            let mut b_found = false;
            // OCCT L151-160: aMEF — edge -> [faces] built via
            // TopExp::MapShapesAndAncestors(aF, TopAbs_EDGE, TopAbs_FACE,
            // aMEF) over every non-avoided face. OCCT aMEF is
            // IndexedDataMap with TopTools_ShapeMapHasher — key identity
            // TShape + Location, insertion order (IndexMap reproduces it).
            let mut a_mef: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            for face in &self.my_shapes {
                if self.my_shapes_to_avoid.contains_key(&(face.ptr_id(), face.location, face.orientation)) { continue; }
                for a_e in face_edges(face) {
                    let ekey = (a_e.ptr_id(), a_e.location);
                    let entry = a_mef.entry(ekey).or_insert_with(|| (a_e.clone(), Vec::new()));
                    entry.1.push(face.clone());
                }
            }
            // OCCT L164-211: per-edge decisions. myShapesToAvoid stores FACES
            // here (a face is only added to it later, in PerformLoops L326-329).
            let keys: Vec<(u64, u32)> = a_mef.keys().copied().collect();
            for ekey in keys {
                let (a_e, a_lf) = &a_mef[&ekey];
                // OCCT L167-170: skip degenerated edges.
                if a_e.as_edge().map_or(false, |ed| ed.degenerated) {
                    continue;
                }
                let a_nb_f = a_lf.len();
                if a_nb_f == 0 { continue; }
                // OCCT L179: aOrE = aE.Orientation().
                let a_or_e = a_e.orientation;
                let a_f1 = &a_lf[0];
                if a_nb_f == 1 {
                    // OCCT L182-190: single face on the edge; the edge is a
                    // free boundary of the face.
                    if a_or_e == topods::Orientation::Internal {
                        continue;
                    }
                    b_found = true;
                    self.my_shapes_to_avoid
                        .entry((a_f1.ptr_id(), a_f1.location, a_f1.orientation))
                        .or_insert_with(|| a_f1.clone());
                } else if a_nb_f == 2 {
                    // OCCT L191-209: two faces on the edge; avoid both copies
                    // when they are IsSame (same TShape + Location), unless
                    // the edge is a seam of the face (BRep_Tool::IsClosed)
                    // or the edge is INTERNAL.
                    let a_f2 = &a_lf[1];
                    if a_f2.is_partner(a_f1) {
                        if edge_closed_on_face(a_e, a_f1) {
                            continue;
                        }
                        if a_or_e == topods::Orientation::Internal {
                            continue;
                        }
                        b_found = true;
                        self.my_shapes_to_avoid
                            .entry((a_f1.ptr_id(), a_f1.location, a_f1.orientation))
                            .or_insert_with(|| a_f1.clone());
                        self.my_shapes_to_avoid
                            .entry((a_f2.ptr_id(), a_f2.location, a_f2.orientation))
                            .or_insert_with(|| a_f2.clone());
                    }
                }
            }
            if !b_found { break; }
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformLoops (BuilderSolid.cxx L237-393).
    /// Groups non-avoided faces into shells via the ShellSplitter, then builds
    /// the internal shells (myLoopsInternal) from the avoided faces.
    fn perform_loops(&mut self) {
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT L230: myLoops.Clear().
        self.my_loops.clear();
        // OCCT L238-257: infinite faces become single-face shells; avoided
        // faces are excluded from the splitter start elements.
        let mut a_start: Vec<Shape> = Vec::new();
        for face in &self.my_shapes {
            if Self::is_infinite_face(face) {
                // OCCT L243-250: MakeShell(aSh); Add(aSh, aF); myLoops.Append.
                self.my_loops.push(vec![face.clone()]);
                continue;
            }
            // OCCT L253-256: if !myShapesToAvoid.Contains(aF).
            if !self.my_shapes_to_avoid.contains_key(&(face.ptr_id(), face.location, face.orientation)) {
                a_start.push(face.clone());
            }
        }
        // OCCT L259-284: ShellSplitter.Perform -> shells appended to myLoops.
        let blocks = make_connexity_blocks(&a_start);
        let shells = make_shells(&blocks, self.ds);
        for shell in shells {
            self.my_loops.push(shell);
        }
        // OCCT L287-331: post-treatment — collect all faces of the loops and
        // of myShapesToAvoid; add the remaining faces to myShapesToAvoid.
        // OCCT aMP (L297) is NCollection_Map<TopoDS_Shape> — default hasher
        // includes orientation. The faces of myLoops may have been flipped by
        // orient_faces_on_shell (ShellSplitter), so aMP must be keyed by
        // (ptr_id, location, orientation): a face whose orientation changed in
        // the splitter output no longer matches its myShapes counterpart and
        // is added to myShapesToAvoid (OCCT L326-329).
        let mut a_mp: HashSet<(u64, u32, Orientation)> = HashSet::new();
        for loop_faces in &self.my_loops {
            for f in loop_faces {
                a_mp.insert((f.ptr_id(), f.location, f.orientation));
            }
        }
        // OCCT L312-317: myShapesToAvoid faces carry their stored orientation.
        for face in &self.my_shapes {
            if self.my_shapes_to_avoid.contains_key(&(face.ptr_id(), face.location, face.orientation)) {
                a_mp.insert((face.ptr_id(), face.location, face.orientation));
            }
        }
        for face in &self.my_shapes {
            if !Self::is_infinite_face(face) {
                if !a_mp.contains(&(face.ptr_id(), face.location, face.orientation)) {
                    self.my_shapes_to_avoid
                        .entry((face.ptr_id(), face.location, face.orientation))
                        .or_insert_with(|| face.clone());
                }
            }
        }
        // OCCT L338-392: build internal shells from the avoided faces.
        // OCCT L339: myLoopsInternal.Clear().
        self.my_loops_internal.clear();
        // OCCT L344-349: edge -> [faces] map from the avoided faces.
        // OCCT aEFMap (L300) is IndexedDataMap with TopTools_ShapeMapHasher —
        // key identity TShape + Location.
        let mut a_ef_map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        // OCCT L368-373: iterate myShapesToAvoid (L368: aNbSh =
        // myShapesToAvoid.Extent(); L371: aFF = myShapesToAvoid(i)).
        for (_, a_ff) in &self.my_shapes_to_avoid {
            for e in face_edges(a_ff) {
                a_ef_map.entry((e.ptr_id(), e.location)).or_default().push(a_ff.clone());
            }
        }
        // OCCT L351-391: grow a shell per avoided face via shared edges.
        // OCCT AddedFacesMap (L296) — NCollection_Map<TopoDS_Shape>, default
        // hasher (TShape + Location + Orientation).
        let mut a_added: HashSet<(u64, u32, Orientation)> = HashSet::new();
        // OCCT L375-416: iterate myShapesToAvoid (L375-381: aNbSh =
        // myShapesToAvoid.Extent(); aFF = myShapesToAvoid(i)).
        for (_, a_ff) in &self.my_shapes_to_avoid {
            if !a_added.insert((a_ff.ptr_id(), a_ff.location, a_ff.orientation)) { continue; }
            let mut a_shell: Vec<Shape> = vec![a_ff.clone()];
            let mut i = 0;
            while i < a_shell.len() {
                let a_f = &a_shell[i];
                for e in face_edges(a_f) {
                    if let Some(a_lf) = a_ef_map.get(&(e.ptr_id(), e.location)) {
                        for a_fl in a_lf {
                            if a_added.insert((a_fl.ptr_id(), a_fl.location, a_fl.orientation)) {
                                a_shell.push(a_fl.clone());
                            }
                        }
                    }
                }
                i += 1;
            }
            // OCCT L414: aShell.Closed(BRep_Tool::IsClosed(aShell)).
            let mut flags = tshape_flags::DEFAULT;
            if crate::bop::algo::builder::shell_is_closed(&a_shell) {
                flags |= tshape_flags::CLOSED;
            }
            let a_shell_shape = Shape::new(
                Arc::new(TShape::Shell(TShellData {
                    my_shapes: vec![],
                    flags,
                    faces: a_shell,
                })),
                0,
                Orientation::Forward,
            );
            self.my_loops_internal.push(a_shell_shape);
        }
    }

    /// OCCT IntTools_Context::IsInfiniteFace — a face without a bounded outer
    /// wire is treated as infinite (unbounded surface).
    fn is_infinite_face(face: &Shape) -> bool {
        match &*face.data {
            TShape::Face(fd) => match &*fd.outer_wire.data {
                TShape::Wire(wd) => wd.edges.is_empty(),
                _ => true,
            },
            _ => true,
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformAreas (BuilderSolid.cxx L397-598).
    /// Classifies shells as growth (solid) or hole (void), then assigns holes
    /// to the innermost growth solid containing them (internal shells). Holes
    /// outside every solid become separate solids.
    fn perform_areas(&mut self) {
        // OCCT L399: myAreas.Clear().
        self.my_solids.clear();
        // OCCT L402-407: the new solids / hole shells / hole-face map.
        let mut new_solids: Vec<Shape> = Vec::new();
        let mut hole_shells: Vec<Vec<Shape>> = Vec::new();
        let mut a_mhf: HashSet<(u64, u32)> = HashSet::new();

        // OCCT L413-442: classify each shell.
        for loop_faces in &self.my_loops {
            if std::env::var("RCAD_BS_DEBUG").is_ok() {
                eprintln!("[BS-SHELL-CAND] n_faces={} faces={:?}", loop_faces.len(),
                    loop_faces.iter().map(|f| {
                        let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                            rcad_kernel::geom::Surface3::Plane(p) => format!("Plane n=({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                            _ => "O".into(),
                        }).unwrap_or_else(|| "?".into());
                        format!("{}:{}", f.ptr_id(), n)
                    }).collect::<Vec<_>>());
            }
            // OCCT L422-427: IsGrowthShell then IsHole.
            let mut b_is_growth = Self::is_growth_shell(loop_faces, &a_mhf);
            if !b_is_growth {
                b_is_growth = !Self::is_hole(loop_faces, &self.ds.locations);
            }
            if b_is_growth {
                // OCCT L430-435: MakeSolid + Add(aSolid, aShell).
                let shell_shape = self.build_shell_shape(loop_faces);
                let solid = self.build_solid_shape(&shell_shape);
                new_solids.push(solid);
            } else {
                // OCCT L437-441: add to hole shells + MapShapes(aShell, FACE, aMHF).
                // OCCT aMHF (L411) is IndexedMap with TopTools_ShapeMapHasher —
                // key identity TShape + Location.
                for f in loop_faces {
                    a_mhf.insert((f.ptr_id(), f.location));
                }
                hole_shells.push(loop_faces.clone());
            }
        }

        // OCCT L444-458: no holes — all growths are the result.
        if hole_shells.is_empty() {
            self.my_solids = new_solids;
            return;
        }

        // OCCT L460-530: classify holes relative to the solids. For each hole
        // find the innermost growth solid containing it (IsInside), keeping the
        // smaller solid when two nested growths both contain the hole.
        // OCCT L462-478: BVH tree of the hole-shell boxes (BRepBndLib::Add);
        // OCCT L493-502: each solid's box; OCCT L499-509: the BoxTreeSelector
        // pre-filters the hole shells — only those whose box interferes with
        // the solid's box reach the IsInside test.
        // OCCT L462-478: BVH tree of the hole-shell boxes (BRepBndLib::Add);
        // OCCT L493-502: each solid's box; OCCT L499-509: the BoxTreeSelector
        // pre-filters the hole shells — only those whose box interferes with
        // the solid's box reach the IsInside test.
        // OCCT aHoleShells holds the hole SHELL shapes; rcad builds each hole
        // shell shape once and reuses it as the aHoleSolidMap key (OCCT L505:
        // IndexedDataMap<TopoDS_Shape, Solid> keyed by the shell TShape).
        let mut hole_shell_shapes: Vec<Shape> = Vec::new();
        for hs in &hole_shells {
            hole_shell_shapes.push(self.build_shell_shape(hs));
        }
        let mut hole_boxes: Vec<Option<(DVec3, DVec3)>> = Vec::new();
        for shell_shape in &hole_shell_shapes {
            hole_boxes.push(crate::bop::algo::builder::shape_bbox(shell_shape));
        }
        let mut a_hole_solid: HashMap<(u64, u32), usize> = HashMap::new();
        for (si, solid) in new_solids.iter().enumerate() {
            let solid_box = crate::bop::algo::builder::shape_bbox(solid);
            for (hi, hs) in hole_shells.iter().enumerate() {
                // OCCT L506-509: BVH pre-filter — skip holes whose box does
                // not interfere with the solid's box (IsInside would be OUT).
                if let (Some((smin, smax)), Some((hmin, hmax))) = (solid_box, hole_boxes[hi]) {
                    if hmax.x < smin.x
                        || hmin.x > smax.x
                        || hmax.y < smin.y
                        || hmin.y > smax.y
                        || hmax.z < smin.z
                        || hmin.z > smax.z
                    {
                        continue;
                    }
                }
                // OCCT L511: if (!IsInside(aHole, aSolid)) continue;
                if !Self::is_inside(hs, solid, self.ds) { continue; }
                // OCCT aHoleSolidMap key — the hole shell shape identity.
                let hkey = (hole_shell_shapes[hi].ptr_id(), hole_shell_shapes[hi].location);
                match a_hole_solid.get(&hkey) {
                    None => {
                        // OCCT L526-528: aHoleSolidMap.Add(aHole, aSolid).
                        a_hole_solid.insert(hkey, si);
                    }
                    Some(&prev) => {
                        // OCCT L520-523: if (IsInside(aSolid, *pSolidWas)) *pSolidWas = aSolid.
                        if Self::is_inside_solid(&new_solids[si], &new_solids[prev], self.ds) {
                            a_hole_solid.insert(hkey, si);
                        }
                    }
                }
            }
        }

        // OCCT L532-548: back map solids -> holes.
        let mut a_solid_holes: HashMap<usize, Vec<usize>> = HashMap::new();
        for (hi, _hs) in hole_shells.iter().enumerate() {
            let hkey = (hole_shell_shapes[hi].ptr_id(), hole_shell_shapes[hi].location);
            if let Some(&si) = a_hole_solid.get(&hkey) {
                a_solid_holes.entry(si).or_default().push(hi);
            }
        }

        // OCCT L550-576: add holes to their solids and append to myAreas.
        for si in 0..new_solids.len() {
            if let Some(holes) = a_solid_holes.get(&si) {
                for &hi in holes {
                    // OCCT L564-568: aBB.Add(aSolid, aHole).
                    let shell_shape = self.build_shell_shape(&hole_shells[hi]);
                    let ts = Arc::make_mut(&mut new_solids[si].data);
                    if let TShape::Solid(sd) = ts {
                        sd.shells.push(shell_shape);
                    }
                }
            }
            self.my_solids.push(new_solids[si].clone());
        }

        // OCCT L578-597: holes outside every solid become separate solids.
        for (hi, hs) in hole_shells.iter().enumerate() {
            let hkey = (hole_shell_shapes[hi].ptr_id(), hole_shell_shapes[hi].location);
            if !a_hole_solid.contains_key(&hkey) {
                let shell_shape = self.build_shell_shape(hs);
                let solid = self.build_solid_shape(&shell_shape);
                self.my_solids.push(solid);
            }
        }
    }

    /// OCCT IsGrowthShell (BuilderSolid.cxx L864-879) — the shell contains a
    /// previously-marked hole face, so it bounds a hole from the outside.
    fn is_growth_shell(faces: &[Shape], the_mhf: &HashSet<(u64, u32)>) -> bool {
        if the_mhf.is_empty() {
            return false;
        }
        faces.iter().any(|f| the_mhf.contains(&(f.ptr_id(), f.location)))
    }

    /// OCCT IsHole (BuilderSolid.cxx L823-831): the shell is classified as a
    /// solid via BRepClass3d_SolidClassifier::PerformInfinitePoint; the point
    /// at infinity is IN => the shell is a hole in space.
    fn is_hole(faces: &[Shape], locations: &[glam::DAffine3]) -> bool {
        let shell_tshape = TShape::Shell(TShellData {
            my_shapes: vec![],
            flags: tshape_flags::DEFAULT | tshape_flags::CLOSED,
            faces: faces.to_vec(),
        });
        let shell = Shape::new(Arc::new(shell_tshape), 0, Orientation::Forward);
        let mut clsf = SolidClassifier::from_shape_with_locations(&shell, locations);
        clsf.perform_infinite_point(f64::MIN_POSITIVE); // OCCT ::RealSmall() = DBL_MIN
        clsf.state() == 3 // TopAbs_IN
    }

    /// OCCT IsInside (BuilderSolid.cxx L835-860) — classify the first face of
    /// the shell relative to the solid via BOPTools_AlgoTools::ComputeState.
    fn is_inside(faces: &[Shape], solid: &Shape, ds: &DS) -> bool {
        let Some(a_f) = faces.first() else {
            // OCCT L869-874: no faces in the shell — PerformInfinitePoint on
            // the solid; State() == IN means the solid is a hole in space.
            let mut clsf = SolidClassifier::from_shape(solid);
            clsf.perform_infinite_point(f64::MIN_POSITIVE); // OCCT ::RealSmall()
            return clsf.state() == 3; // TopAbs_IN
        };
        Self::compute_state_on_solid(a_f, solid, ds) == 3 // TopAbs_IN
    }

    /// IsInside(aSolid, bSolid) — classify the first face of solid a against b.
    fn is_inside_solid(a: &Shape, b: &Shape, ds: &DS) -> bool {
        let Some(a_f) = Self::first_face(a) else {
            return false;
        };
        Self::compute_state_on_solid(&a_f, b, ds) == 3 // TopAbs_IN
    }

    /// OCCT BOPTools_AlgoTools::ComputeState(aF, aSolid, aTol, aBounds)
    /// (BOPTools_AlgoTools.cxx L660-715): try an edge of the face not on the
    /// solid (classify the edge midpoint), else classify a point inside the
    /// face (PointInFace hatcher + PointNearEdge fallback). Delegates to the
    /// shared builder::compute_state_face (the OCCT-aligned implementation).
    fn compute_state_on_solid(a_f: &Shape, solid: &Shape, ds: &DS) -> u8 {
        crate::bop::algo::builder::compute_state_face(a_f, solid, rcad_kernel::CONFUSION, ds)
    }

    /// The first face of a shape (OCCT TopExp_Explorer(shape, TopAbs_FACE)).
    fn first_face(shape: &Shape) -> Option<Shape> {
        let mut stack: Vec<Shape> = vec![shape.clone()];
        while let Some(sh) = stack.pop() {
            match &*sh.data {
                TShape::Solid(sd) => {
                    for x in &sd.shells {
                        stack.push(x.clone());
                    }
                }
                TShape::CompSolid(cd) => {
                    for x in cd {
                        stack.push(x.clone());
                    }
                }
                TShape::Compound(cd) => {
                    for x in cd {
                        stack.push(x.clone());
                    }
                }
                TShape::Shell(sd) => {
                    for x in &sd.faces {
                        stack.push(x.clone());
                    }
                }
                TShape::Face(_) => return Some(sh),
                _ => {}
            }
        }
        None
    }

    /// Build a Shell TShape from a set of faces.
    /// OCCT: the shells come from myLoops, which are already Closed(true)
    /// (BOPAlgo_ShellSplitter::MakeShells L647/L675 — every produced shell,
    /// regular or split, is marked closed).
    fn build_shell_shape(&self, faces: &[Shape]) -> Shape {
        let shell_tshape = TShape::Shell(TShellData {
            my_shapes: vec![],
            flags: tshape_flags::DEFAULT | tshape_flags::CLOSED,
            faces: faces.to_vec(),
        });
        Shape::new(Arc::new(shell_tshape), 0, rcad_kernel::topods::Orientation::Forward)
    }

    /// Build a Solid TShape containing a shell.
    fn build_solid_shape(&self, shell: &Shape) -> Shape {
        let solid_tshape = TShape::Solid(TSolidData {
            my_shapes: vec![],
            flags: tshape_flags::DEFAULT,
            shells: vec![shell.clone()],
            internal_vertices: vec![],
            internal_edges: vec![],
        });
        Shape::new(Arc::new(solid_tshape), 0, rcad_kernel::topods::Orientation::Forward)
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformInternalShapes (BuilderSolid.cxx L602-759).
    fn perform_internal_shapes(&mut self) {
        // OCCT L604-608: user-defined option to avoid internal parts.
        if self.my_avoid_internal_shapes {
            return;
        }
        // OCCT L610-614: no internal parts -> return.
        if self.my_loops_internal.is_empty() {
            return;
        }
        // OCCT L619-629: collect all faces of the internal shells into aMFs.
        // OCCT aMFs (L616) is IndexedMap with TopTools_ShapeMapHasher —
        // key identity TShape + Location.
        let mut a_mfs: Vec<Shape> = Vec::new();
        let mut a_mfs_fence: HashSet<(u64, u32)> = HashSet::new();
        for shell in &self.my_loops_internal {
            // OCCT L648-652: TopoDS_Iterator aIt(aShell) — direct sub-shapes
            // (the faces of the shell).
            if let TShape::Shell(sd) = &*shell.data {
                for f in &sd.faces {
                    if a_mfs_fence.insert((f.ptr_id(), f.location)) {
                        a_mfs.push(f.clone());
                    }
                }
            }
        }
        // OCCT L632-651: no areas — make a solid of the internal faces.
        if self.my_solids.is_empty() {
            let a_lsi = Self::make_internal_shells(&a_mfs);
            let solid_tshape = TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                shells: a_lsi,
                internal_vertices: vec![],
                internal_edges: vec![],
            });
            let solid = Shape::new(Arc::new(solid_tshape), 0, topods::Orientation::Forward);
            self.my_solids.push(solid);
            return;
        }
        // OCCT L673-681: classify the internal faces relative to the areas via
        // BOPAlgo_Tools::ClassifyFaces. OCCT passes myBoxes (the boxes of the
        // areas built in PerformAreas) as theShapeBoxMap; rcad has no box cache,
        // so the empty map makes classify_faces build the boxes itself
        // (equivalent to OCCT L1699-1711). theSolidsIF is an empty map here
        // (OCCT L680).
        let a_mslf_map = crate::bop::algo::builder::Builder::classify_faces(
            self.ds,
            &a_mfs,
            &self.my_solids,
            &HashMap::new(),
            &HashMap::new(),
        );
        // OCCT aMSLF is NCollection_IndexedDataMap — insertion order; IndexMap
        // reproduces it (a HashMap would iterate in random order).
        let mut a_mslf: IndexMap<usize, Vec<Shape>> = IndexMap::new();
        for (si, solid) in self.my_solids.iter().enumerate() {
            if let Some(lf) = a_mslf_map.get(&(solid.ptr_id(), solid.location)) {
                a_mslf.insert(si, lf.clone());
            }
        }
        // OCCT L685-722: update the solids by their internal faces.
        // OCCT aMFDone (L699) — NCollection_Map<TopoDS_Shape,
        // TopTools_ShapeMapHasher>, key TShape + Location.
        let mut a_mf_done: HashSet<(u64, u32)> = HashSet::new();
        for (si, a_lf) in a_mslf {
            if a_lf.is_empty() {
                continue;
            }
            for a_f in &a_lf {
                a_mf_done.insert((a_f.ptr_id(), a_f.location));
            }
            let a_lsi = Self::make_internal_shells(&a_lf);
            let ts = Arc::make_mut(&mut self.my_solids[si].data);
            if let TShape::Solid(sd) = ts {
                sd.shells.extend(a_lsi);
            }
        }
        // OCCT L724-758: warn about the unclassified faces.
        let mut a_mf_unused: Vec<Shape> = Vec::new();
        for a_f in &a_mfs {
            if !a_mf_done.contains(&(a_f.ptr_id(), a_f.location)) {
                a_mf_unused.push(a_f.clone());
            }
        }
        if !a_mf_unused.is_empty() {
            // OCCT L757: AddWarning(BOPAlgo_AlertSolidBuilderUnusedFaces).
            let idxs: Vec<usize> = a_mf_unused.iter().map(|f| f.index).collect();
            self.my_report
                .add_warning(Alert::SolidBuilderUnusedFaces(idxs));
        }
    }

    /// OCCT MakeInternalShells (BuilderSolid.cxx L763-819) — build internal
    /// shells from a set of faces (each face set to INTERNAL orientation).
    fn make_internal_shells(a_mf: &[Shape]) -> Vec<Shape> {
        // OCCT L773-778: edge -> [faces] map. OCCT aMEF (L772) is
        // IndexedDataMap with TopTools_ShapeMapHasher — key TShape + Location.
        let mut a_mef: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        for a_f in a_mf {
            for e in face_edges(a_f) {
                a_mef.entry((e.ptr_id(), e.location)).or_default().push(a_f.clone());
            }
        }
        // OCCT L780-818: grow a shell per face via shared edges.
        let mut a_added: HashSet<(u64, u32)> = HashSet::new();
        let mut shells: Vec<Shape> = Vec::new();
        for a_ff in a_mf {
            if !a_added.insert((a_ff.ptr_id(), a_ff.location)) {
                continue;
            }
            let mut a_shell: Vec<Shape> = Vec::new();
            let mut f0 = a_ff.clone();
            f0.orientation = topods::Orientation::Internal;
            a_shell.push(f0);
            let mut i = 0;
            while i < a_shell.len() {
                let a_f = &a_shell[i];
                for e in face_edges(a_f) {
                    if let Some(a_lf) = a_mef.get(&(e.ptr_id(), e.location)) {
                        for a_fl in a_lf {
                            if a_added.insert((a_fl.ptr_id(), a_fl.location)) {
                                let mut fl = a_fl.clone();
                                fl.orientation = topods::Orientation::Internal;
                                a_shell.push(fl);
                            }
                        }
                    }
                }
                i += 1;
            }
            // OCCT L816: aShell.Closed(BRep_Tool::IsClosed(aShell)).
            let mut flags = tshape_flags::DEFAULT;
            if crate::bop::algo::builder::shell_is_closed(&a_shell) {
                flags |= tshape_flags::CLOSED;
            }
            let shell_tshape = TShape::Shell(TShellData {
                my_shapes: vec![],
                flags,
                faces: a_shell,
            });
            shells.push(Shape::new(
                Arc::new(shell_tshape),
                0,
                topods::Orientation::Forward,
            ));
        }
        shells
    }
}

/// Debug description of a face: surface type + edge count.

/// Extract edge ptr_ids from a Face Shape.
pub(crate) fn face_edge_ptrs(face: &Shape) -> Vec<u64> {    let mut edges = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                for e in &wd.edges {
                    edges.push(e.ptr_id());
                }
            }
            for iw in &fd.inner_wires {
                if let TShape::Wire(wd) = &*iw.data {
                    for e in &wd.edges {
                        edges.push(e.ptr_id());
                    }
                }
            }
        }
        _ => {}
    }
    edges
}

/// Edge endpoint vertex Shapes with composed orientations.
/// OCCT TopoDS_Iterator(aE) with cumOri=true (default) composes the edge
/// orientation into the vertices: each stored vertex keeps its stored
/// orientation composed with the edge's own orientation; a REVERSED edge
/// iterates [last, first] (TopoDS_Iterator.cxx L35-37, L72-80).
fn edge_vertices(e: &Shape) -> Vec<Shape> {
    use rcad_kernel::topods::Orientation;
    let flip_ori = |o: Orientation| -> Orientation {
        match o {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            other => other,
        }
    };
    match &*e.data {
        TShape::Edge(ed) => {
            if e.orientation == Orientation::Reversed {
                vec![
                    Shape::new(ed.last.data.clone(), ed.last.location, flip_ori(ed.last.orientation)),
                    Shape::new(ed.first.data.clone(), ed.first.location, flip_ori(ed.first.orientation)),
                ]
            } else {
                vec![
                    Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
                    Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
                ]
            }
        }
        _ => Vec::new(),
    }
}

/// Extract edge Shapes from a Face Shape (outer + inner wires), composing the
/// face and wire orientations into each edge (OCCT TopExp_Explorer composes
/// the parent orientation at every level: TopExp_Explorer.cxx L152, L110-170).
fn face_edges(face: &Shape) -> Vec<Shape> {    let mut edges = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            // OCCT TopoDS_Iterator(aF) with cumLoc composes the face Location
            // into the edges (TopoDS_Iterator.cxx L76-78): an edge of a located
            // face (e.g. the revolve end cap, the start cap TShape at the
            // rotation Location L1) is keyed by its EFFECTIVE location
            // (face_loc * edge_loc), so the edge→face ancestor map sees the
            // L1-located copies of the same TShape as distinct edges matching
            // the located boundary edges of the neighboring faces.
            let face_loc = face.location;
            let compose = |eloc: u32| -> u32 {
                if face_loc == 0 {
                    eloc
                } else if eloc == 0 {
                    face_loc
                } else {
                    // Both non-identity: the composed transform is one of the
                    // two operands in the single-fold cases (the ring's end cap
                    // at L1 with identity edges); fall back to the face
                    // location index.
                    face_loc
                }
            };
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                let w_or = fd.outer_wire.orientation;
                for e in &wd.edges {
                    let mut e2 = e.clone();
                    e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                    e2.location = compose(e2.location);
                    edges.push(e2);
                }
            }
            for iw in &fd.inner_wires {
                if let TShape::Wire(wd) = &*iw.data {
                    let w_or = iw.orientation;
                    for e in &wd.edges {
                        let mut e2 = e.clone();
                        e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                        e2.location = compose(e2.location);
                        edges.push(e2);
                    }
                }
            }
        }
        _ => {}
    }
    edges
}

/// OCCT BRep_Tool::IsClosed(aE, aF) — the edge has two pcurves on the closed
/// OCCT BRep_Tool::IsClosed(aE, aF) (BRep_Tool.cxx L795-841) — the edge has
/// two pcurves on the closed surface of the face (seam edge). OCCT matches
/// the CurveOnClosedSurface representation by the face's SURFACE handle and
/// returns false for plane faces outright (IsPlane short-circuit, L819-822).
/// rcad keys the representation by the face's BRep index: try the exact index
/// match first (original faces), then fall back to the surface-level match for
/// the BuilderSolid split-image faces whose index does not preserve the
/// original face's (BOPAlgo_BuilderSolid.cxx L196).
pub(crate) fn edge_closed_on_face(a_e: &Shape, a_f: &Shape) -> bool {
    let f_key = (a_f.ptr_id(), a_f.location);
    let ed = match a_e.as_edge() {
        Some(ed) => ed,
        None => return false,
    };
    // OCCT L825-840: any CurveOnClosedSurface on the edge whose surface
    // matches the face's surface — exact index match for original faces.
    if ed.representations.iter().any(|r| {
        matches!(
            r,
            topods::CurveRepresentation::CurveOnClosedSurface { face, .. } if *face == f_key
        )
    }) {
        return true;
    }
    // OCCT L819-822: if (IsPlane(S)) return false.
    match &*a_f.data {
        TShape::Face(fd) => {
            if matches!(fd.surface, Some(rcad_kernel::geom::Surface3::Plane(_))) {
                return false;
            }
        }
        _ => return false,
    }
    // Fallback: split-face images share the original face's closed surface;
    // any CurveOnClosedSurface on the edge matches that surface.
    ed.representations
        .iter()
        .any(|r| matches!(r, topods::CurveRepresentation::CurveOnClosedSurface { .. }))
}

