// OCCT BOPAlgo_BuilderSolid — solid building from shells.
//
// OCCT BOPAlgo_BuilderSolid.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::{Alert, Report};
use crate::bop::algo::shell_splitter::{make_connexity_blocks, make_shells};
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use rcad_kernel::geom::CurveEval;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, ShapeType, TShape, TSolidData, TShellData, TFaceData, TWireData, tshape_flags};
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
    my_shapes_to_avoid: HashSet<u64>,   // OCCT: myShapesToAvoid
    my_loops: Vec<Vec<Shape>>,          // OCCT: myLoops
    my_loops_internal: Vec<Vec<Shape>>, // OCCT: myLoopsInternal
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
            my_shapes_to_avoid: HashSet::new(),
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
        // OCCT L142-218: iterative — mark faces with free edges, repeat.
        loop {
            let mut b_found = false;
            // OCCT L151-160: edge -> [faces] for non-avoided faces.
            // key: edge ptr_id; value: (faces using the edge, first edge shape).
            let mut a_mef: HashMap<u64, (Vec<Shape>, Shape)> = HashMap::new();
            for face in &self.my_shapes {
                if self.my_shapes_to_avoid.contains(&face.ptr_id()) { continue; }
                for e in face_edges(face) {
                    // OCCT L167-169: skip degenerated edges.
                    let degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                    if degen { continue; }
                    let eptr = e.ptr_id();
                    let entry = a_mef.entry(eptr).or_insert_with(|| (Vec::new(), e.clone()));
                    entry.0.push(face.clone());
                }
            }
            // OCCT L163-211: per-edge decisions (iterate in sorted order for
            // determinism; the hasher ignores orientation like TopTools_ShapeMapHasher).
            let mut a_nb_e: Vec<u64> = a_mef.keys().copied().collect();
            a_nb_e.sort();
            for eptr in a_nb_e {
                let (a_lf, a_e) = &a_mef[&eptr];
                let a_nb_f = a_lf.len();
                if a_nb_f == 0 { continue; }
                // OCCT L179: aOrE = aE.Orientation().
                let a_or_e = a_e.orientation;
                let a_f1 = &a_lf[0];
                if a_nb_f == 1 {
                    // OCCT L182-189: aNbF==1, skip INTERNAL edges.
                    if a_or_e == topods::Orientation::Internal { continue; }
                    b_found = true;
                    self.my_shapes_to_avoid.insert(a_f1.ptr_id());
                } else if a_nb_f == 2 {
                    // OCCT L191-209: aNbF==2, same face (IsSame) && !IsClosed
                    // && edge != INTERNAL -> avoid both copies.
                    let a_f2 = &a_lf[1];
                    if a_f2.ptr_id() == a_f1.ptr_id() {
                        if edge_closed_on_face(a_e, a_f1) { continue; }
                        if a_or_e == topods::Orientation::Internal { continue; }
                        b_found = true;
                        self.my_shapes_to_avoid.insert(a_f1.ptr_id());
                        self.my_shapes_to_avoid.insert(a_f2.ptr_id());
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
            if !self.my_shapes_to_avoid.contains(&face.ptr_id()) {
                a_start.push(face.clone());
            }
        }
        // OCCT L259-284: ShellSplitter.Perform -> shells appended to myLoops.
        let blocks = make_connexity_blocks(&a_start);
        let shells = make_shells(&blocks);
        for shell in shells {
            if !shell.is_empty() {
                self.my_loops.push(shell);
            }
        }
        // OCCT L287-331: post-treatment — collect all faces of the loops and
        // of myShapesToAvoid; add the remaining faces to myShapesToAvoid.
        let mut a_mp: HashSet<u64> = HashSet::new();
        for loop_faces in &self.my_loops {
            for f in loop_faces {
                a_mp.insert(f.ptr_id());
            }
        }
        for &fptr in &self.my_shapes_to_avoid {
            a_mp.insert(fptr);
        }
        for face in &self.my_shapes {
            if !Self::is_infinite_face(face) {
                if !a_mp.contains(&face.ptr_id()) {
                    self.my_shapes_to_avoid.insert(face.ptr_id());
                }
            }
        }
        // OCCT L338-392: build internal shells from the avoided faces.
        // OCCT L339: myLoopsInternal.Clear().
        self.my_loops_internal.clear();
        // OCCT L344-349: edge -> [faces] map from the avoided faces.
        let mut a_ef_map: HashMap<u64, Vec<Shape>> = HashMap::new();
        for face in &self.my_shapes {
            if !self.my_shapes_to_avoid.contains(&face.ptr_id()) { continue; }
            for e in face_edges(face) {
                a_ef_map.entry(e.ptr_id()).or_default().push(face.clone());
            }
        }
        // OCCT L351-391: grow a shell per avoided face via shared edges.
        let mut a_added: HashSet<u64> = HashSet::new();
        for a_ff in &self.my_shapes {
            if !self.my_shapes_to_avoid.contains(&a_ff.ptr_id()) { continue; }
            if !a_added.insert(a_ff.ptr_id()) { continue; }
            let mut a_shell: Vec<Shape> = vec![a_ff.clone()];
            let mut i = 0;
            while i < a_shell.len() {
                let a_f = &a_shell[i];
                for e in face_edges(a_f) {
                    if let Some(a_lf) = a_ef_map.get(&e.ptr_id()) {
                        for a_fl in a_lf {
                            if a_added.insert(a_fl.ptr_id()) {
                                a_shell.push(a_fl.clone());
                            }
                        }
                    }
                }
                i += 1;
            }
            self.my_loops_internal.push(a_shell);
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
        let mut a_mhf: HashSet<u64> = HashSet::new();

        // OCCT L413-442: classify each shell.
        for loop_faces in &self.my_loops {
            // OCCT L422-427: IsGrowthShell then IsHole.
            let mut b_is_growth = Self::is_growth_shell(loop_faces, &a_mhf);
            if !b_is_growth {
                b_is_growth = !Self::is_hole(loop_faces);
            }
            if b_is_growth {
                // OCCT L430-435: MakeSolid + Add(aSolid, aShell).
                let shell_shape = self.build_shell_shape(loop_faces);
                let solid = self.build_solid_shape(&shell_shape);
                new_solids.push(solid);
            } else {
                // OCCT L437-441: add to hole shells + MapShapes(aShell, FACE, aMHF).
                for f in loop_faces {
                    a_mhf.insert(f.ptr_id());
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
        let mut a_hole_solid: HashMap<Vec<u64>, usize> = HashMap::new();
        for (si, solid) in new_solids.iter().enumerate() {
            for hs in &hole_shells {
                // OCCT L511: if (!IsInside(aHole, aSolid)) continue;
                if !Self::is_inside(hs, solid) { continue; }
                let hkey = Self::hole_shell_key(hs);
                match a_hole_solid.get(&hkey) {
                    None => {
                        // OCCT L526-528: aHoleSolidMap.Add(aHole, aSolid).
                        a_hole_solid.insert(hkey, si);
                    }
                    Some(&prev) => {
                        // OCCT L520-523: if (IsInside(aSolid, *pSolidWas)) *pSolidWas = aSolid.
                        if Self::is_inside_solid(&new_solids[si], &new_solids[prev]) {
                            a_hole_solid.insert(hkey, si);
                        }
                    }
                }
            }
        }

        // OCCT L532-548: back map solids -> holes.
        let mut a_solid_holes: HashMap<usize, Vec<usize>> = HashMap::new();
        for (hi, hs) in hole_shells.iter().enumerate() {
            let hkey = Self::hole_shell_key(hs);
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
            let hkey = Self::hole_shell_key(hs);
            if !a_hole_solid.contains_key(&hkey) {
                let shell_shape = self.build_shell_shape(hs);
                let solid = self.build_solid_shape(&shell_shape);
                self.my_solids.push(solid);
            }
        }
    }

    /// OCCT IsGrowthShell (BuilderSolid.cxx L864-879) — the shell contains a
    /// previously-marked hole face, so it bounds a hole from the outside.
    fn is_growth_shell(faces: &[Shape], the_mhf: &HashSet<u64>) -> bool {
        if the_mhf.is_empty() {
            return false;
        }
        faces.iter().any(|f| the_mhf.contains(&f.ptr_id()))
    }

    /// OCCT IsHole (BuilderSolid.cxx L823-831) — the shell is a "hole in
    /// space": the point at infinity is IN the shell-solid (material outside).
    fn is_hole(faces: &[Shape]) -> bool {
        Self::is_hole_in_space(faces)
    }

    /// OCCT IsInside (BuilderSolid.cxx L835-860) — classify the first face of
    /// the shell relative to the solid via BOPTools_AlgoTools::ComputeState.
    fn is_inside(faces: &[Shape], solid: &Shape) -> bool {
        let Some(a_f) = faces.first() else {
            return false;
        };
        Self::compute_state_on_solid(a_f, solid) == 3 // TopAbs_IN
    }

    /// IsInside(aSolid, bSolid) — classify the first face of solid a against b.
    fn is_inside_solid(a: &Shape, b: &Shape) -> bool {
        let Some(a_f) = Self::first_face(a) else {
            return false;
        };
        Self::compute_state_on_solid(&a_f, b) == 3 // TopAbs_IN
    }

    /// Classify a point relative to a solid (OCCT BRepClass3d_SolidClassifier).
    fn point_classify(solid: &Shape, p: DVec3) -> u8 {
        let mut clsf =
            crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(solid);
        clsf.perform(p, 1e-7);
        clsf.my_state
    }

    /// OCCT BOPTools_AlgoTools::ComputeState(aF, aSolid, aTol, aBounds)
    /// (BOPTools_AlgoTools.cxx L660-715): try an edge of the face not on the
    /// solid (classify the edge midpoint), else classify a point inside the face.
    fn compute_state_on_solid(a_f: &Shape, solid: &Shape) -> u8 {
        let solid_edges = Self::solid_edges(solid);
        for e in face_edges(a_f) {
            if e.as_edge().map_or(true, |ed| ed.degenerated) {
                continue;
            }
            if !solid_edges.contains(&e.ptr_id()) {
                let p = Self::edge_midpoint(&e);
                return Self::point_classify(solid, p);
            }
        }
        let p = Self::face_centroid(a_f);
        Self::point_classify(solid, p)
    }

    /// All edge ptr_ids of a solid (OCCT TopExp::MapShapes(aSolid, TopAbs_EDGE)).
    fn solid_edges(solid: &Shape) -> HashSet<u64> {
        let mut set: HashSet<u64> = HashSet::new();
        let mut stack: Vec<Shape> = vec![solid.clone()];
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
                TShape::Face(_) => {
                    for e in face_edges(&sh) {
                        set.insert(e.ptr_id());
                    }
                }
                _ => {}
            }
        }
        set
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

    /// Middle point of an edge (OCCT IntTools_Tools::IntermediatePoint on the
    /// 3D curve).
    fn edge_midpoint(edge: &Shape) -> DVec3 {
        if let TShape::Edge(ed) = &*edge.data {
            if let Some(curve) = &ed.curve {
                let t = crate::bop::int_tools::face_make_curve::intermediate_point(
                    ed.range[0],
                    ed.range[1],
                );
                return curve.point_at(t);
            }
            if let TShape::Vertex(vd) = &*ed.first.data {
                return vd.point;
            }
        }
        DVec3::ZERO
    }

    /// Canonical key of a shell face set (OCCT TopoDS_Shape map key equivalent).
    fn hole_shell_key(faces: &[Shape]) -> Vec<u64> {
        let mut v: Vec<u64> = faces.iter().map(|f| f.ptr_id()).collect();
        v.sort();
        v
    }

    /// OCCT BOPAlgo_BuilderSolid::IsHole (BuilderSolid.cxx L823-831) via
    /// BRepClass3d_SClassifier::PerformInfinitePoint (SClassifier.cxx L82-199).
    ///
    /// Casts a ray from a face's inner point along the reversed oriented normal
    /// (OCCT L141: aLin(aPoint, -aDN)) and finds the closest intersection with
    /// the shell's faces (OCCT L143-180). The ray starts on the probe face
    /// itself at w=0; its transition is In (the ray is opposite to the face
    /// normal). A closest-intersection transition of Out means the ray exits
    /// the material at the closest face, so the point at infinity is IN → hole
    /// (OCCT L184-189); a transition of In means the infinite point is OUT
    /// (OCCT L192-195).
    fn is_hole_in_space(faces: &[Shape]) -> bool {
        if faces.is_empty() {
            return false;
        }
        // OCCT L125-198: try each face as the probe (up to 10 random points per
        // face; rcad uses the face centroid). The first definitive answer wins.
        for f in faces {
            let p = Self::face_centroid(f);
            let n = match Self::face_outward_normal(f) {
                Some(n) => n,
                None => continue,
            };
            if n.length_squared() < 1e-12 {
                continue;
            }
            // OCCT L141: the ray direction is the reversed normal.
            let ray_dir = -n;
            // Find the closest intersection of the ray with all the faces
            // (OCCT L143-180, the minimal WParameter wins).
            let mut parmin = f64::MAX;
            let mut best = 0u8; // 0=no valid transition, 1=In, 2=Out
            for g in faces {
                let g_n = match Self::face_outward_normal(g) {
                    Some(v) => v,
                    None => continue,
                };
                let g_o = match Self::face_plane_origin(g) {
                    Some(v) => v,
                    None => continue,
                };
                let denom = ray_dir.dot(g_n);
                if denom.abs() < 1e-12 {
                    continue;
                }
                let w = (g_o - p).dot(g_n) / denom;
                // OCCT L168-171: parmin = WParameter(imin) is assigned for the
                // minimal parameter of every face unconditionally, so a later
                // face with an equal parameter (w <= parmin) overwrites the
                // transition. The probe face itself always sits at w=0 (the ray
                // starts on its surface); a second copy of the same face with the
                // reversed orientation also lands at w=0 and its Out transition
                // then wins, which is how a degenerate double-face shell
                // (a face in both orientations) is detected as a hole.
                if w <= parmin {
                    parmin = w;
                    // OCCT int_cs transition: cos_dir = nSurf · dirCurve; <0 -> In,
                    // >0 -> Out (IntCurveSurface_InterUtils.pxx L856-895).
                    best = if denom > 0.0 { 2 } else { 1 };
                }
            }
            if best == 2 {
                // OCCT L184-189: transition Out -> the infinite point is IN.
                return true;
            } else if best == 1 {
                // OCCT L192-195: transition In -> the infinite point is OUT.
                return false;
            }
        }
        false
    }

    /// Outward normal of a face: the surface normal flipped by the face
    /// orientation (OCCT BRepClass3d_SClassifier::FaceNormal L606-627).
    fn face_outward_normal(f: &Shape) -> Option<DVec3> {
        match &*f.data {
            TShape::Face(fd) => {
                let surf = fd.surface.as_ref()?;
                let n = match surf {
                    rcad_kernel::geom::Surface3::Plane(pl) => pl.normal,
                    _ => return None, // curved surfaces need point evaluation
                };
                Some(if f.orientation == rcad_kernel::topods::Orientation::Reversed {
                    -n
                } else {
                    n
                })
            }
            _ => None,
        }
    }

    /// Origin point of a planar face's surface.
    fn face_plane_origin(f: &Shape) -> Option<DVec3> {
        match &*f.data {
            TShape::Face(fd) => match fd.surface.as_ref()? {
                rcad_kernel::geom::Surface3::Plane(pl) => Some(pl.origin),
                _ => None,
            },
            _ => None,
        }
    }

    /// Compute face centroid from its outer wire vertices.
    fn face_centroid(face: &Shape) -> DVec3 {
        match &*face.data {
            TShape::Face(fd) => {
                let mut pts: Vec<DVec3> = Vec::new();
                if let TShape::Wire(wd) = &*fd.outer_wire.data {
                    for e in &wd.edges {
                        if let TShape::Edge(ed) = &*e.data {
                            if let TShape::Vertex(vd) = &*ed.first.data { pts.push(vd.point); }
                            if let TShape::Vertex(vd) = &*ed.last.data { pts.push(vd.point); }
                        }
                    }
                }
                if pts.is_empty() { DVec3::ZERO } else { pts.iter().sum::<DVec3>() / pts.len() as f64 }
            }
            _ => DVec3::ZERO,
        }
    }

    /// Build a Shell TShape from a set of faces.
    fn build_shell_shape(&self, faces: &[Shape]) -> Shape {
        let shell_tshape = TShape::Shell(TShellData {
            my_shapes: vec![],
            flags: tshape_flags::DEFAULT,
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
        // OCCT L604-608: myAvoidInternalShapes — rcad has no such option.
        // OCCT L610-614: no internal parts -> return.
        if self.my_loops_internal.is_empty() {
            return;
        }
        // OCCT L619-629: collect all faces of the internal shells into aMFs.
        let mut a_mfs: Vec<Shape> = Vec::new();
        let mut a_mfs_fence: HashSet<u64> = HashSet::new();
        for shell in &self.my_loops_internal {
            for f in shell {
                if a_mfs_fence.insert(f.ptr_id()) {
                    a_mfs.push(f.clone());
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
        // OCCT L661-681: classify the internal faces relative to the areas
        // (BOPAlgo_Tools::ClassifyFaces — semantic point-in-solid per face).
        let mut a_mslf: HashMap<usize, Vec<Shape>> = HashMap::new();
        for (si, solid) in self.my_solids.iter().enumerate() {
            for a_f in &a_mfs {
                let pt = Self::face_centroid(a_f);
                if Self::point_classify(solid, pt) == 3 {
                    a_mslf.entry(si).or_default().push(a_f.clone());
                }
            }
        }
        // OCCT L685-722: update the solids by their internal faces.
        let mut a_mf_done: HashSet<u64> = HashSet::new();
        for (si, a_lf) in a_mslf {
            if a_lf.is_empty() {
                continue;
            }
            for a_f in &a_lf {
                a_mf_done.insert(a_f.ptr_id());
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
            if !a_mf_done.contains(&a_f.ptr_id()) {
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
        // OCCT L773-778: edge -> [faces] map.
        let mut a_mef: HashMap<u64, Vec<Shape>> = HashMap::new();
        for a_f in a_mf {
            for e in face_edges(a_f) {
                a_mef.entry(e.ptr_id()).or_default().push(a_f.clone());
            }
        }
        // OCCT L780-818: grow a shell per face via shared edges.
        let mut a_added: HashSet<u64> = HashSet::new();
        let mut shells: Vec<Shape> = Vec::new();
        for a_ff in a_mf {
            if !a_added.insert(a_ff.ptr_id()) {
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
                    if let Some(a_lf) = a_mef.get(&e.ptr_id()) {
                        for a_fl in a_lf {
                            if a_added.insert(a_fl.ptr_id()) {
                                let mut fl = a_fl.clone();
                                fl.orientation = topods::Orientation::Internal;
                                a_shell.push(fl);
                            }
                        }
                    }
                }
                i += 1;
            }
            let shell_tshape = TShape::Shell(TShellData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
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

/// Extract edge ptr_ids from a Face Shape.
pub(crate) fn face_edge_ptrs(face: &Shape) -> Vec<u64> {
    let mut edges = Vec::new();
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

/// Extract edge Shapes from a Face Shape (outer + inner wires), composing the
/// face and wire orientations into each edge (OCCT TopExp_Explorer composes
/// the parent orientation at every level: TopExp_Explorer.cxx L152, L110-170).
fn face_edges(face: &Shape) -> Vec<Shape> {
    let mut edges = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                let w_or = fd.outer_wire.orientation;
                for e in &wd.edges {
                    let mut e2 = e.clone();
                    e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                    edges.push(e2);
                }
            }
            for iw in &fd.inner_wires {
                if let TShape::Wire(wd) = &*iw.data {
                    let w_or = iw.orientation;
                    for e in &wd.edges {
                        let mut e2 = e.clone();
                        e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
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
/// surface (seam edge). BOPAlgo_Builder_2.cxx L397.
fn edge_closed_on_face(a_e: &Shape, a_f: &Shape) -> bool {
    let f_index = a_f.index;
    a_e.as_edge()
        .map(|ed| {
            ed.representations.iter().any(|r| {
                matches!(
                    r,
                    topods::CurveRepresentation::CurveOnClosedSurface { face, .. }
                        if *face == f_index
                )
            })
        })
        .unwrap_or(false)
}
