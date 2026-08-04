// OCCT BOPAlgo_BuilderSolid — solid building from shells.
//
// OCCT BOPAlgo_BuilderSolid.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use crate::bop::algo::shell_splitter::make_connexity_blocks_from_shapes;
use crate::bop::int_tools::context::IntToolsContext;
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
    pub my_solids: Vec<Shape>,          // OCCT: mySolids
    my_shapes_to_avoid: HashSet<u64>,   // OCCT: myShapesToAvoid
    my_loops: Vec<Vec<Shape>>,          // OCCT: myLoops
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
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::Perform (BOPAlgo_BuilderSolid.cxx L76-125).
    pub fn perform(&mut self) {
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

    /// OCCT BOPAlgo_BuilderSolid::PerformShapesToAvoid (BuilderSolid.cxx L129-220).
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT L138: myShapesToAvoid.Clear()
        self.my_shapes_to_avoid.clear();
        // OCCT L142-218: iterative — remove faces with free edges, repeat
        loop {
            // OCCT L151-158: build edge→[faces] map for non-avoided faces
            let mut a_mef: HashMap<u64, Vec<usize>> = HashMap::new();
            for (fi, face) in self.my_shapes.iter().enumerate() {
                if self.my_shapes_to_avoid.contains(&face.ptr_id()) { continue; }
                for eptr in face_edge_ptrs(face) {
                    a_mef.entry(eptr).or_default().push(fi);
                }
            }
            // OCCT L162-180: find faces with edges used only once
            let mut b_found = false;
            for face in &self.my_shapes {
                if self.my_shapes_to_avoid.contains(&face.ptr_id()) { continue; }
                let has_free = face_edge_ptrs(face).iter().any(|eptr| {
                    a_mef.get(eptr).map_or(true, |v| v.len() == 1)
                });
                if has_free {
                    self.my_shapes_to_avoid.insert(face.ptr_id());
                    b_found = true;
                }
            }
            if !b_found { break; }
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformLoops (L237-393).
    /// Groups faces into connected shells using edge adjacency.
    fn perform_loops(&mut self) {
        if self.my_shapes.is_empty() {
            return;
        }
        // OCCT L237-276: build shells via ShellSplitter (MakeConnexityBlocks)
        let mut blocks: Vec<Vec<usize>> = Vec::new();
        make_connexity_blocks_from_shapes(&self.my_shapes, self.ds, &mut blocks);
        self.my_loops.clear();
        for block in &blocks {
            let faces: Vec<Shape> = block.iter()
                .map(|&idx| self.my_shapes[idx].clone())
                .collect();
            if !faces.is_empty() {
                self.my_loops.push(faces);
            }
        }
        // OCCT L285-331: post-treatment — find unprocessed faces
        let mut processed: HashSet<u64> = HashSet::new();
        for loop_faces in &self.my_loops {
            for f in loop_faces {
                processed.insert(f.ptr_id());
            }
        }
        // OCCT L331: add unprocessed faces (with free edges) to myShapesToAvoid
        for face in &self.my_shapes {
            if !processed.contains(&face.ptr_id()) {
                self.my_shapes_to_avoid.insert(face.ptr_id());
            }
        }
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformAreas (L397-530).
    /// Classifies shells as growth (solid) or hole (void).
    fn perform_areas(&mut self) {
        // OCCT L400-404: growth solids and hole shells
        let mut new_solids: Vec<Shape> = Vec::new();
        let mut hole_shells: Vec<Vec<Shape>> = Vec::new();
        let mut hole_face_ptrs: HashSet<u64> = HashSet::new();

        for loop_faces in &self.my_loops {
            // OCCT L422: check if growth shell
            let is_growth = !self.is_hole_shell(loop_faces, &hole_face_ptrs);
            if is_growth {
                // OCCT L431-435: build Solid from Shell
                let shell_shape = self.build_shell_shape(loop_faces);
                let solid = self.build_solid_shape(&shell_shape);
                new_solids.push(solid);
            } else {
                // OCCT L438-440: add to hole shells, track face ptrs
                for f in loop_faces {
                    hole_face_ptrs.insert(f.ptr_id());
                }
                hole_shells.push(loop_faces.clone());
            }
        }

        // OCCT L444-457: no holes — all growths are the result
        if hole_shells.is_empty() {
            self.my_solids = new_solids;
            return;
        }

        // OCCT L460+: classify holes relative to growth solids
        // Simplified: each hole becomes internal to the nearest growth
        if !new_solids.is_empty() {
            for hs in &hole_shells {
                let shell_shape = self.build_shell_shape(hs);
                // Add hole shell as internal of first growth solid
                if let Some(first) = new_solids.first_mut() {
                    let ts = Arc::make_mut(&mut first.data);
                    if let TShape::Solid(sd) = ts {
                        sd.internal_vertices.push(shell_shape);
                    }
                }
            }
        }
        self.my_solids = new_solids;
    }

    /// OCCT BOPAlgo_BuilderSolid::PerformAreas (BuilderSolid.cxx L422-427).
    /// A shell is a hole when: (1) fast check IsGrowthShell fails — the shell
    /// does not contain a previously-found hole face; and (2) IsHole succeeds —
    /// the shell is a "hole in space" (its faces are oriented with the material
    /// outside, so the point at infinity is IN the shell-solid).
    fn is_hole_shell(&self, faces: &[Shape], hole_face_ptrs: &HashSet<u64>) -> bool {
        // OCCT IsGrowthShell (L864-879): if the shell contains any hole-face
        // marker it bounds a hole from the outside → it is a growth, not a hole.
        if faces.iter().any(|f| hole_face_ptrs.contains(&f.ptr_id())) {
            return false;
        }
        // OCCT IsHole (L823-831): the shell is a hole-in-space (inverted).
        self.is_hole_in_space(faces)
    }

    /// OCCT BOPAlgo_BuilderSolid::IsHole (BuilderSolid.cxx L823-831) via
    /// BRepClass3d_SClassifier::PerformInfinitePoint (L82-199). From a point on
    /// a face, cast a ray along the reversed outward normal into the material;
    /// at the closest crossing, an EXIT means the material is inside (growth),
    /// an ENTER means the material is outside (hole-in-space).
    fn is_hole_in_space(&self, faces: &[Shape]) -> bool {
        if faces.is_empty() {
            return false;
        }
        for f in faces {
            let p = Self::face_centroid(f);
            let n = match Self::face_outward_normal(f) {
                Some(n) => n,
                None => continue,
            };
            let d = -n;
            let mut best_t = f64::MAX;
            let mut best_denom = 0.0f64;
            for f2 in faces {
                if f2.ptr_id() == f.ptr_id() {
                    continue;
                }
                let n2 = match Self::face_outward_normal(f2) {
                    Some(n2) => n2,
                    None => continue,
                };
                let denom = d.dot(n2);
                if denom.abs() < 1e-12 {
                    continue;
                }
                // t for plane (x - origin)·n2 = 0 along ray p + t·d.
                let t = match Self::face_plane_origin(f2) {
                    Some(o2) => (o2 - p).dot(n2) / denom,
                    None => continue,
                };
                if t > 1e-7 && t < best_t {
                    best_t = t;
                    best_denom = denom;
                }
            }
            if best_t < f64::MAX {
                // denom > 0 → ray moves toward +n2 (outside) → EXIT → growth.
                // denom < 0 → ray moves toward -n2 (material) → ENTER → hole.
                return best_denom < 0.0;
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

    /// OCCT BOPAlgo_BuilderSolid::PerformInternalShapes.
    fn perform_internal_shapes(&mut self) {
        // OCCT BOPAlgo_BuilderSolid::PerformInternalShapes (BuilderSolid.cxx L602-660).
        // rcad: internal shapes use solid_classifier_is_above.
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
