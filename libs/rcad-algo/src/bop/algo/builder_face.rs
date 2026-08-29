// OCCT BOPAlgo_BuilderFace — face splitting with section edges.
//
// OCCT BOPAlgo_BuilderFace.cxx
// Performs: PerformShapesToAvoid -> PerformLoops -> PerformAreas -> PerformInternalShapes

use crate::bop::algo::Report;
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use glam::DVec2;
use indexmap::IndexMap;
use rcad_kernel::geom::{Curve2d, Curve2dEval, SurfaceEval};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{TShape, TWireData, tshape_flags};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// OCCT BOPAlgo_BuilderFace — splits a face using section edges.
pub struct BuilderFace<'a> {
    ds: &'a DS,
    // BOPAlgo_Algo (inherited)
    my_report: Report,
    my_run_parallel: bool,
    // BOPAlgo_BuilderFace
    pub my_face: Option<Shape>,         // OCCT: myFace
    pub my_face_index: Option<usize>,   // rcad: DS index for my_face
    pub my_edges: Vec<Shape>,           // OCCT: myShapes (section edges)
    pub my_areas: Vec<Shape>,           // OCCT: myAreas (result faces)
    pub my_loops: Vec<Shape>,           // OCCT: myLoops (result wires)
    pub my_loops_internal: Vec<Shape>,  // OCCT: myLoopsInternal (internal wires)
    my_shapes_to_avoid: HashSet<(u64, u32, rcad_kernel::topods::Orientation)>, // OCCT: myShapesToAvoid (NCollection_Map — default hasher TShape+Location+Orientation)
    my_avoid_internal_shapes: bool,      // OCCT: myAvoidInternalShapes (BuilderArea)
    my_context: IntToolsContext,         // OCCT: myContext
}

impl<'a> BuilderFace<'a> {
    pub fn new(ds: &'a DS) -> Self {
        BuilderFace {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_face: None,
            my_face_index: None,
            my_edges: Vec::new(),
            my_areas: Vec::new(),
            my_loops: Vec::new(),
            my_loops_internal: Vec::new(),
            my_shapes_to_avoid: HashSet::new(),
            my_avoid_internal_shapes: false,
            my_context: IntToolsContext::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderFace::Perform (BOPAlgo_BuilderFace.cxx L117-148).
    pub fn perform(&mut self) {
        // OCCT L121: GetReport()->Clear().
        self.my_report.clear();
        if self.my_face.is_none() {
            return;
        }
        // OCCT L124: PerformShapesToAvoid
        self.perform_shapes_to_avoid();
        if self.has_errors() { return; }
        // OCCT L130: PerformLoops — build closed wires from edges
        self.perform_loops();
        if self.has_errors() { return; }
        // OCCT L136: PerformAreas — classify areas as IN/OUT
        self.perform_areas();
        if self.has_errors() { return; }
        // OCCT L147: PerformInternalShapes
        self.perform_internal_shapes();
        // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L327-373) matches the
        // edge pcurve by the SURFACE geometry, so the split faces (areas) —
        // which share the source surface — resolve their edge pcurves written
        // under the source face key.  rcad's CurveOnSurface matches by face
        // TShape pointer, so each area face needs its own pcurve entry: copy
        // the source face key's pcurve to the area face key (the shared
        // surface makes the pcurve identical).
        if let Some(src) = self.my_face.as_ref() {
            for area in &self.my_areas {
                let area_key = (area.ptr_id(), area.location);
                let mut es: Vec<Shape> = Vec::new();
                if let TShape::Face(fd) = &*area.data {
                    if let TShape::Wire(wd) = &*fd.outer_wire.data {
                        es.extend(wd.edges.iter().cloned());
                    }
                    for w in &fd.inner_wires {
                        if let TShape::Wire(wd) = &*w.data {
                            es.extend(wd.edges.iter().cloned());
                        }
                    }
                }
                for e in es {
                    // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the
                    // pcurve key location is L.Predivided(E.Location()) — the
                    // face location divided by the edge's location. A LOCATED
                    // split edge (e.g. the x=1 cap edges of a prism at the
                    // translation location) stores its pcurve under the
                    // composed key; the raw face key would miss it and the
                    // area face's boundary (Green) integral would fall back to
                    // a zero area.
                    let src_key = (
                        src.ptr_id(),
                        crate::bop::algo::compose_face_edge_pcurve_location(src.location, e.location, &self.ds.locations),
                    );
                    let raw = Arc::as_ptr(&e.data) as *mut TShape;
                    unsafe {
                        if let TShape::Edge(ed) = &mut *raw {
                            let v = ed.pcurves.get(&src_key).cloned();
                            if let Some(v) = v {
                                ed.pcurves.entry(area_key).or_insert(v);
                            }
                            // OCCT BRep_Tool::CurveOnSurface matches by the
                            // surface geometry (BRep_Tool.cxx L350), so the
                            // seam edge's CurveOnClosedSurface representation
                            // resolves on the area faces too (they share the
                            // source surface).  rcad keys by (face TShape,
                            // location), so the representation is copied to
                            // the area face key as well.
                            let reps = ed.representations.clone();
                            for r in reps {
                                let r2 = match r {
                                    rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face, pcurve, range } if face == src_key =>
                                        Some(rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face: area_key, pcurve, range }),
                                    rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, pcurve1, pcurve2, range } if face == src_key =>
                                        Some(rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face: area_key, pcurve1, pcurve2, range }),
                                    _ => None,
                                };
                                if let Some(r2) = r2 {
                                    use rcad_kernel::topods::CurveRepresentation;
                                    let dup = ed.representations.iter().any(|rr| match (&rr, &r2) {
                                        (CurveRepresentation::CurveOnClosedSurface { face: f1, .. }, CurveRepresentation::CurveOnClosedSurface { face: f2, .. }) => f1 == f2,
                                        (CurveRepresentation::CurveOnSurface { face: f1, .. }, CurveRepresentation::CurveOnSurface { face: f2, .. }) => f1 == f2,
                                        _ => false,
                                    });
                                    if !dup {
                                        ed.representations.push(r2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn has_errors(&self) -> bool { self.my_report.has_errors() }

    /// OCCT BOPAlgo_BuilderFace::GetReport — the accumulated alerts.
    pub fn report(&self) -> &Report { &self.my_report }

    /// OCCT BOPAlgo_BuilderFace::PerformShapesToAvoid (BuilderFace.cxx L152-235).
    /// Iteratively marks edges with a free boundary (a vertex used by at most
    /// one non-avoided edge, or by two IsSame copies of one edge) as "to avoid".
    fn perform_shapes_to_avoid(&mut self) {
        // OCCT L160: myShapesToAvoid.Clear()
        self.my_shapes_to_avoid.clear();
        // OCCT L164-234: iterate until no more edges are found.
        loop {
            let mut b_found = false;
            // OCCT L173-182: aMVE — vertex → [edges] (skipping avoided edges).
            // OCCT aMVE is NCollection_IndexedDataMap with TopTools_ShapeMapHasher
            // (L156-157) — insertion order, key identity TShape + Location.
            let mut a_mve: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            for a_e in &self.my_edges {
                if self.my_shapes_to_avoid.contains(&(a_e.ptr_id(), a_e.location, a_e.orientation)) { continue; }
                for a_v in Self::edge_vertices(a_e, &self.ds.locations) {
                    if std::env::var("RCAD_AVOID_DEBUG").is_ok() {
                        let p = a_v.as_vertex().map(|vd| {
                            let m = self.ds.locations.get(
                                crate::bop::algo::compose_edge_vertex_location(a_e.location, a_v.location, &self.ds.locations) as usize,
                            ).copied().unwrap_or(glam::DAffine3::IDENTITY);
                            let w = m.transform_point3(vd.point);
                            format!("({:.2},{:.2},{:.2})", w.x, w.y, w.z)
                        }).unwrap_or_default();
                        eprintln!(
                            "[AVOID] face={:?} e={:+x}@{} v={:+x}@{} wpt={} -> mve",
                            self.my_face_index.map(|i| i), a_e.ptr_id() & 0xffff, a_e.location, a_v.ptr_id() & 0xffff, a_v.location, p);
                    }
                    let entry = a_mve
                        .entry((a_v.ptr_id(), a_v.location))
                        .or_insert_with(|| (a_v.clone(), Vec::new()));
                    entry.1.push(a_e.clone());
                }
            }
            // OCCT L186-228: for each vertex decide.
            for ((_vptr, _vloc), (a_v, a_le)) in &a_mve {
                let a_nb_e = a_le.len();
                if a_nb_e == 0 { continue; }
                let a_e1 = &a_le[0];
                if a_nb_e == 1 {
                    // OCCT L198-210: single edge at the vertex.
                    if a_e1.as_edge().map_or(true, |ed| ed.degenerated) {
                        continue;
                    }
                    if a_v.orientation == rcad_kernel::topods::Orientation::Internal {
                        continue;
                    }
                    b_found = true;
                    self.my_shapes_to_avoid.insert((a_e1.ptr_id(), a_e1.location, a_e1.orientation));
                } else if a_nb_e == 2 {
                    // OCCT L211-227: two edges at the vertex.
                    let a_e2 = &a_le[1];
                    // OCCT L214: aE2.IsSame(aE1) — same TShape AND Location
                    // (IsSame semantics; Orientation ignored).
                    if a_e2.is_partner(a_e1) {
                        // OCCT L216-221: TopExp::Vertices(aE1, aV1x, aV2x) —
                        // if both endpoints are the same vertex (degenerated
                        // ring), skip.
                        let vv = Self::edge_vertices(a_e1, &self.ds.locations);
                        if vv.len() >= 2 && vv[0].is_partner(&vv[1]) {
                            // Degenerated ring — both ends are the same vertex.
                            continue;
                        }
                        b_found = true;
                        self.my_shapes_to_avoid.insert((a_e1.ptr_id(), a_e1.location, a_e1.orientation));
                        self.my_shapes_to_avoid.insert((a_e2.ptr_id(), a_e2.location, a_e2.orientation));
                    }
                }
            }
            if !b_found { break; }
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

    /// Get edge endpoint vertex Shapes.
    /// OCCT TopoDS_Iterator(aE) iterates the edge's stored vertices IN STORAGE
    /// ORDER [first, last] (TopoDS_Iterator.cxx L57-70), composing the edge's
    /// orientation into each child (updateCurrentShape L72-78) and the edge's
    /// location (cumLoc, L80). The order is NOT reversed for a REVERSED edge —
    /// only the child orientations flip.
    pub(crate) fn edge_vertices(e: &Shape, locations: &[glam::DAffine3]) -> Vec<Shape> {
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
                let vf_loc =
                    crate::bop::algo::compose_edge_vertex_location(e.location, ed.first.location, locations);
                let vl_loc =
                    crate::bop::algo::compose_edge_vertex_location(e.location, ed.last.location, locations);
                if e.orientation == Orientation::Reversed {
                    vec![
                        Shape::new(ed.first.data.clone(), vf_loc, flip_ori(ed.first.orientation)),
                        Shape::new(ed.last.data.clone(), vl_loc, flip_ori(ed.last.orientation)),
                    ]
                } else {
                    vec![
                        Shape::new(ed.first.data.clone(), vf_loc, ed.first.orientation),
                        Shape::new(ed.last.data.clone(), vl_loc, ed.last.orientation),
                    ]
                }
            }
            _ => Vec::new(),
        }
    }

    /// OCCT BOPAlgo_BuilderFace::PerformLoops (BOPAlgo_BuilderFace.cxx L239-383).
    /// Builds closed wires from section edges by connecting edges at shared vertices.
    fn perform_loops(&mut self) {        // OCCT L256: aWES.SetFace(myFace)
        // OCCT L258-266: add edges to wire edge set (excluding shapes to avoid)
        let edges: Vec<Shape> = self.my_edges.iter()
            .filter(|e| !self.my_shapes_to_avoid.contains(&(e.ptr_id(), e.location, e.orientation)))
            .cloned()
            .collect();

        // OCCT L268-271: BOPAlgo_WireSplitter(aWSp) with the wire edge set.
        let a_face = self.my_face.clone().unwrap_or_else(Shape::null);
        let a_face_index = self.my_face_index.unwrap_or(usize::MAX);
        if std::env::var("RCAD_WS_DEBUG").is_ok() {
            for e in &edges {
                let hc = crate::bop::algo::wire_splitter::edge_has_pcurve(e, a_face_index, &self.ds);
                let pc = edge_pcurve_on_face(e, a_face_index, &self.ds);
                let pd = pc.as_ref().map(|(c, _, _)| {
                    let name = match c {
                        Curve2d::Line(_) => "Line", Curve2d::Circle(_) => "Circle", Curve2d::Ellipse(_) => "Ellipse",
                        Curve2d::Parabola(_) => "Parabola", Curve2d::Hyperbola(_) => "Hyperbola",
                        Curve2d::CircleInvolute(_) => "CircleInvolute", Curve2d::ArchimedeanSpiral(_) => "ArchimedeanSpiral",
                        Curve2d::LogarithmicSpiral(_) => "LogSpiral", Curve2d::SineWave(_) => "SineWave",
                        Curve2d::BSpline(_) => "BSpline", Curve2d::Bezier(_) => "Bezier", Curve2d::Trimmed(_) => "Trimmed",
                        _ => "Other",
                    };
                    match c {
                        Curve2d::Line(l) => format!("L2D o=({:.3},{:.3}) d=({:.3},{:.3}) t=({:.3},{:.3})", l.origin.x, l.origin.y, l.direction.x, l.direction.y, pc.as_ref().map(|p| p.1).unwrap_or(0.0), pc.as_ref().map(|p| p.2).unwrap_or(0.0)),
                        Curve2d::Trimmed(t) => match &*t.curve {
                            Curve2d::Line(l) => format!("TL2D o=({:.3},{:.3}) d=({:.3},{:.3}) t=({:.3},{:.3})", l.origin.x, l.origin.y, l.direction.x, l.direction.y, t.t_min, t.t_max),
                            _ => format!("TO2D({})", name),
                        },
                        _ => format!("O2D({})", name),
                    }
                }).unwrap_or_else(|| "no-pc".into());
                eprintln!("[WS-IN] face={} edge_ptr={} pcurve={} {}", a_face_index, e.ptr_id(), hc, pd);
            }
        }
        let wires = crate::bop::algo::wire_splitter::split_into_wires(&a_face, a_face_index, &edges, &self.ds);
        if std::env::var("RCAD_WS_DEBUG").is_ok() {
            eprintln!("[WS-OUT] face={} n_wires={} wire_edge_counts={:?}",
                a_face_index, wires.len(), wires.iter().map(|w| wire_edges(w).len()).collect::<Vec<_>>());
            for (wi, w) in wires.iter().enumerate() {
                let mut desc: Vec<String> = Vec::new();
                for e in wire_edges(w) {
                    let pc = edge_pcurve_on_face(&e, a_face_index, &self.ds);
                    let pd = pc.as_ref().map(|(c, _, _)| match c {
                        Curve2d::Line(l) => format!("L({:.3},{:.3})", l.direction.x, l.direction.y),
                        rcad_kernel::geom::Curve2d::BSpline(b) => {
                            if b.control_points.len() >= 2 {
                                let d = b.control_points[1] - b.control_points[0];
                                format!("BS({:.3},{:.3})", d.x, d.y)
                            } else { "BS?".into() }
                        }
                        _ => "O".into(),
                    }).unwrap_or_else(|| "no-pc".into());
                    let verts = Self::edge_vertices(&e, &self.ds.locations);
                    let vs = verts.iter().map(|v| {
                        let p = v.as_vertex().map(|vd| format!("({:.1},{:.1},{:.1})", vd.point.x, vd.point.y, vd.point.z)).unwrap_or_default();
                        format!("{}:{}", p, v.location)
                    }).collect::<Vec<_>>().join(",");
                    desc.push(format!("{}:{}{}[{}]", e.ptr_id(), pd, if e.orientation == rcad_kernel::topods::Orientation::Reversed { "R" } else { "F" }, vs));
                }
                eprintln!("[WS-LOOP] face={} wire={} edges=[{}]", a_face_index, wi, desc.join(" "));
            }
        }

        

        // OCCT L277-283: store result wires into myLoops
        self.my_loops = wires;

        // OCCT L284-321: Post-treatment — find unprocessed edges.
        // a. collect all edges that are in loops (OCCT L287-298).
        // OCCT aMEP (L285) is NCollection_Map<TopoDS_Shape> — default hasher
        // TShape + Location + Orientation.
        let mut a_mep: HashSet<(u64, u32, rcad_kernel::topods::Orientation)> = HashSet::new();
        for a_w in &self.my_loops {
            for e in wire_edges(a_w) {
                a_mep.insert((e.ptr_id(), e.location, e.orientation));
            }
        }
        // b. collect all edges that are to avoid (OCCT L304-310).
        for &key in &self.my_shapes_to_avoid {
            a_mep.insert(key);
        }
        // c. add all edges that are not processed to myShapesToAvoid (OCCT L312-321).
        for e in &self.my_edges {
            if !a_mep.contains(&(e.ptr_id(), e.location, e.orientation)) {
                self.my_shapes_to_avoid.insert((e.ptr_id(), e.location, e.orientation));
            }
        }

        // OCCT L327-382: 2. Internal Wires — build wires from the avoided
        // edges, connecting them at shared vertices.
        self.my_loops_internal.clear();
        // OCCT L330-335: aVEMap — vertex -> [avoid edges] via
        // MapShapesAndAncestors(aEE, VERTEX, EDGE). aVEMap (L244-245) uses
        // TopTools_ShapeMapHasher — key identity TShape + Location.
        let avoid_edges: Vec<Shape> = self.my_edges.iter()
            .filter(|e| self.my_shapes_to_avoid.contains(&(e.ptr_id(), e.location, e.orientation)))
            .cloned()
            .collect();
        let mut a_ve_map: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
        for a_ee in &avoid_edges {
            for a_v in Self::edge_vertices(a_ee, &self.ds.locations) {
                let l = a_ve_map
                    .entry((a_v.ptr_id(), a_v.location))
                    .or_default();
                // OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-107) appends
                // EVERY ancestor edge occurrence — the FORWARD and REVERSED
                // copies of the same TShape are distinct ancestors, both kept.
                l.push(a_ee.clone());
            }
        }
        // OCCT L337-382: per-start wire growth with the aMAdded fence; the
        // loop stops (bFlag) once every avoided edge is collected. OCCT
        // aMAdded (L246) is NCollection_Map<TopoDS_Shape> with the DEFAULT
        // hasher — TopoDS_Shape::HashCode includes TShape + Location +
        // Orientation, so the fence is orientation-sensitive.
        let a_nb_ea = avoid_edges.len();
        let mut a_m_added: HashSet<(u64, u32, rcad_kernel::topods::Orientation)> = HashSet::new();
        let mut b_flag = true;
        for a_ee in &avoid_edges {
            if !b_flag { break; }
            if !a_m_added.insert((a_ee.ptr_id(), a_ee.location, a_ee.orientation)) { continue; }
            // OCCT L350-353: make new wire and add the start edge.
            let mut a_w: Vec<Shape> = vec![a_ee.clone()];
            // OCCT L355-379: grow through the vertices of the added edges.
            let mut i = 0;
            while i < a_w.len() && b_flag {
                let a_e = &a_w[i];
                for a_v in Self::edge_vertices(a_e, &self.ds.locations) {
                    if let Some(a_le) = a_ve_map.get(&(a_v.ptr_id(), a_v.location)) {
                        for a_ex in a_le {
                            if a_m_added.insert((a_ex.ptr_id(), a_ex.location, a_ex.orientation)) {
                                a_w.push(a_ex.clone());
                                if a_m_added.len() == a_nb_ea {
                                    b_flag = false;
                                }
                            }
                        }
                    }
                }
                i += 1;
            }
            // OCCT L380-381: aW.Closed(BRep_Tool::IsClosed(aW));
            // myLoopsInternal.Append(aW) — the wire shape carries the CLOSED
            // flag (BOPAlgo_WireSplitter::MakeWire semantics).
            self.my_loops_internal
                .push(crate::bop::algo::wire_splitter::make_wire(&a_w, &self.ds.locations));
        }
    }

    /// OCCT BOPAlgo_BuilderFace::PerformAreas (BOPAlgo_BuilderFace.cxx L387-613).
    fn perform_areas(&mut self) {
        self.my_areas.clear();
        let (a_surf_opt, a_tol) = match self.my_face.as_ref().and_then(|f| match &*f.data {
            TShape::Face(fd) => Some((fd.surface.clone(), fd.tolerance)),
            _ => None,
        }) { Some(v) => v, None => return };
        let a_surf = match a_surf_opt { Some(s) => s, None => return };
        // OCCT L401-414: empty loops — an infinite face becomes a face without
        // wires.
        if self.my_loops.is_empty() {
            let is_infinite = self
                .my_face
                .as_ref()
                .map(Self::is_infinite_face)
                .unwrap_or(false);
            if is_infinite {
                let natural_restriction = self
                    .my_face
                    .as_ref()
                    .and_then(|f| match &*f.data {
                        TShape::Face(fd) => Some(fd.natural_restriction),
                        _ => None,
                    })
                    .unwrap_or(false);
                let face_tshape = TShape::Face(rcad_kernel::topods::TFaceData {
                    my_shapes: vec![],
                    flags: tshape_flags::DEFAULT,
                    surface: Some(a_surf.clone()),
                    surface_location: 0,
                    outer_wire: Shape::null(),
                    inner_wires: vec![],
                    sample_point: None,
                    uv_domain: None,
                    internal_vertices: vec![],
                    tolerance: a_tol,
                    natural_restriction,
                });
                self.my_areas.push(Shape::new(
                    std::sync::Arc::new(face_tshape),
                    0,
                    rcad_kernel::topods::Orientation::Forward,
                ));
            }
            return;
        }

        // OCCT L417-423: growth faces + hole faces + hole edge map
        let mut a_new_faces: Vec<Shape> = Vec::new();
        let mut a_hole_faces: Vec<Shape> = Vec::new();
        let mut a_mhe: HashSet<(u64, u32)> = HashSet::new();

        // OCCT L427-458: classify each loop
        for a_wire in &self.my_loops {
            let loop_edges = wire_edges(a_wire);
            // OCCT L437-439: aBB.MakeFace(aFace, aS, aLoc, aTol);
            // aBB.Add(aFace, aWire) — the loop wire shape is added as-is.
            let face_tshape = TShape::Face(rcad_kernel::topods::TFaceData {
                my_shapes: vec![], flags: tshape_flags::DEFAULT,
                surface: Some(a_surf.clone()), surface_location: 0,
                outer_wire: a_wire.clone(), inner_wires: vec![],
                sample_point: None, uv_domain: None,
                internal_vertices: vec![], tolerance: a_tol,
                natural_restriction: false,
            });
            let a_face = Shape::new(
                std::sync::Arc::new(face_tshape),
                0, rcad_kernel::topods::Orientation::Forward,
            );

            // OCCT L441-447: IsGrowthWire + FClass2d::IsHole
            let b_is_growth = {
                // OCCT L441: IsGrowthWire(aWire, aMHE) — returns true when the
                // wire contains any hole-face edge (BOPAlgo_BuilderFace.cxx
                // L898-913: theMHE.Contains(aIt.Value())). Only when it has no
                // hole edge is the FClass2d classification run.
                let has_hole_edge = loop_edges.iter().any(|e| a_mhe.contains(&(e.ptr_id(), e.location)));
                if has_hole_edge {
                    true
                } else {
                    // OCCT L445-446: FClass2d(aFace).IsHole() — aFace is the
                    // temporary face built from the analyzed loop wire.
                    let fi = self.my_face_index.unwrap_or(0);
                    let is_hole = self.my_context.fclass2d_is_hole(self.ds, fi, &loop_edges);
                    if std::env::var("RCAD_AREA_DEBUG").is_ok() {
                        eprintln!("[AREA] face={} loop n_edges={} is_hole={}", fi, loop_edges.len(), is_hole);
                    }
                    !is_hole
                }
            };

            // OCCT L450-458: save growth vs hole
            if b_is_growth {
                a_new_faces.push(a_face);
            } else {
                a_hole_faces.push(a_face);
                for e in loop_edges { a_mhe.insert((e.ptr_id(), e.location)); }
            }
        }

        // OCCT L461-466: no holes
        if a_hole_faces.is_empty() {
            self.my_areas = a_new_faces;
            return;
        }

        // OCCT L468-540: classify the holes relative to the growth faces.
        // aHoleFaceMap: hole index -> growth index (most specific / innermost).
        let fi = self.my_face_index.unwrap_or(0);
        let a_nb_h = a_hole_faces.len();
        // OCCT L487/L540: aHoleFaceMap and aFaceHolesMap are IndexedDataMap —
        // the back-map iteration (L544) follows insertion order.
        let mut a_hole_face_map: IndexMap<usize, usize> = IndexMap::new();
        // OCCT L470-491: prepare the hole UV boxes (BOPTools_Box2dTree). rcad
        // scans linearly; the tree is a performance structure, the overlap
        // predicate is preserved.
        let hole_boxes: Vec<Option<[f64; 4]>> = a_hole_faces
            .iter()
            .map(|hf| wire_uv_bounds(self.ds, fi, &face_wire_edges(hf)))
            .collect();
        for (gfi, face) in a_new_faces.iter().enumerate() {
            let g_edges = face_wire_edges(face);
            let g_box = match wire_uv_bounds(self.ds, fi, &g_edges) {
                Some(b) => b,
                None => continue,
            };
            for (hi, hf) in a_hole_faces.iter().enumerate() {
                // OCCT: selector — candidate holes whose UV box intersects.
                let h_box = match hole_boxes[hi] {
                    Some(b) => b,
                    None => continue,
                };
                if !boxes_overlap(g_box, h_box) {
                    continue;
                }
                // OCCT L518: if (!IsInside(aHole, aFace, myContext)) continue;
                let h_edges = face_wire_edges(hf);
                if !is_inside_wire(self.ds, fi, &h_edges, &g_edges) {
                    continue;
                }
                // OCCT L522-533: keep the most specific (innermost) face.
                match a_hole_face_map.get(&hi) {
                    Some(&old_gfi) => {
                        let old_g_edges = face_wire_edges(&a_new_faces[old_gfi]);
                        if is_inside_wire(self.ds, fi, &g_edges, &old_g_edges) {
                            a_hole_face_map.insert(hi, gfi);
                        }
                    }
                    None => {
                        a_hole_face_map.insert(hi, gfi);
                    }
                }
            }
        }
        // OCCT L536-553: back map — face -> its holes.
        let mut a_face_holes: HashMap<usize, Vec<usize>> = HashMap::new();
        for (hi, gfi) in &a_hole_face_map {
            a_face_holes.entry(*gfi).or_default().push(*hi);
        }
        // OCCT L556-580: unused holes — a new face if the original face is
        // unbounded (BRepBndLib box has open sides; rcad: infinite face).
        if a_nb_h != a_hole_face_map.len() {
            let is_open = self
                .my_face
                .as_ref()
                .map(Self::is_infinite_face)
                .unwrap_or(false);
            if is_open {
                let mut an_unused: Vec<Shape> = Vec::new();
                for (hi, hf) in a_hole_faces.iter().enumerate() {
                    if !a_hole_face_map.contains_key(&hi) {
                        if let TShape::Face(hfd) = &*hf.data {
                            an_unused.push(hfd.outer_wire.clone());
                        }
                    }
                }
                let a_face = Shape::new(
                    std::sync::Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
                        my_shapes: vec![], flags: tshape_flags::DEFAULT,
                        surface: Some(a_surf.clone()), surface_location: 0,
                        outer_wire: Shape::null(),
                        inner_wires: an_unused,
                        sample_point: None, uv_domain: None,
                        internal_vertices: vec![], tolerance: a_tol,
                        natural_restriction: false,
                    })),
                    0, rcad_kernel::topods::Orientation::Forward,
                );
                a_new_faces.push(a_face);
            }
        }
        // OCCT L583-613: add the holes to the faces and append to myAreas.
        let mut result_faces: Vec<Shape> = Vec::new();
        for (gfi, face) in a_new_faces.iter().enumerate() {
            let mut a_face = face.clone();
            if let Some(holes) = a_face_holes.get(&gfi) {
                let mut inner_wires: Vec<Shape> = Vec::new();
                for &hi in holes {
                    if let TShape::Face(hfd) = &*a_hole_faces[hi].data {
                        inner_wires.push(hfd.outer_wire.clone());
                    }
                }
                if let TShape::Face(fd) = &*a_face.data {
                    a_face = Shape::new(
                        std::sync::Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
                            my_shapes: vec![], flags: tshape_flags::DEFAULT,
                            surface: fd.surface.clone(), surface_location: 0,
                            outer_wire: fd.outer_wire.clone(),
                            inner_wires, sample_point: None, uv_domain: None,
                            internal_vertices: vec![], tolerance: a_tol,
                            natural_restriction: false,
                        })),
                        0, rcad_kernel::topods::Orientation::Forward,
                    );
                }
            }
            result_faces.push(a_face);
        }
        // OCCT L543-613: internal wires
        self.my_areas = result_faces;
    }

    /// OCCT BOPAlgo_BuilderFace::PerformInternalShapes (BuilderFace.cxx L618-750).
    /// Classifies the internal loops (myLoopsInternal) relatively the area faces
    /// and adds them as internal wires to the faces that contain them.
    fn perform_internal_shapes(&mut self) {
        // OCCT L620-622: myAvoidInternalShapes.
        if self.my_avoid_internal_shapes {
            return;
        }
        // OCCT L624-626: myLoopsInternal.IsEmpty().
        if self.my_loops_internal.is_empty() {
            return;
        }
        let fi = self.my_face_index.unwrap_or(0);
        // OCCT L638-664: map of internal edges + their UV boxes. OCCT
        // anEdgesMap (L637) is IndexedMap with TopTools_ShapeMapHasher —
        // key identity TShape + Location.
        let mut edges_map: Vec<Shape> = Vec::new();
        let mut edges_idx: HashMap<(u64, u32), usize> = HashMap::new();
        let mut edge_boxes: Vec<[f64; 4]> = Vec::new();
        for wire in &self.my_loops_internal {
            for e in wire_edges(wire) {
                if edges_idx.contains_key(&(e.ptr_id(), e.location)) {
                    continue;
                }
                let box_e = match edge_pcurve_on_face(&e, fi, self.ds) {
                    // OCCT L645-649: BRepTools::AddUVBounds(myFace, aE, aBoxE)
                    // — exact UV bounds via GeomBndLib_Curve2d.
                    Some((pc, t1, t2)) => match add_uv_bounds(self.ds, fi, &pc, t1, t2) {
                        Some(b) => b,
                        None => continue,
                    },
                    None => continue,
                };
                edges_idx.insert((e.ptr_id(), e.location), edges_map.len());
                edges_map.push(e.clone());
                edge_boxes.push(box_e);
            }
        }
        if edges_map.is_empty() {
            return;
        }
        // OCCT L673-740: classify the edges relatively the area faces.
        let a_medone: HashSet<usize> = HashSet::new();
        let mut a_medone = a_medone;
        let mut a_face_holes: HashMap<usize, Vec<Shape>> = HashMap::new();
        for (ai, face) in self.my_areas.iter().enumerate() {
            let f_edges = face_wire_edges(face);
            let f_box = match wire_uv_bounds(self.ds, fi, &f_edges) {
                Some(b) => b,
                None => continue,
            };
            // OCCT L687-708: collect the edges inside the face.
            let mut edges_inside: Vec<Shape> = Vec::new();
            for (ei, e) in edges_map.iter().enumerate() {
                if a_medone.contains(&ei) {
                    continue;
                }
                // OCCT: BOPTools_Box2dTreeSelector — candidate via UV box.
                if !boxes_overlap(f_box, edge_boxes[ei]) {
                    continue;
                }
                // OCCT L703: if (IsInside(aE, aF, myContext)).
                if is_inside_wire(self.ds, fi, &[e.clone()], &f_edges) {
                    edges_inside.push(e.clone());
                    a_medone.insert(ei);
                }
            }
            if edges_inside.is_empty() {
                continue;
            }
            // OCCT L712: MakeInternalWires(anEdgesInside, aLSI).
            let a_lsi = make_internal_wires(&edges_inside, &self.ds.locations);
            a_face_holes.entry(ai).or_default().extend(a_lsi);
            // OCCT L730-736: early exit when all edges are classified.
            if a_medone.len() == edges_map.len() {
                break;
            }
        }
        // Add the internal wires to the faces (OCCT aBB.Add(aF, aWI)).
        let mut result_faces: Vec<Shape> = Vec::new();
        for (ai, face) in self.my_areas.iter().enumerate() {
            let mut a_face = face.clone();
            if let Some(wires) = a_face_holes.get(&ai) {
                let mut inner_wires: Vec<Shape> = Vec::new();
                if let TShape::Face(fd) = &*a_face.data {
                    inner_wires = fd.inner_wires.clone();
                }
                for w in wires {
                    inner_wires.push(w.clone());
                }
                if let TShape::Face(fd) = &*a_face.data {
                    a_face = Shape::new(
                        std::sync::Arc::new(TShape::Face(rcad_kernel::topods::TFaceData {
                            my_shapes: vec![], flags: tshape_flags::DEFAULT,
                            surface: fd.surface.clone(), surface_location: 0,
                            outer_wire: fd.outer_wire.clone(),
                            inner_wires, sample_point: None, uv_domain: None,
                            internal_vertices: vec![], tolerance: fd.tolerance,
                            natural_restriction: fd.natural_restriction,
                        })),
                        0, rcad_kernel::topods::Orientation::Forward,
                    );
                }
            }
            result_faces.push(a_face);
        }
        // OCCT L742-750: unused edges — warning (rcad: edges left unused are
        // simply not added; the alert is not modelled).
        self.my_areas = result_faces;
    }
}

/// OCCT TopoDS_Iterator(theW) — the wire's edges with the wire's orientation
/// composed (cumOri). The wires built by this pipeline are FORWARD, so the
/// stored edge orientations are returned unchanged.
pub(crate) fn wire_edges(w: &Shape) -> Vec<Shape> {
    use rcad_kernel::topods::Orientation;
    match &*w.data {
        TShape::Wire(wd) => {
            if w.orientation == Orientation::Reversed {
                wd.edges
                    .iter()
                    .map(|e| {
                        let mut c = e.clone();
                        c.orientation = match c.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            other => other,
                        };
                        c
                    })
                    .collect()
            } else {
                wd.edges.clone()
            }
        }
        _ => Vec::new(),
    }
}

/// OCCT MakeInternalWires (BOPAlgo_BuilderFace.cxx L782-838) — groups the
/// connected internal edges into wires (edges set with INTERNAL orientation);
/// each wire gets aW.Closed(BRep_Tool::IsClosed(aW)) (L835).
/// aMVE (L788-789) and aAddedMap (L786-787) use TopTools_ShapeMapHasher —
/// key identity TShape + Location, orientation ignored.
fn make_internal_wires(edges: &[Shape], locations: &[glam::DAffine3]) -> Vec<Shape> {
    use rcad_kernel::topods::Orientation;
    // aMVE: vertex -> edges.
    let mut a_mve: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();
    for e in edges {
        for v in BuilderFace::edge_vertices(e, locations) {
            a_mve
                .entry((v.ptr_id(), v.location))
                .or_default()
                .push(e.clone());
        }
    }
    let mut a_added: HashSet<(u64, u32)> = HashSet::new();
    let mut wires: Vec<Shape> = Vec::new();
    for e in edges {
        if !a_added.insert((e.ptr_id(), e.location)) {
            continue;
        }
        // OCCT L810: aEE.Orientation(TopAbs_INTERNAL) — the start edge is
        // added with INTERNAL orientation.
        let mut a_w: Vec<Shape> = vec![{
            let mut e0 = e.clone();
            e0.orientation = Orientation::Internal;
            e0
        }];
        // Grow the wire through shared vertices.
        let mut i = 0;
        while i < a_w.len() {
            let cur = &a_w[i];
            for v in BuilderFace::edge_vertices(cur, locations) {
                if let Some(le) = a_mve.get(&(v.ptr_id(), v.location)) {
                    for ex in le {
                        if a_added.insert((ex.ptr_id(), ex.location)) {
                            let mut a_el = ex.clone();
                            a_el.orientation = Orientation::Internal;
                            a_w.push(a_el);
                        }
                    }
                }
            }
            i += 1;
        }
        // OCCT L835: aW.Closed(BRep_Tool::IsClosed(aW)).
        wires.push(crate::bop::algo::wire_splitter::make_wire(&a_w, locations));
    }
    wires
}

// ---- hole-classification helpers (OCCT BOPAlgo_BuilderFace.cxx L468-613) ----

/// OCCT BRep_Tool::CurveOnSurface(aE, aF) — the edge's pcurve on the face,
/// keyed by the face's TShape identity (ptr_id, location, with the
/// DS-canonical fallback).
fn edge_pcurve_on_face(e: &Shape, face_index: usize, ds: &DS) -> Option<(Curve2d, f64, f64)> {
    // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the pcurve key
    // location is L.Predivided(E.Location()).
    let (fid, floc) = ds.face_key(face_index)?;
    let fkey = (
        fid,
        crate::bop::algo::compose_face_edge_pcurve_location(floc, e.location, &ds.locations),
    );
    match &*e.data {
        TShape::Edge(ed) => {
            // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L354-361): a closed
            // surface seam edge (CurveOnClosedSurface) returns the second pcurve
            // for a REVERSED edge and the first otherwise — the two wire
            // instances of the seam map to u=2*PI and u=0.
            if let Some(r) = ed.representations.iter().find_map(|r| match r {
                rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    face,
                    pcurve1,
                    pcurve2,
                    range,
                } if *face == fkey => Some((pcurve1.clone(), pcurve2.clone(), *range)),
                _ => None,
            }) {
                let (pc1, pc2, range) = r;
                // The BuilderFace face is normalized to FORWARD (SetFace),
                // so the face_reversed term of wire_splitter::edge_pcurve is
                // always false here.
                let pc = if e.orientation == rcad_kernel::topods::Orientation::Reversed {
                    pc2
                } else {
                    pc1
                };
                return Some((pc, range[0], range[1]));
            }
            if let Some(v) = ed.pcurves.get(&fkey) {
                return Some(v.clone());
            }
            if let Some(idx) = ds.map_shape_index.get(&(e.ptr_id(), e.location)) {
                if let Some(ed2) = ds.shape_info(*idx).shape.as_edge() {
                    if let Some(v) = ed2.pcurves.get(&fkey) {
                        return Some(v.clone());
                    }
                }
            }
            None
        }
        _ => None,
    }
}

/// OCCT BRepTools::AddUVBounds (BRepTools.cxx L172-362) — the exact UV bounds
/// of the edge pcurve on the face: GeomBndLib_Curve2d::Box (via
/// BndLib_Add2dCurve::Add with tolerance 0), then clamped to the face UV
/// domain for non-periodic surfaces. Returns [u_min, u_max, v_min, v_max].
fn add_uv_bounds(ds: &DS, face_index: usize, pc: &Curve2d, t1: f64, t2: f64) -> Option<[f64; 4]> {
    let face = ds.shape(face_index);
    let (uv_domain, surf) = match &*face.data {
        TShape::Face(fd) => (fd.uv_domain, fd.surface.clone()?),
        _ => return None,
    };
    // OCCT L185: BndLib_Add2dCurve::Add(aC2D, aT1, aT2, 0., aBoxC).
    let mut b = rcad_kernel::curve2d_bounding_box(pc, t1, t2, 0.0);
    if let Some(uvd) = uv_domain {
        // OCCT L202-281 (U) / L283-361 (V): clamp to the face bounds when the
        // surface is not periodic. The OCCT B-spline periodicity verification
        // (L210-268) is approximated by the surface's declared periodicity.
        if !surf.is_u_periodic() {
            let (a_umin, a_umax) = (uvd[0], uvd[1]);
            if b[0] < a_umin && a_umin < b[1] {
                b[0] = a_umin;
            }
            if b[0] < a_umax && a_umax < b[1] {
                b[1] = a_umax;
            }
        }
        if !surf.is_v_periodic() {
            let (a_vmin, a_vmax) = (uvd[2], uvd[3]);
            if b[2] < a_vmin && a_vmin < b[3] {
                b[2] = a_vmin;
            }
            if b[2] < a_vmax && a_vmax < b[3] {
                b[3] = a_vmax;
            }
        }
    }
    Some(b)
}

/// The edges of a face's outer wire plus all inner wires.
/// OCCT BOPAlgo_BuilderFace.cxx L850-851: `TopExp::MapShapes(theF, TopAbs_EDGE,
/// aFaceEdgesMap)` — the whole face, including holes. The caller
/// perform_internal_shapes runs after PerformAreas has attached holes, so
/// only outer-wire edges would miss shared-edge detection and shrink the UV
/// box (AddUVBounds(aF, aBoxF), L685-686).
fn face_wire_edges(face: &Shape) -> Vec<Shape> {
    match &*face.data {
        TShape::Face(fd) => {
            let mut edges: Vec<Shape> = Vec::new();
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                edges.extend(wd.edges.iter().cloned());
            }
            for w in &fd.inner_wires {
                if let TShape::Wire(wd) = &*w.data {
                    edges.extend(wd.edges.iter().cloned());
                }
            }
            edges
        }
        _ => Vec::new(),
    }
}

/// Sample a wire's pcurves into a closed 2D polygon and its UV box
/// (OCCT BRepTools::AddUVBounds + the wire's 2D boundary).
fn wire_uv_polygon(ds: &DS, face_index: usize, edges: &[Shape]) -> Option<(Vec<DVec2>, [f64; 4])> {
    let mut pts: Vec<DVec2> = Vec::new();
    let mut b: [f64; 4] = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    for e in edges {
        let (pc, t1, t2) = match edge_pcurve_on_face(e, face_index, ds) {
            Some(v) => v,
            None => continue,
        };
        let n = match pc {
            Curve2d::Line(_) => 1,
            Curve2d::Circle(_) | Curve2d::Ellipse(_) => 6,
            _ => 3,
        };
        for k in 0..=n {
            let t = t1 + (t2 - t1) * (k as f64) / (n as f64);
            let p = Curve2dEval::point_at(&pc, t);
            b[0] = b[0].min(p.x);
            b[1] = b[1].max(p.x);
            b[2] = b[2].min(p.y);
            b[3] = b[3].max(p.y);
            pts.push(p);
        }
    }
    if pts.len() < 3 {
        return None;
    }
    Some((pts, b))
}

/// The UV bounds of a wire (min/max over the sampled pcurves).
fn wire_uv_bounds(ds: &DS, face_index: usize, edges: &[Shape]) -> Option<[f64; 4]> {
    wire_uv_polygon(ds, face_index, edges).map(|(_, b)| b)
}

fn boxes_overlap(a: [f64; 4], b: [f64; 4]) -> bool {
    a[0] <= b[1] && b[0] <= a[1] && a[2] <= b[3] && b[2] <= a[3]
}

/// OCCT IsInside (BOPAlgo_BuilderFace.cxx L842-897) — the wire `the_wire` is
/// inside the face `the_f` when the first non-degenerated edge of the wire, not
/// shared with the face, has its pcurve midpoint classified IN the face by the
/// face classifier (IntTools_FClass2d::Perform, L890-891). rcad classifies the
/// midpoint with FClass2d (brep_top_adaptor), not a sampled UV polygon.
///
/// OCCT theF is the temporary face built from the GROWTH loop wire
/// (BOPAlgo_BuilderFace.cxx L437-445), and IntTools_Context::FClass2d(aF)
/// (IntTools_Context.cxx L225-242) builds the classifier from that face's
/// (single) loop — so `against_edges` (the growth loop's edges) must feed the
/// classifier, not the original DS face. For the innermost-hole check (L548)
/// theF is the previously-assigned growth face, so the same `against_edges`
/// semantic applies.
fn is_inside_wire(ds: &DS, face_index: usize, wire_edges: &[Shape], against_edges: &[Shape]) -> bool {
    // OCCT BOPAlgo_BuilderFace::IsInside (BuilderFace.cxx L842-894).
    // OCCT L850-851: aFaceEdgesMap = TopExp::MapShapes(theF, EDGE) with the
    // identity hasher (TShape + Location, orientation ignored).
    let a_face_edges: std::collections::HashSet<(u64, u32)> = against_edges
        .iter()
        .map(|e| (e.ptr_id(), e.location))
        .collect();
    // OCCT L860-893: iterate the wire edges; the first edge that is not
    // degenerated, is not contained in the face and has a pcurve decides.
    for e in wire_edges {
        // OCCT L864-868: skip degenerated edges.
        let degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(true);
        if degen {
            continue;
        }
        // OCCT L870-875: the face contains the edge from the wire, thus the
        // wire cannot be inside that face.
        if a_face_edges.contains(&(e.ptr_id(), e.location)) {
            return false;
        }
        // OCCT L877-883: get the 2D curve of the edge on the face; skip if
        // the curve is null.
        let (pc, t1, t2) = match edge_pcurve_on_face(e, face_index, ds) {
            Some(v) => v,
            None => continue,
        };
        // OCCT L885-891: classify the middle point; the first classification
        // decides the result.
        let p = Curve2dEval::point_at(&pc, 0.5 * (t1 + t2));
        let fclass2d = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new_for_loop(
            ds,
            face_index,
            ds.face_tolerance(face_index),
            against_edges,
        );
        return fclass2d.perform(ds, p, true)
            == crate::topalgo::brep_top_adaptor::fclass2d::State::In;
    }
    false
}

// ============================================================================
// Tests
// ============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::topods::{TFaceData, Orientation};

    fn wire(edges: Vec<Shape>) -> Shape {
        Shape::new(
            Arc::new(TShape::Wire(TWireData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                edges,
            })),
            0,
            Orientation::Forward,
        )
    }

    /// Regression test for perform_internal_shapes (audit item #3).
    ///
    /// OCCT IsInside (BOPAlgo_BuilderFace.cxx L850-851) builds the face edge
    /// set with `TopExp::MapShapes(theF, TopAbs_EDGE)` — the WHOLE face,
    /// including the inner wires (holes) attached by PerformAreas; and the UV
    /// box comes from `BRepTools::AddUVBounds(aF, aBoxF)` (L685-686), also
    /// over all wires. rcad's face_wire_edges used to return only the outer
    /// wire's edges, so a hole edge shared with an internal edge was not
    /// detected (is_inside_wire returned the wrong classification) and the UV
    /// box was too small. The fix returns outer + inner edges.
    #[test]
    fn face_wire_edges_includes_inner_wires() {
        let e_outer = Shape::null();
        let e_inner = Shape::null();
        let face = Shape::new(
            Arc::new(TShape::Face(TFaceData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                surface: None,
                surface_location: 0,
                outer_wire: wire(vec![e_outer.clone()]),
                inner_wires: vec![wire(vec![e_inner.clone()])],
                sample_point: None,
                uv_domain: None,
                internal_vertices: vec![],
                tolerance: 0.0,
                natural_restriction: false,
            })),
            0,
            Orientation::Forward,
        );
        let edges = face_wire_edges(&face);
        assert_eq!(edges.len(), 2, "must include outer + inner wire edges");
        assert_eq!(edges[0].ptr_id(), e_outer.ptr_id());
        assert_eq!(edges[1].ptr_id(), e_inner.ptr_id());
    }
}
