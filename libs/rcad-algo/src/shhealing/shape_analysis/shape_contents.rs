//! OCCT ShapeAnalysis package class (TKShHealing):
//! `ShapeAnalysis_ShapeContents` (`ShapeAnalysis_ShapeContents.hxx`
//! L17-217 + `.cxx` L1-358).
//!
//! Counts the topological and geometric contents of a shape (solids,
//! shells, faces, wires, edges, vertices, big splines, C0 objects, offset
//! objects, seams, shared sub-shapes, free faces/wires/edges) with the
//! optional per-case sections.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `BRep_Tool::Surface/Curve/CurveOnSurface/
//!    IsClosed` and `TopExp_Explorer` read the TShape graph through
//!    `rcad_kernel::BRep` (the edge.rs bridge #1).
//! 2. `occ::handle<NCollection_HSequence<TopoDS_Shape>>` -> `Vec<Shape>`.
//! 3. `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> mapsh` ->
//!    `HashSet<(u64, u32)>` keyed by the IsSame identity (TShape pointer +
//!    location); the extent reads are the set lengths.
//! 4. `shape.Location(TopLoc_Location())` — the OCCT copy assignment with
//!    the identity location -> the rcad `location = 0` on the local clone.
//! 5. `Geom_Surface::IsCNu(N)/IsCNv(N)` and `Geom_Curve::IsCN(N)` — only
//!    the BSpline carriers can decline; the rcad expanded-knot vector
//!    recomputes the interior-knot multiplicities (the local
//!    `bspline_is_cn`, the Geom_BSplineSurface::IsCNu bridge).
//! 6. `Geom_ElementarySurface::Position().Direct()` — the rcad elementary
//!    surface carriers carry no handedness (right-handed frames by
//!    construction), so Direct() reads as true; the left-handed swept
//!    cylinder (the `y_dir` override) cannot be queried and stays in the
//!    OCCT Direct branch (architecture difference).
//! 7. `TopExp_Explorer(S, T, A)` (the avoid form) -> the local
//!    `topexp_explorer_avoid` (the walk stops before entering shapes of
//!    the avoided type).

use std::collections::HashSet;

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, CurveRepresentation, ShapeType, TShape};

use crate::shhealing::shape_build::brep_tool::topexp_explorer;

/// OCCT BRep_Tool::Surface(face, loc) — the surface and the location.
fn brep_tool_surface_loc(fac: &Shape) -> (Option<rcad_kernel::geom::Surface3>, u32) {
    match fac.data.as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fac.location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::Curve(edge, first, last) — the no-location variant.
fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed
            .curve
            .as_ref()
            .map(|c| (c.clone(), ed.range[0], ed.range[1])),
        _ => None,
    }
}

/// OCCT BRep_Tool::Curve(edge, loc, first, last) — the out-location
/// variant (the edge wrapper location; edge.rs bridge #5).
fn brep_tool_curve_loc(edg: &Shape) -> (Option<Curve3>, u32, f64, f64) {
    match edg.data.as_ref() {
        TShape::Edge(ed) => (
            ed.curve.clone(),
            edg.location,
            ed.range[0],
            ed.range[1],
        ),
        _ => (None, 0, 0.0, 0.0),
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edge, face, first, last) — the stored
/// overload matched by the face TShape pointer.
fn brep_tool_curve_on_surface_face(
    the_edge: &Shape,
    the_face: &Shape,
) -> Option<(Curve2d, f64, f64)> {
    let fptr = std::sync::Arc::as_ptr(&the_face.data) as u64;
    let reps = match the_edge.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return None,
    };
    for r in reps {
        let (key, pc, range) = match r {
            CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        if key.0 == fptr {
            return Some((pc.clone(), range[0], range[1]));
        }
    }
    None
}

/// OCCT BRep_Tool::IsClosed(edge, face) — the seam test: a
/// BRep_CurveOnClosedSurface representation for the face (the edge.rs
/// re-host over the pool registry).
fn brep_tool_is_closed_edge_face(brep: &BRep, edg: &Shape, fac: &Shape) -> bool {
    let fsurf = match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    };
    let expected_loc =
        rcad_kernel::topods::compose_pcurve_location(fac.location, edg.location, &brep.locations);
    for r in edge_representations(edg) {
        if let CurveRepresentation::CurveOnClosedSurface { face: (fptr, lhash), .. } = r {
            if *lhash != expected_loc {
                continue;
            }
            // OCCT compares the surface handle; rcad compares the surface
            // value when the pointer row resolves.
            match (face_surface_by_ptr(brep, *fptr), fsurf.as_ref()) {
                (Some(s), Some(fs)) => {
                    if rcad_kernel::topods::surface_same(&s, fs) {
                        return true;
                    }
                }
                (None, _) => return true,
                _ => {}
            }
        }
    }
    false
}

/// The edge's curve representations (empty for non-edges).
fn edge_representations(edg: &Shape) -> &[CurveRepresentation] {
    match edg.data.as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => &[],
    }
}

/// The face surface registered in the pool under a TShape pointer (the
/// Geom_Surface handle identity stand-in; the kernel face_surface_by_ptr
/// precedent, topods.rs L2236).
fn face_surface_by_ptr(brep: &BRep, fptr: u64) -> Option<rcad_kernel::geom::Surface3> {
    let ts = brep
        .tshapes
        .iter()
        .find(|ts| std::sync::Arc::as_ptr(ts) as u64 == fptr)?;
    match ts.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT Geom_BSplineSurface::IsCNu(N) / IsCNv(N) (and Geom_BSplineCurve::
/// IsCN(N)) over the rcad expanded knot vector: every interior knot
/// multiplicity must not exceed `degree - N`.
fn bspline_is_cn(knots: &[f64], degree: usize, n: i32) -> bool {
    if knots.is_empty() {
        return false;
    }
    let eps = f64::EPSILON * 100.0;
    let mut i = 0usize;
    let len = knots.len();
    while i < len {
        // Group the equal knots (the expanded multiplicity).
        let mut j = i;
        while j + 1 < len && (knots[j + 1] - knots[i]).abs() <= eps {
            j += 1;
        }
        let mult = j - i + 1;
        let is_interior = i > 0 && j + 1 < len;
        if is_interior && mult as i32 > degree as i32 - n {
            return false;
        }
        i = j + 1;
    }
    true
}

/// OCCT TopExp_Explorer(S, T, A) — the avoid form: collects the shapes of
/// type `to_find` without descending into the shapes of type `to_avoid`.
fn topexp_explorer_avoid(
    brep: &mut BRep,
    shape: &Shape,
    to_find: ShapeType,
    to_avoid: ShapeType,
) -> Vec<Shape> {
    let mut out = Vec::new();
    let mut stack: Vec<Shape> = vec![shape.clone()];
    while let Some(cur) = stack.pop() {
        let st = cur.shape_type();
        if st == to_find {
            out.push(cur.clone());
        }
        if st == to_avoid {
            // Do not descend into the avoided sub-shape.
            continue;
        }
        // Descend into the children (composed orientations; the OCCT
        // explorer walk).
        let children: Vec<Shape> = crate::shhealing::shape_build::brep_tool::raw_subshapes(
            brep,
            &cur,
        );
        for child in children {
            stack.push(child);
        }
    }
    out
}

/// OCCT ShapeAnalysis_ShapeContents (hxx L30-217).
pub struct ShapeAnalysisShapeContents {
    // The counters (OCCT hxx L175-207).
    /// OCCT `myNbSolids`.
    my_nb_solids: i32,
    /// OCCT `myNbShells`.
    my_nb_shells: i32,
    /// OCCT `myNbFaces`.
    my_nb_faces: i32,
    /// OCCT `myNbWires`.
    my_nb_wires: i32,
    /// OCCT `myNbEdges`.
    my_nb_edges: i32,
    /// OCCT `myNbVertices`.
    my_nb_vertices: i32,
    /// OCCT `myNbSolidsWithVoids`.
    my_nb_solids_with_voids: i32,
    /// OCCT `myNbBigSplines`.
    my_nb_big_splines: i32,
    /// OCCT `myNbC0Surfaces`.
    my_nb_c0_surfaces: i32,
    /// OCCT `myNbC0Curves`.
    my_nb_c0_curves: i32,
    /// OCCT `myNbOffsetSurf`.
    my_nb_offset_surf: i32,
    /// OCCT `myNbIndirectSurf`.
    my_nb_indirect_surf: i32,
    /// OCCT `myNbOffsetCurves`.
    my_nb_offset_curves: i32,
    /// OCCT `myNbTrimmedCurve2d`.
    my_nb_trimmed_curve2d: i32,
    /// OCCT `myNbTrimmedCurve3d`.
    my_nb_trimmed_curve3d: i32,
    /// OCCT `myNbBSplibeSurf` (the OCCT spelling).
    my_nb_bsplibe_surf: i32,
    /// OCCT `myNbBezierSurf`.
    my_nb_bezier_surf: i32,
    /// OCCT `myNbTrimSurf`.
    my_nb_trim_surf: i32,
    /// OCCT `myNbWireWitnSeam` (the OCCT spelling).
    my_nb_wire_witn_seam: i32,
    /// OCCT `myNbWireWithSevSeams`.
    my_nb_wire_with_sev_seams: i32,
    /// OCCT `myNbFaceWithSevWires`.
    my_nb_face_with_sev_wires: i32,
    /// OCCT `myNbNoPCurve`.
    my_nb_no_pcurve: i32,
    /// OCCT `myNbFreeFaces`.
    my_nb_free_faces: i32,
    /// OCCT `myNbFreeWires`.
    my_nb_free_wires: i32,
    /// OCCT `myNbFreeEdges`.
    my_nb_free_edges: i32,
    /// OCCT `myNbSharedSolids`.
    my_nb_shared_solids: i32,
    /// OCCT `myNbSharedShells`.
    my_nb_shared_shells: i32,
    /// OCCT `myNbSharedFaces`.
    my_nb_shared_faces: i32,
    /// OCCT `myNbSharedWires`.
    my_nb_shared_wires: i32,
    /// OCCT `myNbSharedFreeWires`.
    my_nb_shared_free_wires: i32,
    /// OCCT `myNbSharedFreeEdges`.
    my_nb_shared_free_edges: i32,
    /// OCCT `myNbSharedEdges`.
    my_nb_shared_edges: i32,
    /// OCCT `myNbSharedVertices`.
    my_nb_shared_vertices: i32,

    // The sections (OCCT hxx L211-216).
    /// OCCT `myBigSplineSec`.
    my_big_spline_sec: Vec<Shape>,
    /// OCCT `myIndirectSec`.
    my_indirect_sec: Vec<Shape>,
    /// OCCT `myOffsetSurfaceSec`.
    my_offset_surface_sec: Vec<Shape>,
    /// OCCT `myTrimmed3dSec`.
    my_trimmed3d_sec: Vec<Shape>,
    /// OCCT `myOffsetCurveSec`.
    my_offset_curve_sec: Vec<Shape>,
    /// OCCT `myTrimmed2dSec`.
    my_trimmed2d_sec: Vec<Shape>,

    // The flags (OCCT hxx L168-206).
    /// OCCT `myBigSplineMode`.
    my_big_spline_mode: bool,
    /// OCCT `myIndirectMode`.
    my_indirect_mode: bool,
    /// OCCT `myOffsetSurfaceMode`.
    my_offset_surface_mode: bool,
    /// OCCT `myTrimmed3dMode`.
    my_trimmed3d_mode: bool,
    /// OCCT `myOffsetCurveMode`.
    my_offset_curve_mode: bool,
    /// OCCT `myTrimmed2dMode`.
    my_trimmed2d_mode: bool,
}

impl Default for ShapeAnalysisShapeContents {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeAnalysisShapeContents {
    /// OCCT ShapeAnalysis_ShapeContents() (cxx L46-55).
    pub fn new() -> Self {
        let mut this = ShapeAnalysisShapeContents {
            my_nb_solids: 0,
            my_nb_shells: 0,
            my_nb_faces: 0,
            my_nb_wires: 0,
            my_nb_edges: 0,
            my_nb_vertices: 0,
            my_nb_solids_with_voids: 0,
            my_nb_big_splines: 0,
            my_nb_c0_surfaces: 0,
            my_nb_c0_curves: 0,
            my_nb_offset_surf: 0,
            my_nb_indirect_surf: 0,
            my_nb_offset_curves: 0,
            my_nb_trimmed_curve2d: 0,
            my_nb_trimmed_curve3d: 0,
            my_nb_bsplibe_surf: 0,
            my_nb_bezier_surf: 0,
            my_nb_trim_surf: 0,
            my_nb_wire_witn_seam: 0,
            my_nb_wire_with_sev_seams: 0,
            my_nb_face_with_sev_wires: 0,
            my_nb_no_pcurve: 0,
            my_nb_free_faces: 0,
            my_nb_free_wires: 0,
            my_nb_free_edges: 0,
            my_nb_shared_solids: 0,
            my_nb_shared_shells: 0,
            my_nb_shared_faces: 0,
            my_nb_shared_wires: 0,
            my_nb_shared_free_wires: 0,
            my_nb_shared_free_edges: 0,
            my_nb_shared_edges: 0,
            my_nb_shared_vertices: 0,
            my_big_spline_sec: Vec::new(),
            my_indirect_sec: Vec::new(),
            my_offset_surface_sec: Vec::new(),
            my_trimmed3d_sec: Vec::new(),
            my_offset_curve_sec: Vec::new(),
            my_trimmed2d_sec: Vec::new(),
            my_big_spline_mode: false,
            my_indirect_mode: false,
            my_offset_surface_mode: false,
            my_trimmed3d_mode: false,
            my_offset_curve_mode: false,
            my_trimmed2d_mode: false,
        };
        this.clear_flags();
        this
    }

    /// OCCT Clear() (cxx L57-100).
    pub fn clear(&mut self) {
        self.my_nb_solids = 0;
        self.my_nb_shells = 0;
        self.my_nb_faces = 0;
        self.my_nb_wires = 0;
        self.my_nb_edges = 0;
        self.my_nb_vertices = 0;
        self.my_nb_solids_with_voids = 0;
        self.my_nb_big_splines = 0;
        self.my_nb_c0_surfaces = 0;
        self.my_nb_c0_curves = 0;
        self.my_nb_offset_surf = 0;
        self.my_nb_indirect_surf = 0;
        self.my_nb_offset_curves = 0;
        self.my_nb_trimmed_curve2d = 0;
        self.my_nb_trimmed_curve3d = 0;
        self.my_nb_bsplibe_surf = 0;
        self.my_nb_bezier_surf = 0;
        self.my_nb_trim_surf = 0;
        self.my_nb_wire_witn_seam = 0;
        self.my_nb_wire_with_sev_seams = 0;
        self.my_nb_face_with_sev_wires = 0;
        self.my_nb_no_pcurve = 0;
        self.my_nb_free_faces = 0;
        self.my_nb_free_wires = 0;
        self.my_nb_free_edges = 0;

        self.my_nb_shared_solids = 0;
        self.my_nb_shared_shells = 0;
        self.my_nb_shared_faces = 0;
        self.my_nb_shared_wires = 0;
        self.my_nb_shared_free_wires = 0;
        self.my_nb_shared_free_edges = 0;
        self.my_nb_shared_edges = 0;
        self.my_nb_shared_vertices = 0;

        self.my_big_spline_sec.clear();
        self.my_indirect_sec.clear();
        self.my_offset_surface_sec.clear();
        self.my_trimmed3d_sec.clear();
        self.my_offset_curve_sec.clear();
        self.my_trimmed2d_sec.clear();
    }

    /// OCCT ClearFlags() (cxx L102-110).
    pub fn clear_flags(&mut self) {
        self.my_big_spline_mode = false;
        self.my_indirect_mode = false;
        self.my_offset_surface_mode = false;
        self.my_trimmed3d_mode = false;
        self.my_offset_curve_mode = false;
        self.my_trimmed2d_mode = false;
    }

    /// OCCT Perform(Shape) (cxx L112-358).
    pub fn perform(&mut self, brep: &mut BRep, shape: &Shape) {
        self.clear();
        //  On y va
        //  On note pour les SOLIDES : ceux qui ont des trous (plus d un SHELL)
        let mut mapsh: HashSet<(u64, u32)> = HashSet::new();

        for exp in topexp_explorer(brep, shape, ShapeType::Solid) {
            // OCCT: TopoDS_Solid sol = TopoDS::Solid(exp.Current());
            // sol.Location(TopLoc_Location()); — the identity-location
            // clone (bridge #4).
            let mut sol = exp;
            sol.location = 0;
            mapsh.insert((std::sync::Arc::as_ptr(&sol.data) as u64, sol.location));
            let mut nbs = 0;
            for _shel in topexp_explorer(brep, &sol, ShapeType::Shell) {
                nbs += 1;
            }
            if nbs > 1 {
                self.my_nb_solids_with_voids += 1;
            }
            self.my_nb_solids += 1;
        }
        self.my_nb_shared_solids = mapsh.len() as i32;

        //  Pour les SHELLS, on compte les faces dans les SHELLS
        //  Ensuite une soustraction, et on a les faces libres
        mapsh.clear();
        let mut nbfaceshell = 0;
        for exp in topexp_explorer(brep, shape, ShapeType::Shell) {
            self.my_nb_shells += 1;
            let mut she = exp;
            she.location = 0;
            mapsh.insert((std::sync::Arc::as_ptr(&she.data) as u64, she.location));
            for _shel in topexp_explorer(brep, &she, ShapeType::Face) {
                nbfaceshell += 1;
            }
        }
        self.my_nb_shared_shells = mapsh.len() as i32;
        //  On note pour les FACES pas mal de choses (surface, topologie)
        //  * Surface BSpline > 8192 poles
        //  * Surface BSpline "OnlyC0" (not yet impl)
        //  * Surface Offset
        //  * Surface Elementaire INDIRECTE
        //  * Presence de COUTURES; en particulier WIRE A PLUS D UNE COUTURE
        //  * Edge : OffsetCurve

        mapsh.clear();
        for exp in topexp_explorer(brep, shape, ShapeType::Face) {
            let face = exp;
            self.my_nb_faces += 1;
            let (mut surf, _loc) = brep_tool_surface_loc(&face);
            let mut face = face;
            face.location = 0;
            mapsh.insert((std::sync::Arc::as_ptr(&face.data) as u64, face.location));
            // OCCT L171-177: the Geom_RectangularTrimmedSurface unwrap.
            if let Some(Surface3::Trimmed(trsu)) = surf.as_ref() {
                self.my_nb_trim_surf += 1;
                surf = Some(trsu.basis.as_ref().clone());
            }
            // #10 rln 27/02/98 BUC50003 entity 56
            // C0 if at least in one direction (U or V)
            let is_cn = match surf.as_ref() {
                None => false,
                Some(s) => match s {
                    Surface3::BSpline(bs) => {
                        bspline_is_cn(&bs.knots_u, bs.degree_u, 1)
                            && bspline_is_cn(&bs.knots_v, bs.degree_v, 1)
                    }
                    // The non-BSpline carriers are CN at all orders (the
                    // OCCT Geom_Surface default).
                    _ => true,
                },
            };
            if surf.is_some() && !is_cn {
                self.my_nb_c0_surfaces += 1;
            }

            if let Some(Surface3::BSpline(bsps)) = surf.as_ref() {
                self.my_nb_bsplibe_surf += 1;
                let nbu = bsps.control_points.len();
                let nbv = bsps.control_points.first().map_or(0, |row| row.len());
                if nbu * nbv > 8192 {
                    self.my_nb_big_splines += 1;
                    if self.my_big_spline_mode {
                        self.my_big_spline_sec.push(face.clone());
                    }
                }
            }
            if let Some(els) = surf.as_ref() {
                if matches!(
                    els,
                    Surface3::Plane(_)
                        | Surface3::Cylinder(_)
                        | Surface3::Cone(_)
                        | Surface3::Sphere(_)
                        | Surface3::Torus(_)
                ) {
                    // OCCT L201: if (!els->Position().Direct()) — the rcad
                    // elementary carriers are right-handed by construction
                    // (bridge #6), so the indirect branch stays off.
                    let direct = true;
                    if !direct {
                        self.my_nb_indirect_surf += 1;
                        if self.my_indirect_mode {
                            self.my_indirect_sec.push(face.clone());
                        }
                    }
                }
            }
            if matches!(surf.as_ref(), Some(Surface3::Offset(_))) {
                self.my_nb_offset_surf += 1;
                if self.my_offset_surface_mode {
                    self.my_offset_surface_sec.push(face.clone());
                }
            } else if matches!(surf.as_ref(), Some(Surface3::Bezier(_))) {
                self.my_nb_bezier_surf += 1;
            }

            let mut maxseam = 0;
            let mut nbwires = 0;
            for wires in topexp_explorer(brep, &face, ShapeType::Wire) {
                let wire = wires;
                let mut nbseam = 0;
                nbwires += 1;
                for edg in topexp_explorer(brep, &wire, ShapeType::Edge) {
                    let edge = edg;
                    if brep_tool_is_closed_edge_face(brep, &edge, &face) {
                        nbseam += 1;
                    }
                    // OCCT L231-232: double first, last (the out range).
                    let c3d = brep_tool_curve(&edge);
                    let (c3d, _range) = match c3d {
                        Some((c, f, l)) => (Some(c), (f, l)),
                        None => (None, (0.0, 0.0)),
                    };
                    if let Some(c3d) = c3d.as_ref() {
                        if matches!(c3d, Curve3::Trimmed(_)) {
                            self.my_nb_trimmed_curve3d += 1;
                            if self.my_trimmed3d_mode {
                                self.my_trimmed3d_sec.push(face.clone());
                            }
                        }
                    }
                    let c2d = brep_tool_curve_on_surface_face(&edge, &face);
                    match c2d {
                        None => {
                            self.my_nb_no_pcurve += 1;
                        }
                        Some((c2d, _, _)) => {
                            if matches!(c2d, Curve2d::Offset(_)) {
                                self.my_nb_offset_curves += 1;
                                if self.my_offset_curve_mode {
                                    self.my_offset_curve_sec.push(face.clone());
                                }
                            } else if matches!(c2d, Curve2d::Trimmed(_)) {
                                self.my_nb_trimmed_curve2d += 1;
                                if self.my_trimmed2d_mode {
                                    self.my_trimmed2d_sec.push(face.clone());
                                }
                            }
                        }
                    }
                }
                if nbseam > maxseam {
                    maxseam = nbseam;
                }
            }
            if maxseam == 1 {
                self.my_nb_wire_witn_seam += 1;
            } else if maxseam > 1 {
                self.my_nb_wire_with_sev_seams += 1;
            }
            if nbwires > 1 {
                self.my_nb_face_with_sev_wires += 1;
            }
        }
        self.my_nb_shared_faces = mapsh.len() as i32;

        mapsh.clear();
        for exp in topexp_explorer(brep, shape, ShapeType::Wire) {
            let mut wire = exp;
            wire.location = 0;
            mapsh.insert((std::sync::Arc::as_ptr(&wire.data) as u64, wire.location));
            self.my_nb_wires += 1;
        }
        self.my_nb_shared_wires = mapsh.len() as i32;

        //  Ne pas oublier les FACES :
        self.my_nb_free_faces = self.my_nb_faces - nbfaceshell;

        mapsh.clear();
        for exp in topexp_explorer(brep, shape, ShapeType::Edge) {
            let mut edge = exp;
            edge.location = 0;
            mapsh.insert((std::sync::Arc::as_ptr(&edge.data) as u64, edge.location));
            self.my_nb_edges += 1;
            let (c3d, _loc, _first, _last) = brep_tool_curve_loc(&edge);
            if let Some(c3d) = c3d.as_ref() {
                if matches!(c3d, Curve3::Offset(_)) {
                    self.my_nb_offset_curves += 1;
                    if self.my_offset_curve_mode {
                        self.my_offset_curve_sec.push(edge.clone());
                    }
                }
                // OCCT L322: !c3d->IsCN(1) — only the BSpline carrier can
                // decline (the Geom_Curve default is CN).
                let is_cn = match c3d {
                    Curve3::BSpline(bs) => {
                        bspline_is_cn(&bs.knots, bs.degree, 1)
                    }
                    _ => true,
                };
                if !is_cn {
                    self.my_nb_c0_curves += 1;
                }
            }
        }
        self.my_nb_shared_edges = mapsh.len() as i32;

        mapsh.clear();
        for exp in topexp_explorer(brep, shape, ShapeType::Vertex) {
            let mut vert = exp;
            vert.location = 0;
            self.my_nb_vertices += 1;
            mapsh.insert((std::sync::Arc::as_ptr(&vert.data) as u64, vert.location));
        }
        self.my_nb_shared_vertices = mapsh.len() as i32;

        mapsh.clear();
        for exp in topexp_explorer_avoid(brep, shape, ShapeType::Edge, ShapeType::Face) {
            let mut edge = exp;
            edge.location = 0;
            self.my_nb_free_edges += 1;
            mapsh.insert((std::sync::Arc::as_ptr(&edge.data) as u64, edge.location));
        }
        self.my_nb_shared_free_edges = mapsh.len() as i32;

        mapsh.clear();
        for exp in topexp_explorer_avoid(brep, shape, ShapeType::Wire, ShapeType::Face) {
            let mut wire = exp;
            wire.location = 0;
            self.my_nb_free_wires += 1;
            mapsh.insert((std::sync::Arc::as_ptr(&wire.data) as u64, wire.location));
        }
        self.my_nb_shared_free_wires = mapsh.len() as i32;
    }

    // -- the hxx inline accessors (L110-206) -------------------------------

    /// OCCT BigSplineSec() (hxx L137).
    pub fn big_spline_sec(&self) -> &[Shape] {
        &self.my_big_spline_sec
    }

    /// OCCT IndirectSec() (hxx L142).
    pub fn indirect_sec(&self) -> &[Shape] {
        &self.my_indirect_sec
    }

    /// OCCT OffsetSurfaceSec() (hxx L147).
    pub fn offset_surface_sec(&self) -> &[Shape] {
        &self.my_offset_surface_sec
    }

    /// OCCT Trimmed3dSec() (hxx L152).
    pub fn trimmed3d_sec(&self) -> &[Shape] {
        &self.my_trimmed3d_sec
    }

    /// OCCT OffsetCurveSec() (hxx L157).
    pub fn offset_curve_sec(&self) -> &[Shape] {
        &self.my_offset_curve_sec
    }

    /// OCCT Trimmed2dSec() (hxx L162).
    pub fn trimmed2d_sec(&self) -> &[Shape] {
        &self.my_trimmed2d_sec
    }

    /// OCCT ModifyBigSplineMode() (hxx L52) — the modifiable flag.
    pub fn modify_big_spline_mode(&mut self) -> &mut bool {
        &mut self.my_big_spline_mode
    }

    /// OCCT ModifyIndirectMode() (hxx L56).
    pub fn modify_indirect_mode(&mut self) -> &mut bool {
        &mut self.my_indirect_mode
    }

    /// OCCT ModifyOffsetSurfaceMode() (hxx L60).
    pub fn modify_offset_surface_mode(&mut self) -> &mut bool {
        &mut self.my_offset_surface_mode
    }

    /// OCCT ModifyTrimmed3dMode() (hxx L64).
    pub fn modify_trimmed3d_mode(&mut self) -> &mut bool {
        &mut self.my_trimmed3d_mode
    }

    /// OCCT ModifyOffsetCurveMode() (hxx L68).
    pub fn modify_offset_curve_mode(&mut self) -> &mut bool {
        &mut self.my_offset_curve_mode
    }

    /// OCCT ModifyTrimmed2dMode() (hxx L72).
    pub fn modify_trimmed2d_mode(&mut self) -> &mut bool {
        &mut self.my_trimmed2d_mode
    }

    /// OCCT NbSolids() (hxx L75).
    pub fn nb_solids(&self) -> i32 {
        self.my_nb_solids
    }

    /// OCCT NbShells() (hxx L77).
    pub fn nb_shells(&self) -> i32 {
        self.my_nb_shells
    }

    /// OCCT NbFaces() (hxx L79).
    pub fn nb_faces(&self) -> i32 {
        self.my_nb_faces
    }

    /// OCCT NbWires() (hxx L81).
    pub fn nb_wires(&self) -> i32 {
        self.my_nb_wires
    }

    /// OCCT NbEdges() (hxx L83).
    pub fn nb_edges(&self) -> i32 {
        self.my_nb_edges
    }

    /// OCCT NbVertices() (hxx L85).
    pub fn nb_vertices(&self) -> i32 {
        self.my_nb_vertices
    }

    /// OCCT NbSolidsWithVoids() (hxx L87).
    pub fn nb_solids_with_voids(&self) -> i32 {
        self.my_nb_solids_with_voids
    }

    /// OCCT NbBigSplines() (hxx L89).
    pub fn nb_big_splines(&self) -> i32 {
        self.my_nb_big_splines
    }

    /// OCCT NbC0Surfaces() (hxx L91).
    pub fn nb_c0_surfaces(&self) -> i32 {
        self.my_nb_c0_surfaces
    }

    /// OCCT NbC0Curves() (hxx L93).
    pub fn nb_c0_curves(&self) -> i32 {
        self.my_nb_c0_curves
    }

    /// OCCT NbOffsetSurf() (hxx L95).
    pub fn nb_offset_surf(&self) -> i32 {
        self.my_nb_offset_surf
    }

    /// OCCT NbIndirectSurf() (hxx L97).
    pub fn nb_indirect_surf(&self) -> i32 {
        self.my_nb_indirect_surf
    }

    /// OCCT NbOffsetCurves() (hxx L99).
    pub fn nb_offset_curves(&self) -> i32 {
        self.my_nb_offset_curves
    }

    /// OCCT NbTrimmedCurve2d() (hxx L101).
    pub fn nb_trimmed_curve2d(&self) -> i32 {
        self.my_nb_trimmed_curve2d
    }

    /// OCCT NbTrimmedCurve3d() (hxx L103).
    pub fn nb_trimmed_curve3d(&self) -> i32 {
        self.my_nb_trimmed_curve3d
    }

    /// OCCT NbBSplibeSurf() (hxx L105).
    pub fn nb_bsplibe_surf(&self) -> i32 {
        self.my_nb_bsplibe_surf
    }

    /// OCCT NbBezierSurf() (hxx L107).
    pub fn nb_bezier_surf(&self) -> i32 {
        self.my_nb_bezier_surf
    }

    /// OCCT NbTrimSurf() (hxx L109).
    pub fn nb_trim_surf(&self) -> i32 {
        self.my_nb_trim_surf
    }

    /// OCCT NbWireWitnSeam() (hxx L111).
    pub fn nb_wire_witn_seam(&self) -> i32 {
        self.my_nb_wire_witn_seam
    }

    /// OCCT NbWireWithSevSeams() (hxx L113).
    pub fn nb_wire_with_sev_seams(&self) -> i32 {
        self.my_nb_wire_with_sev_seams
    }

    /// OCCT NbFaceWithSevWires() (hxx L115).
    pub fn nb_face_with_sev_wires(&self) -> i32 {
        self.my_nb_face_with_sev_wires
    }

    /// OCCT NbNoPCurve() (hxx L117).
    pub fn nb_no_pcurve(&self) -> i32 {
        self.my_nb_no_pcurve
    }

    /// OCCT NbFreeFaces() (hxx L119).
    pub fn nb_free_faces(&self) -> i32 {
        self.my_nb_free_faces
    }

    /// OCCT NbFreeWires() (hxx L121).
    pub fn nb_free_wires(&self) -> i32 {
        self.my_nb_free_wires
    }

    /// OCCT NbFreeEdges() (hxx L123).
    pub fn nb_free_edges(&self) -> i32 {
        self.my_nb_free_edges
    }

    /// OCCT NbSharedSolids() (hxx L125).
    pub fn nb_shared_solids(&self) -> i32 {
        self.my_nb_shared_solids
    }

    /// OCCT NbSharedShells() (hxx L127).
    pub fn nb_shared_shells(&self) -> i32 {
        self.my_nb_shared_shells
    }

    /// OCCT NbSharedFaces() (hxx L129).
    pub fn nb_shared_faces(&self) -> i32 {
        self.my_nb_shared_faces
    }

    /// OCCT NbSharedWires() (hxx L131).
    pub fn nb_shared_wires(&self) -> i32 {
        self.my_nb_shared_wires
    }

    /// OCCT NbSharedFreeWires() (hxx L133).
    pub fn nb_shared_free_wires(&self) -> i32 {
        self.my_nb_shared_free_wires
    }

    /// OCCT NbSharedFreeEdges() (hxx L135).
    pub fn nb_shared_free_edges(&self) -> i32 {
        self.my_nb_shared_free_edges
    }

    /// OCCT NbSharedEdges() (hxx L137).
    pub fn nb_shared_edges(&self) -> i32 {
        self.my_nb_shared_edges
    }

    /// OCCT NbSharedVertices() (hxx L139).
    pub fn nb_shared_vertices(&self) -> i32 {
        self.my_nb_shared_vertices
    }
}
