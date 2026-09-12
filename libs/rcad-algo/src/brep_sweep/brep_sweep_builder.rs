//! OCCT BRepSweep_Builder (TKPrim/BRepSweep) — implements the abstract
//! Builder with the BRep Builder.
//!
//! Sources:
//! - BRepSweep_Builder.hxx L29-65
//! - BRepSweep_Builder.cxx L23-79
//! - BRepSweep_Builder.lxx L19-22 (the Builder() accessor)
//!
//! The BRep_Builder methods consumed through `myBuilder.Builder()` inside
//! BRepSweep_Translation / BRepSweep_Rotation (MakeVertex, MakeEdge,
//! MakeFace, UpdateEdge, UpdateVertex, Degenerated, Continuity) are re-hosted
//! on [`BRepSweepBRepBuilder`] with their OCCT anchors
//! (BRep_Builder.cxx / TopoDS_Builder.cxx).
//!
//! Architecture difference (the pool-free style): OCCT edits TopoDS TShapes
//! in place through BRep_Builder; the rcad TShape payload travels inside the
//! Shape's Arc and the in-place edits go through the kernel
//! edge_mut_inplace / wire_mut unsafe pattern (identity is preserved across
//! the aliased myShapes slots and parent containers).

use rcad_kernel::geom::{Curve2d, Curve3, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{tshape_flags, GeomAbsShape, Orientation, ShapeType, TShape, TVertexData, TEdgeData, TWireData, TFaceData, TShellData, TSolidData, PointRepresentation, CurveRepresentation};
use std::collections::HashMap;
use std::sync::Arc;

use super::tool_rehost::{
    edge_data_mut, face_data_mut, shell_data_mut, solid_data_mut, vertex_data_mut, wire_data_mut,
};

/// OCCT BRepSweep_Builder (BRepSweep_Builder.hxx L29-65).
///
/// Architecture difference: OCCT stores a `BRep_Builder myBuilder` member;
/// the rcad BRep_Builder re-host is a stateless operation set, so the member
/// is carried as a unit struct to keep the member/accessor form.
pub struct BRepSweepBuilder {
    /// OCCT: `BRep_Builder myBuilder` (BRepSweep_Builder.hxx L64).
    pub my_builder: BRepSweepBRepBuilder,
}

impl Default for BRepSweepBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepSweepBuilder {
    /// OCCT BRepSweep_Builder::BRepSweep_Builder(aBuilder) (cxx L23-26).
    pub fn new() -> Self {
        BRepSweepBuilder {
            my_builder: BRepSweepBRepBuilder,
        }
    }

    /// OCCT BRepSweep_Builder::Builder() (lxx L19-22).
    pub fn builder(&self) -> &BRepSweepBRepBuilder {
        &self.my_builder
    }

    /// OCCT BRepSweep_Builder::MakeCompound(aCompound) (cxx L30-33) — returns
    /// an empty Compound in the slot.
    pub fn make_compound(&self, a_compound: &mut Shape) {
        *a_compound = self.my_builder.make_compound();
    }

    /// OCCT BRepSweep_Builder::MakeCompSolid(aCompSolid) (cxx L37-40).
    pub fn make_compsolid(&self, a_compsolid: &mut Shape) {
        *a_compsolid = self.my_builder.make_compsolid();
    }

    /// OCCT BRepSweep_Builder::MakeSolid(aSolid) (cxx L44-47).
    pub fn make_solid(&self, a_solid: &mut Shape) {
        *a_solid = self.my_builder.make_solid();
    }

    /// OCCT BRepSweep_Builder::MakeShell(aShell) (cxx L51-54).
    pub fn make_shell(&self, a_shell: &mut Shape) {
        *a_shell = self.my_builder.make_shell();
    }

    /// OCCT BRepSweep_Builder::MakeWire(aWire) (cxx L58-61).
    pub fn make_wire(&self, a_wire: &mut Shape) {
        *a_wire = self.my_builder.make_wire();
    }

    /// OCCT BRepSweep_Builder::Add(aShape1, aShape2, Orient) (cxx L65-72) —
    /// adds the Shape 2 in the Shape 1, set to <Orient> orientation.
    pub fn add_oriented(&self, a_shape1: &Shape, a_shape2: &Shape, orient: Orientation) {
        let mut a_comp = a_shape2.clone();
        a_comp.orientation = orient;
        self.add(a_shape1, &a_comp);
    }

    /// OCCT BRepSweep_Builder::Add(aShape1, aShape2) (cxx L76-79) —
    /// myBuilder.Add(aShape1, aShape2) (the TopoDS_Builder::Add walk, see
    /// [`BRepSweepBRepBuilder::add`]).
    pub fn add(&self, a_shape1: &Shape, a_shape2: &Shape) {
        self.my_builder.add(a_shape1, a_shape2);
    }
}

// ---------------------------------------------------------------------------
// BRepSweepBRepBuilder — the BRep_Builder / TopoDS_Builder re-host consumed
// through BRepSweep_Builder::Builder().
// ---------------------------------------------------------------------------

/// OCCT BRep_Builder (the operations consumed by the sweep).
pub struct BRepSweepBRepBuilder;

impl BRepSweepBRepBuilder {
    /// OCCT BRep_Builder::MakeCompound(C) — an empty Compound TShape.
    pub fn make_compound(&self) -> Shape {
        Shape {
            data: Arc::new(TShape::Compound(Vec::new())),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeCompSolid(C) — an empty CompSolid TShape.
    pub fn make_compsolid(&self) -> Shape {
        Shape {
            data: Arc::new(TShape::CompSolid(Vec::new())),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeSolid(S) — an empty Solid TShape.
    pub fn make_solid(&self) -> Shape {
        Shape {
            data: Arc::new(TShape::Solid(TSolidData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                shells: Vec::new(),
                internal_vertices: Vec::new(),
                internal_edges: Vec::new(),
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeShell(S) — an empty Shell TShape.
    pub fn make_shell(&self) -> Shape {
        Shape {
            data: Arc::new(TShape::Shell(TShellData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                faces: Vec::new(),
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeWire(W) — an empty Wire TShape.
    pub fn make_wire(&self) -> Shape {
        Shape {
            data: Arc::new(TShape::Wire(TWireData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                edges: Vec::new(),
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeVertex(V, P, Tol) (BRep_Builder.hxx MakeVertex
    /// forms; the point/tolerance land through UpdateVertex(V, P, Tol)
    /// L1203-1213 — the vertex carries the point and the tolerance).
    pub fn make_vertex(&self, a_pnt: glam::DVec3, tol: f64) -> Shape {
        Shape {
            data: Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                point: a_pnt,
                tolerance: tol,
                points: Vec::new(),
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeEdge(E) (BRep_Builder.cxx L623-630) — an empty
    /// Edge TShape.
    pub fn make_edge(&self) -> Shape {
        Shape {
            data: Arc::new(TShape::Edge(TEdgeData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                curve: None,
                first: Shape::null(),
                last: Shape::null(),
                range: [0.0, 0.0],
                degenerated: false,
                pcurves: indexmap::IndexMap::new(),
                representations: Vec::new(),
                vertex_params: HashMap::new(),
                tolerance: 0.0,
                same_parameter: false,
                same_range: false,
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::MakeEdge(E, C, Tol) = MakeEdge(E) +
    /// UpdateEdge(E, C, Tol) (the UpdateCurves walk BRep_Builder.cxx
    /// L66-102): the 3D curve lands in the representations; the edge
    /// tolerance takes the max.
    pub fn make_edge_curve(&self, the_curve: &Curve3, tol: f64) -> Shape {
        let e = self.make_edge();
        {
            let ed = edge_data_mut(&e);
            ed.curve = Some(the_curve.clone());
            ed.tolerance = ed.tolerance.max(tol);
        }
        e
    }

    /// OCCT BRep_Builder::UpdateEdge(E, Tol) (BRep_Builder.cxx
    /// UpdateEdge(E,Tol)) — the tolerance max.
    pub fn update_edge_tol(&self, the_e: &Shape, the_tol: f64) {
        let ed = edge_data_mut(the_e);
        ed.tolerance = ed.tolerance.max(the_tol);
    }

    /// OCCT BRep_Builder::MakeFace(F, S, Tol) (BRep_Builder.cxx L499-512) —
    /// the face carries the surface and the tolerance (no wire yet).
    pub fn make_face(&self, the_surface: &Surface3, tol: f64) -> Shape {
        Shape {
            data: Arc::new(TShape::Face(TFaceData {
                my_shapes: Vec::new(),
                flags: tshape_flags::DEFAULT,
                surface: Some(the_surface.clone()),
                surface_location: 0,
                outer_wire: Shape::null(),
                inner_wires: Vec::new(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: Vec::new(),
                tolerance: tol,
                natural_restriction: false,
            })),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        }
    }

    /// OCCT BRep_Builder::UpdateEdge(E, C2d, F, Tol) (BRep_Builder.cxx
    /// L692-748): the pcurve representation keyed by the face (the
    /// L.Predivided(E.Location()) form — identity locations in this
    /// pipeline).
    pub fn update_edge_pcurve(&self, the_e: &Shape, the_c2d: &Curve2d, the_f: &Shape, the_tol: f64) {
        let key = crate::brep_algo::tool::shape_key(the_f);
        let ed = edge_data_mut(the_e);
        // OCCT: the pcurve range is the edge's current range (UpdateCurves
        // keeps the first GCurve range; a fresh edge carries the 3D range).
        let (f0, l0) = (ed.range[0], ed.range[1]);
        ed.pcurves
            .insert(key, (the_c2d.clone(), f0, l0));
        ed.representations
            .push(CurveRepresentation::CurveOnSurface {
                face: key,
                pcurve: the_c2d.clone(),
                range: [f0, l0],
            });
        ed.tolerance = ed.tolerance.max(the_tol);
    }

    /// OCCT BRep_Builder::UpdateEdge(E, C1, C2, F, Tol) (BRep_Builder.cxx
    /// L750-812): the seam pair — the two pcurves land in one
    /// BRep_CurveOnClosedSurface representation.
    pub fn update_edge_two_pcurves(
        &self,
        the_e: &Shape,
        the_c1: &Curve2d,
        the_c2: &Curve2d,
        the_f: &Shape,
        the_tol: f64,
    ) {
        let key = crate::brep_algo::tool::shape_key(the_f);
        let ed = edge_data_mut(the_e);
        let (f0, l0) = (ed.range[0], ed.range[1]);
        ed.pcurves
            .insert(key, (the_c1.clone(), f0, l0));
        ed.representations
            .push(CurveRepresentation::CurveOnClosedSurface {
                face: key,
                pcurve1: the_c1.clone(),
                pcurve2: the_c2.clone(),
                range: [f0, l0],
            });
        ed.tolerance = ed.tolerance.max(the_tol);
    }

    /// OCCT BRep_Builder::Degenerated(E, ToDegenerate) — the flag write.
    pub fn degenerated(&self, the_e: &Shape, to_degenerate: bool) {
        edge_data_mut(the_e).degenerated = to_degenerate;
    }

    /// OCCT TopoDS_Builder::Add(aShape, aComponent) (TopoDS_Builder.cxx
    /// L37-100): the Free(false) freeze of the component, the compatibility
    /// table, the relative orientation, and the myShapes append.  The
    /// relative-location branch (aLoc.Inverted() Move) is the identity no-op:
    /// the sweep never locates parent containers (arch. diff. #1).  The
    /// typed per-parent views (wire.edges, shell.faces, ...) follow the
    /// myShapes append (the kernel BRepBuilder precedent).
    pub fn add(&self, a_shape: &Shape, a_component: &Shape) {
        // OCCT L39: aComponent.TShape()->Free(false).
        // (the FREE flag clear; the sweep never re-adds a frozen component to
        // another parent, and the kernel precedent keeps the flag write.)
        {
            let flags = match unsafe { &mut *(Arc::as_ptr(&a_component.data) as *mut TShape) } {
                TShape::Vertex(v) => &mut v.flags,
                TShape::Edge(e) => &mut e.flags,
                TShape::Wire(w) => &mut w.flags,
                TShape::Face(f) => &mut f.flags,
                TShape::Shell(sh) => &mut sh.flags,
                TShape::Solid(so) => &mut so.flags,
                TShape::CompSolid(_) | TShape::Compound(_) => return,
            };
            *flags &= !tshape_flags::FREE;
        }
        // OCCT L44: if (aShape.Free()).
        let parent_free = match a_shape.data.as_ref() {
            TShape::Vertex(v) => v.flags & tshape_flags::FREE != 0,
            TShape::Edge(e) => e.flags & tshape_flags::FREE != 0,
            TShape::Wire(w) => w.flags & tshape_flags::FREE != 0,
            TShape::Face(f) => f.flags & tshape_flags::FREE != 0,
            TShape::Shell(s) => s.flags & tshape_flags::FREE != 0,
            TShape::Solid(s) => s.flags & tshape_flags::FREE != 0,
            TShape::CompSolid(_) | TShape::Compound(_) => true,
        };
        if !parent_free {
            // OCCT: throw TopoDS_FrozenShape("TopoDS_Builder::Add").
            panic!("TopoDS_FrozenShape: TopoDS_Builder::Add");
        }
        // OCCT L54-71: the compatibility table (child type -> parent types).
        let i_c = a_component.shape_type();
        let i_s = a_shape.shape_type();
        if !add_compatible(i_c, i_s) {
            // OCCT: throw TopoDS_UnCompatibleShapes("TopoDS_Builder::Add").
            panic!("TopoDS_UnCompatibleShapes: TopoDS_Builder::Add");
        }
        // OCCT L75-79: aChild = aComponent; the relative orientation (the
        // parent REVERSED reverses the child).
        let mut a_child = a_component.clone();
        if a_shape.orientation == Orientation::Reversed {
            a_child.orientation = Orientation::Reversed.compose(a_child.orientation);
        }
        // OCCT L82-89: the relative location — the identity no-op here (arch.
        // diff. #1: no locations on the sweep containers).
        // OCCT L92-95: the myShapes append + Modified(true).
        match i_s {
            ShapeType::Edge => {
                // The edge branch (kernel BRepBuilder::add_to_edge): the
                // FORWARD child lands in `first`, the REVERSED one in `last`
                // (the OCCT reads them back from the myShapes order; the rcad
                // first/last mirror the front/back).
                let ed = edge_data_mut(a_shape);
                if a_child.orientation == Orientation::Reversed {
                    ed.last = a_child.clone();
                } else {
                    ed.first = a_child.clone();
                }
                ed.my_shapes.push(a_child);
            }
            ShapeType::Wire => {
                let wd = wire_data_mut(a_shape);
                wd.edges.push(a_child.clone());
                wd.my_shapes.push(a_child);
            }
            ShapeType::Face => {
                // OCCT BRep_Builder::Add(F, W) / Add(F, V): the face's child
                // list is type-agnostic (TopoDS_Iterator yields every entry),
                // the first WIRE being the outer wire; rcad's typed slots need
                // the dispatch, otherwise a vertex child would take (or land
                // beside) the wire slot.
                let fd = face_data_mut(a_shape);
                match a_child.shape_type() {
                    ShapeType::Wire => {
                        if fd.outer_wire.is_null() {
                            fd.outer_wire = a_child.clone();
                        } else {
                            fd.inner_wires.push(a_child.clone());
                        }
                    }
                    ShapeType::Vertex => fd.internal_vertices.push(a_child.clone()),
                    _ => {}
                }
                fd.my_shapes.push(a_child);
            }
            ShapeType::Shell => {
                let sd = shell_data_mut(a_shape);
                sd.faces.push(a_child.clone());
                sd.my_shapes.push(a_child);
            }
            ShapeType::Solid => {
                let sd = solid_data_mut(a_shape);
                sd.shells.push(a_child.clone());
                sd.my_shapes.push(a_child);
            }
            ShapeType::CompSolid => {
                match unsafe { &mut *(Arc::as_ptr(&a_shape.data) as *mut TShape) } {
                    TShape::CompSolid(cs) => cs.push(a_child.clone()),
                    _ => unreachable!(),
                }
            }
            ShapeType::Compound => {
                match unsafe { &mut *(Arc::as_ptr(&a_shape.data) as *mut TShape) } {
                    TShape::Compound(cd) => cd.push(a_child.clone()),
                    _ => unreachable!(),
                }
            }
            _ => {}
        }
    }

    /// OCCT BRep_Builder::UpdateVertex(V, Par, E, Tol) (BRep_Builder.cxx
    /// L1220-1313): the vertex is searched among the edge's children; the
    /// matched child's orientation selects the update: FORWARD sets First on
    /// the curve representations, REVERSED sets Last, INTERNAL/unmatched
    /// stores the point-on-curve parameter (the kernel
    /// BRepBuilder::update_vertex_on_edge translation).
    pub fn update_vertex_param_on_edge(
        &self,
        v: &Shape,
        par: f64,
        e: &Shape,
        tol: f64,
    ) {
        // OCCT L1224-1227: the infinite-parameter DomainError
        // (Precision::IsPositiveInfinite / IsNegativeInfinite,
        // Precision.hxx L357-367).
        if rcad_kernel::precision::is_positive_infinite_value(par)
            || rcad_kernel::precision::is_negative_infinite_value(par)
        {
            panic!("Standard_DomainError: BRep_Builder::Infinite parameter");
        }
        // OCCT L1237: L = E.Location().Predivided(V.Location()) — identity.
        // OCCT L1240-1242: ori = INTERNAL until matched.
        let mut ori = Orientation::Internal;
        let children = match e.data.as_ref() {
            TShape::Edge(ed) => ed.my_shapes.clone(),
            _ => Vec::new(),
        };
        let degenerated = match e.data.as_ref() {
            TShape::Edge(ed) => ed.degenerated,
            _ => false,
        };
        // OCCT L1249-1253: the degenerated no-vertices case uses the vertex
        // orientation (RLE, june 94).
        if children.is_empty() && degenerated {
            ori = v.orientation;
        }
        for vcur in &children {
            if v.is_same(vcur) {
                ori = vcur.orientation;
                if ori == v.orientation {
                    break;
                }
            }
        }
        {
            let ed = edge_data_mut(e);
            match ori {
                Orientation::Forward => {
                    // OCCT: GC->First(Par) on every GCurve.
                    ed.range[0] = par;
                    for entry in ed.pcurves.values_mut() {
                        entry.1 = par;
                    }
                }
                Orientation::Reversed => {
                    // OCCT: GC->Last(Par) on every GCurve.
                    ed.range[1] = par;
                    for entry in ed.pcurves.values_mut() {
                        entry.2 = par;
                    }
                }
                _ => {
                    // OCCT: UpdatePoints(lpr, Par, ...) — the point-on-curve
                    // parameter (the kernel precedent stores it on the edge's
                    // vertex_params; the rcad read path
                    // BRep_Tool::Parameter(V, E) reads it back).
                    let id = v.ptr_id();
                    ed.vertex_params.insert(id, par);
                }
            }
        }
        // OCCT L1306-1310: TV->UpdateTolerance(Tol) — the max-update.
        {
            let vd = vertex_data_mut(v);
            vd.tolerance = vd.tolerance.max(tol);
        }
    }

    /// OCCT BRep_Builder::UpdateVertex(V, U, V2, F, Tol) (BRep_Builder.cxx
    /// L1417-1442): the UV point-on-surface lands in the vertex points list
    /// and the tolerance takes the max.  Architecture difference: the rcad
    /// PointOnSurface representation carries a pool index slot that has no
    /// pool-free meaning (stored as 0).
    pub fn update_vertex_uv_on_face(
        &self,
        ve: &Shape,
        u: f64,
        v: f64,
        _f: &Shape,
        tol: f64,
    ) {
        let vd = vertex_data_mut(ve);
        // OCCT: UpdatePoints(lpr, U, V, S, L).
        vd.points.push(PointRepresentation::PointOnSurface {
            face: 0,
            u,
            v,
            tolerance: tol,
        });
        // OCCT: TV->UpdateTolerance(Tol).
        vd.tolerance = vd.tolerance.max(tol);
    }

    /// OCCT BRep_Builder::Continuity(E, F1, F2, C) (BRep_Builder.cxx
    /// L1012-1043): the regularity lands as a BRep_CurveOn2Surfaces
    /// representation on the edge (the UpdateCurves regularity walk
    /// L376-405).
    pub fn continuity(&self, the_e: &Shape, f1: &Shape, f2: &Shape, c: GeomAbsShape) {
        // OCCT L1017-1019: S1 = Surface(F1), S2 = Surface(F2) (missing
        // surfaces store nothing — the OCCT null-handle path).
        let s1 = match f1.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        let s2 = match f2.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        };
        let (Some(s1), Some(s2)) = (s1, s2) else {
            return;
        };
        // OCCT L1037-1038: l1/l2 = L.Predivided(E.Location()) — the identity
        // location form in this pipeline.
        let l1 = f1.location;
        let l2 = f2.location;
        // OCCT L1040: UpdateCurves(TE->ChangeCurves(), S1, S2, l1, l2, C).
        let ed = edge_data_mut(the_e);
        for cr in ed.representations.iter_mut() {
            // OCCT L387: cr->IsRegularity(S1, S2, L1, L2).
            if cr.is_regularity_on(&s1, &s2, l1, l2) {
                // OCCT L397-398: cr->Continuity(C).
                if let CurveRepresentation::CurveOn2Surfaces { continuity, .. } = cr {
                    *continuity = c;
                }
                return;
            }
        }
        // OCCT L402-403: a new BRep_CurveOn2Surfaces(S1, S2, L1, L2, C).
        ed.representations
            .push(CurveRepresentation::CurveOn2Surfaces {
                surface1: s1,
                surface2: s2,
                location1: l1,
                location2: l2,
                continuity: c,
            });
    }
}

/// OCCT TopoDS_Builder::Add compatibility table (TopoDS_Builder.cxx L46-69):
/// (child type, parent type) -> the child can be added.
fn add_compatible(child: ShapeType, parent: ShapeType) -> bool {
    use ShapeType::*;
    match (child, parent) {
        (_, Compound) => true,
        (Solid, CompSolid) => true,
        (Shell, Solid) => true,
        (Face, Shell) => true,
        (Wire, Face) => true,
        (Edge, Wire) | (Edge, Solid) => true,
        (Vertex, Edge) | (Vertex, Face) | (Vertex, Solid) => true,
        _ => false,
    }
}

