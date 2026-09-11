//! OCCT BRepBuilderAPI_Transform (TKTopAlgo) — apply a gp_Trsf to a shape.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepBuilderAPI/BRepBuilderAPI_Transform.cxx
//!         $OCCT_SRC/src/ModelingData/TKBRep/BRepTools/BRepTools_TrsfModification.cxx
//!
//! OCCT BRepBuilderAPI_Transform::Perform (cxx L44-64) splits in two:
//! - `myUseModif` (copyGeom || negative || |scale| != 1): the shape is
//!   rebuilt through BRepTools_Modifier with a BRepTools_TrsfModification —
//!   every vertex point / edge curve / face surface is replaced by its
//!   transformed image (NewPoint cxx L301-309 / NewCurve cxx L275-297 /
//!   NewSurface cxx L65-92), tolerances scale by |scale| and mirrored faces
//!   come out REVERSED (RevFace, NewSurface L82 + Rebuild L264-267/L644).
//! - otherwise (a rigid motion with copyGeom=false): the geometry is NOT
//!   touched — the transform is carried as a TopLoc_Location on the result
//!   wrapper (cxx L58-63, `theShape.Moved(myLocation)`).
//!
//! rcad architecture adaptation: the flat-pool BRep carries per-shape
//! location indices, but the consumer read paths that the transformed shape
//! feeds (the feat BRep_Tool re-hosts, brep_algo/tool.rs arch. difference #1)
//! read TShape geometry flat, so a located-but-untransformed shape would be
//! unobservable.  Both OCCT branches therefore materialise the transformation
//! on the TShape geometry (the observable geometry of OCCT's located shape
//! equals the materialised one); the branch decision still gates the
//! TrsfModification extras (tolerance scaling, RevFace) exactly as OCCT does.
//! The per-TShape geometry pass is `BRep::apply_transform` (the kernel's
//! TrsfModification equivalent — its comments carry the OCCT anchors for the
//! range / vertex-parameter / pcurve reparameterisation steps).

use std::collections::HashSet;
use std::sync::Arc;

use rcad_kernel::math::gp::Trsf;
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{BRep, Orientation, TShape};

/// OCCT TopLoc_Location::ScalePrec() — TopLoc_Location.hxx L157.
const SCALE_PREC: f64 = 1.0e-14;

/// OCCT BRepBuilderAPI_Transform(shape, trsf).Shape() with the default
/// copyGeom=false / copyMesh=false (cxx L32-40) applied to a whole BRep pool
/// (the flat-pool equivalent of performing on the pool's root shape).
pub fn transform_brep(brep: &mut BRep, trsf: &Trsf) {
    // OCCT Perform cxx L48-49.
    let my_use_modif =
        trsf.is_negative() || ((trsf.scale.abs() - 1.0).abs() > SCALE_PREC);

    // Geometry materialisation (both branches — see the module note).
    brep.apply_transform(trsf.to_daffine3());

    // OCCT Perform cxx L50-57 (myUseModif): the TrsfModification extras.
    if my_use_modif {
        trsf_modification_extras(&brep.tshapes, trsf);
    }
    // OCCT Perform cxx L58-63 (else): myLocation = myTrsf; the rigid transform
    // rides on the shape location.  Materialised above (arch. adaptation).
}

/// OCCT BRepBuilderAPI_Transform::Perform(theShape, ...) restricted to the
/// subgraph of `shape` (the TShapes reachable from it are transformed in
/// place, preserving TShape identity — the shared-handle semantics of both
/// the OCCT rebuild and `Moved`).
pub fn perform_shape(shape: &Shape, trsf: &Trsf) {
    // OCCT Perform cxx L48-49.
    let my_use_modif =
        trsf.is_negative() || ((trsf.scale.abs() - 1.0).abs() > SCALE_PREC);

    // Collect the reachable TShape handles (read-only walk over the graph).
    let mut arcs: Vec<Arc<TShape>> = Vec::new();
    let mut visited: HashSet<u64> = HashSet::new();
    collect_subgraph(shape, &mut arcs, &mut visited);

    // Geometry materialisation on the collected handles.  A scratch pool over
    // the same Arcs lets BRep::apply_transform drive its (already
    // OCCT-anchored) per-TShape pass without duplicating it; the Arc contents
    // are mutated in place, so the original shapes observe the transform.
    let mut scratch = BRep::new();
    scratch.tshapes = arcs;
    scratch.apply_transform(trsf.to_daffine3());
    let arcs = scratch.tshapes;

    // OCCT Perform cxx L50-57 (myUseModif): the TrsfModification extras.
    if my_use_modif {
        trsf_modification_extras(&arcs, trsf);
    }
}

/// BRepTools_TrsfModification extras over a set of TShape handles:
/// - tolerance scaling by |scale| (NewPoint cxx L304-305, NewCurve
///   cxx L283-284, NewSurface cxx L79-80 — an unconditional `Tol *=`),
/// - RevFace: mirrored faces come out REVERSED (NewSurface cxx L82 sets
///   RevFace = IsNegative; Rebuild cxx L264-267 picks ResOr = REVERSED and
///   L644 assigns it — RevWires stays false (cxx L81), so only the face
///   wrapper orientation is assigned, never the children's).
fn trsf_modification_extras(arcs: &[Arc<TShape>], trsf: &Trsf) {
    let abs_scale = trsf.scale.abs();
    let rev_face = trsf.is_negative();
    for arc in arcs {
        // Shared-handle semantics: mutate in place (the clone_arguments_private
        // model — OCCT TShape is a shared handle; single-threaded).
        // SAFETY: no other &TShape for this Arc is alive inside the loop.
        let ptr = Arc::as_ptr(arc) as *mut TShape;
        let ts = unsafe { &mut *ptr };
        match ts {
            TShape::Vertex(vd) => {
                vd.tolerance *= abs_scale;
            }
            TShape::Edge(ed) => {
                ed.tolerance *= abs_scale;
            }
            TShape::Face(fd) => {
                fd.tolerance *= abs_scale;
            }
            TShape::Shell(sd) => {
                if rev_face {
                    for f in &mut sd.faces {
                        f.orientation = Orientation::Reversed;
                    }
                }
            }
            TShape::Compound(children) | TShape::CompSolid(children) => {
                if rev_face {
                    for c in &mut *children {
                        if c.shape_type() == rcad_kernel::topo::topods::ShapeType::Face {
                            c.orientation = Orientation::Reversed;
                        }
                    }
                }
            }
            _ => {}
        }
    }
}

/// Read-only DFS over the shape graph collecting the reachable TShape
/// handles (TShape identity by Arc pointer, OCCT TopTools_ShapeMapHasher
/// semantics).  Container children and the edge vertex links are followed.
pub(crate) fn collect_subgraph(
    s: &Shape,
    arcs: &mut Vec<Arc<TShape>>,
    visited: &mut HashSet<u64>,
) {
    if s.is_null() || !visited.insert(s.ptr_id()) {
        return;
    }
    arcs.push(s.data.clone());
    match s.data.as_ref() {
        TShape::Vertex(_) => {}
        TShape::Edge(ed) => {
            collect_subgraph(&ed.first, arcs, visited);
            collect_subgraph(&ed.last, arcs, visited);
        }
        TShape::Wire(wd) => {
            for e in &wd.edges {
                collect_subgraph(e, arcs, visited);
            }
        }
        TShape::Face(fd) => {
            collect_subgraph(&fd.outer_wire, arcs, visited);
            for w in &fd.inner_wires {
                collect_subgraph(w, arcs, visited);
            }
            for v in &fd.internal_vertices {
                collect_subgraph(v, arcs, visited);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                collect_subgraph(f, arcs, visited);
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                collect_subgraph(sh, arcs, visited);
            }
            for e in &sd.internal_edges {
                collect_subgraph(e, arcs, visited);
            }
            for v in &sd.internal_vertices {
                collect_subgraph(v, arcs, visited);
            }
        }
        TShape::CompSolid(children) | TShape::Compound(children) => {
            for c in children {
                collect_subgraph(c, arcs, visited);
            }
        }
    }
}
