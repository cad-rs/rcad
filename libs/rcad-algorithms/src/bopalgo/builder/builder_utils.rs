use crate::bopalgo::builder::SourceSide;
use crate::bopalgo::builder::types::{BooleanOpType, WireEdgeSource, WireOrientation, WireSegment};
use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::history::{
    BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use crate::inttools::context::Context;
use crate::tolerance::*;
use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods;
use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::bopalgo::builder::angle_2d::angle_2d;
use crate::bopalgo::builder::edge_builders::is_split_to_reverse;
use crate::bopalgo::builder::wire_splitter::{
    are_verts_coincident, edge_angle_2d, edge_uv_tangent, is_edge_isoline, world_to_uv,
};

///  compare two Curve3 for identity (same TShape).
pub(crate) fn curve_eq(a: &Curve3, b: &Curve3) -> bool {
    match (a, b) {
        (Curve3::Circle(ca), Curve3::Circle(cb)) => {
            (ca.center - cb.center).length_squared() < TOLERANCE_ABS_SQ
                && (ca.normal - cb.normal).length_squared() < TOLERANCE_ABS_SQ
                && (ca.radius - cb.radius).abs() < TOLERANCE_ABS
        }
        (Curve3::Line(la), Curve3::Line(lb)) => {
            (la.origin - lb.origin).length_squared() < TOLERANCE_ABS_SQ
                && (la.direction - lb.direction).length_squared() < TOLERANCE_ABS_SQ
        }
        _ => false,
    }
}

pub(crate) fn hash_point(p: DVec3) -> u64 {
    // Quantize to tolerance grid for spatial hashing
    let scale = 1.0 / TOLERANCE_ABS;
    let ix = (p.x * scale).round() as i64;
    let iy = (p.y * scale).round() as i64;
    let iz = (p.z * scale).round() as i64;
    // FNV-1a style hash
    let mut h: u64 = 14695981039346656037;
    for v in [ix, iy, iz] {
        h ^= v as u64;
        h = h.wrapping_mul(1099511628211);
    }
    h
}

/// Annotate a `BooleanHistory` with per-edge and per-vertex origins by
/// matching result BRep positions against the DS vertex/edge pool.
///
/// Both `edge_origins` and `vertex_origins` are filled in-place.
/// OCCT PostTreat equivalent: builds shape-to-origin maps for history tracking.
///
/// OCCT ref: BOPAlgo_Builder_3.cxx  ?`BOPAlgo_Builder::PostTreat`
/// (L1-250: builds `myLocModified` and `myLocGenerated` maps from DS images).
///
/// OCCT PostTreat algorithm (line-by-line mapping):
/// L20-40:  For each original shape, iterate sub-shapes (vertices, edges, faces).
/// L42-80:  Check `myImages[ei]` on each edge  ?if non-empty, record as Modified.
/// L82-110: For edges without images but present in result  ?record as Preserved.
/// L112-130: Generated edges (intersection edges)  ?record in myGenerated.
/// L132-170: For faces, check if wire edges were split  ?Modified; if not in
/// result  ?IsDeleted.
/// L172-200: Generated faces  ?myGenerated.
/// L202-230: Vertex tracking (fromA/fromB/intersection).
/// L232-250: Compute IsDeleted for entities absent from the result shape.
///
/// Differences from OCCT PostTreat:
/// - OCCT's PostTreat builds two maps: *myLocModified* (original -> last-modified
/// shape, for tracking splits and merges) and *myLocGenerated* (original -> list of
/// generated sub-shapes).  rcad's `annotate_history_from_ds` builds a simpler
/// `BooleanHistory` with flat `VertexOrigin`/`EdgeOrigin` arrays indexed by result
/// BRep position.
/// - OCCT PostTreat processes vertices, edges, and faces by iterating the DS images
/// (`myImages`, `myOrigins`, `myShapesSD`) and copying images from the source DS.
/// rcad uses spatial proximity (vertex point comparison) to match result vertices
/// to DS vertices, then traces edge origin from matched endpoints.
/// - OCCT PostTreat sets `myModified` for faces that were split (maps old -> new faces
/// via `myImages`).  rcad builds `FaceOrigin` separately (in `aggregate_face_origin`).
/// - OCCT PostTreat is called once at the end of `BOPAlgo_Builder::Build`.  rcad calls
/// `annotate_history_from_ds` inside `boolean_op_with_retry` after result assembly.
///
/// See also `BooleanHistory::update_with_post_treat()` for a more
/// implementation that uses `ds.my_images` instead of spatial proximity.
///
///  core concept (history tracking from DS) matches OCCT's
/// image-map-based approach, adapted for rcad's flat-array data model.
///  TopExp::MapShapes(myShape, myMapShape)  ?build result S index map.
/// OCCT maps TopoDS_Shape  ?identity for myMapShape lookup.
/// rcad: maps result vertex index  ?DS vertex index, result edge index  ?(DS vertices).
/// Used by PrepareHistory to determine Modified/Generated/Deleted provenance.
#[allow(dead_code)]
pub(crate) fn map_result_shapes(brep: &topods::BRep, ds: &DS) -> (Vec<usize>, Vec<(usize, usize)>) {
    // Collect flat vertex list from topods in ShapeRef.index order
    let topo_vertices: Vec<DVec3> = brep
        .tshapes
        .iter()
        .filter_map(|ts| match &**ts {
            topods::TShape::Vertex(v) => Some(v.point),
            _ => None,
        })
        .collect();
    let mut result_to_ds: Vec<usize> = vec![usize::MAX; topo_vertices.len()];
    for (ri, pt) in topo_vertices.iter().enumerate() {
        for (di, dv) in ds.vertices.iter().enumerate() {
            if (dv.point - *pt).length_squared()
                < crate::tolerance::TOLERANCE_ABS * crate::tolerance::TOLERANCE_ABS * 4.0
            {
                result_to_ds[ri] = di;
                break;
            }
        }
    }
    // Edge pairs from topods edges: map ShapeRef.index -> flat position
    let topo_edges: Vec<(usize, usize)> = brep
        .tshapes
        .iter()
        .filter_map(|ts| match &**ts {
            topods::TShape::Edge(e) => Some((e.first.index, e.last.index)),
            _ => None,
        })
        .collect();
    let edge_pairs: Vec<(usize, usize)> = topo_edges
        .iter()
        .map(|&(s, e)| {
            let ds_s = result_to_ds.get(s).copied().unwrap_or(usize::MAX);
            let ds_e = result_to_ds.get(e).copied().unwrap_or(usize::MAX);
            (ds_s, ds_e)
        })
        .collect();
    (result_to_ds, edge_pairs)
}

///  PrepareHistory (Builder_4.cxx L164-252).
/// OCCT iterates source shapes  ?LocModified  ?AddModified / AddGenerated / Remove.
/// rcad: uses pre-built result_to_ds map to annotate vertex/edge provenance.
#[allow(dead_code)]
pub(crate) fn annotate_history_from_ds(brep: &topods::BRep, history: &mut BooleanHistory, ds: &DS) {
    let (result_to_ds, _) = map_result_shapes(brep, ds);

    let a_vc = ds.a_vertex_count;
    let n_result_verts = brep
        .tshapes
        .iter()
        .filter(|ts| std::matches!(ts.as_ref(), topods::TShape::Vertex(_)))
        .count();
    let mut vertex_origins: Vec<VertexOrigin> = Vec::with_capacity(n_result_verts);
    for ri in 0..n_result_verts {
        let di = result_to_ds[ri];
        let origin = if di == usize::MAX {
            VertexOrigin::Intersection
        } else if di < a_vc {
            VertexOrigin::FromA(di)
        } else {
            VertexOrigin::FromB(di - a_vc)
        };
        vertex_origins.push(origin);
    }
    history.vertex_origins = vertex_origins;

    // --- edge origins ---
    let a_vc = ds.a_vertex_count;
    let n_result_edges = brep
        .tshapes
        .iter()
        .filter(|ts| std::matches!(ts.as_ref(), topods::TShape::Edge(_)))
        .count();
    let mut edge_origins: Vec<EdgeOrigin> = Vec::with_capacity(n_result_edges);
    let a_ec = ds.a_edge_count;
    let total_ds_edges = ds.edge_count();

    for ts in &brep.tshapes {
        if let topods::TShape::Edge(ed) = &**ts {
            let ds_s = result_to_ds
                .get(ed.first.index)
                .copied()
                .unwrap_or(usize::MAX);
            let ds_e = result_to_ds
                .get(ed.last.index)
                .copied()
                .unwrap_or(usize::MAX);

            let origin = if ds_s == usize::MAX || ds_e == usize::MAX {
                EdgeOrigin::Generated
            } else if ds_s < a_vc && ds_e < a_vc {
                // Both endpoints are A vertices =look for a DS edge in A range.
                let found = (0..a_ec.min(total_ds_edges)).find(|&dei| {
                    (ds.edge_start_vertex_ds(dei) == ds_s && ds.edge_end_vertex_ds(dei) == ds_e)
                        || (ds.edge_start_vertex_ds(dei) == ds_e && ds.edge_end_vertex_ds(dei) == ds_s)
                });
                match found {
                    Some(dei) => EdgeOrigin::FromA(dei),
                    None => EdgeOrigin::SplitFromA(ds_s.min(a_vc - 1)),
                }
            } else if ds_s >= a_vc && ds_e >= a_vc {
                // Both endpoints are B vertices =look for a DS edge in B range.
                let found = (a_ec..total_ds_edges).find(|&dei| {
                    (ds.edge_start_vertex_ds(dei) == ds_s && ds.edge_end_vertex_ds(dei) == ds_e)
                        || (ds.edge_start_vertex_ds(dei) == ds_e && ds.edge_end_vertex_ds(dei) == ds_s)
                });
                match found {
                    Some(dei) => EdgeOrigin::FromB(dei - a_ec),
                    None => {
                        EdgeOrigin::SplitFromB(ds_s.min(ds.vertex_count().saturating_sub(1)) - a_vc)
                    }
                }
            } else {
                EdgeOrigin::Generated
            };
            edge_origins.push(origin);
        }
    }
    history.edge_origins = edge_origins;
}

pub(crate) fn aggregate_face_region_origin(face_origins: &[FaceOrigin]) -> ShellOrigin {
    let mut has_a = false;
    let mut has_b = false;
    let mut has_generated = false;
    for origin in face_origins {
        match origin {
            FaceOrigin::FromA(_) => has_a = true,
            FaceOrigin::FromB(_) => has_b = true,
            FaceOrigin::Generated => has_generated = true,
        }
    }

    match (has_a, has_b, has_generated) {
        (true, false, false) => ShellOrigin::FromA,
        (false, true, false) => ShellOrigin::FromB,
        (false, false, true) => ShellOrigin::Generated,
        _ => ShellOrigin::Mixed,
    }
}

pub(crate) fn aggregate_shell_region_origin(shell_origins: &[ShellOrigin]) -> SolidOrigin {
    let mut has_a = false;
    let mut has_b = false;
    let mut has_generated = false;
    let mut has_mixed = false;
    for origin in shell_origins {
        match origin {
            ShellOrigin::FromA => has_a = true,
            ShellOrigin::FromB => has_b = true,
            ShellOrigin::Generated => has_generated = true,
            ShellOrigin::Mixed => has_mixed = true,
        }
    }

    if has_mixed {
        return SolidOrigin::Mixed;
    }

    match (has_a, has_b, has_generated) {
        (true, false, false) => SolidOrigin::FromA,
        (false, true, false) => SolidOrigin::FromB,
        (false, false, true) => SolidOrigin::Generated,
        _ => SolidOrigin::Mixed,
    }
}

///  PrepareHistory shell/solid provenance (Builder_4.cxx L164-252).
/// OCCT iterates source shapes  ?LocModified  ?AddModified/AddGenerated/Remove.
/// rcad: aggregates per-face origins to shell/solid level via face_region  ?shell  ?solid.
pub(crate) fn annotate_shell_and_solid_history(brep: &topods::BRep, history: &mut BooleanHistory) {
    let mut face_cursor = 0;
    let mut shell_origins = Vec::new();
    let mut solid_origins = Vec::new();

    for ts in &brep.tshapes {
        if let topods::TShape::Solid(sd) = &**ts {
            let solid_shell_start = shell_origins.len();
            for shell_sr in &sd.shells {
                if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                    let shell_face_count = shd.faces.len();
                    let shell_face_origins = history
                        .face_origins
                        .get(face_cursor..face_cursor + shell_face_count)
                        .unwrap_or(&[]);
                    shell_origins.push(aggregate_face_region_origin(shell_face_origins));
                    face_cursor += shell_face_count;
                }
            }
            solid_origins.push(aggregate_shell_region_origin(
                &shell_origins[solid_shell_start..],
            ));
        }
    }

    if face_cursor != history.face_origins.len() {
        // Face count mismatch: BRep has more/fewer faces than history tracks.
        // This happens when compound reconstruction adds/removes faces or when
        // the face order in BRep differs from the emission order.  OCCT's
        // history tracking works with TopoDS shape identity  ?rcad's index-based
        // tracking is inherently more fragile.  Pad shell_origins to match.
        eprintln!(
            "[HISTORY] face_cursor={} != history={}",
            face_cursor,
            history.face_origins.len()
        );
    }
    history.shell_origins = shell_origins;
    history.solid_origins = solid_origins;
}

// =============================================================================
// OCCT 1:1  ? IsInternalFace (BOPTools_AlgoTools.cxx L791-872)
// =============================================================================

///  ?MEF (Map Edge= aces) = 椤?椤?閳??
/// OCCT BOPAlgo_FillIn3DParts::MapEdgesAndFaces (BOPAlgo_Tools.cxx L1479-1503)
/// IsTangentFace (BOPTools_AlgoTools).
/// Checks if two faces are tangent (parallel normals + close distance).
pub fn is_tangent_face(
    fi_a: usize,
    fi_b: usize,
    ds: &crate::bopds::ds::DS,
    angle_tol: f64,
    dist_tol: f64,
) -> bool {
    let n_dot = ds.face_normal(fi_a).dot(ds.face_normal(fi_b)).abs();
    if n_dot < angle_tol.cos() {
        return false;
    }
    let sample_a = if !ds.face_boundary_verts(fi_a).is_empty() {
        ds.vertex_point(ds.face_boundary_verts(fi_a)[0])
    } else {
        return false;
    };
    let dist = match ds.face_surface(fi_b).unwrap() {
        rcad_kernel::geom::Surface3::Plane(p) => (sample_a - p.origin).dot(p.normal).abs(),
        rcad_kernel::geom::Surface3::Sphere(s) => ((sample_a - s.center).length() - s.radius).abs(),
        _ => return false,
    };
    dist < dist_tol
}

pub(crate) fn build_edge_bounds(
    face_indices: &[usize],
    ds: &DS,
) -> std::collections::BTreeSet<usize> {
    let mut bounds: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for &fi in face_indices {
        for &ei in ds.face_boundary_edges(fi) {
            bounds.insert(ei);
        }
    }
    bounds
}

///  PointInFace  椤??= ?FaceSampleData =UV domain  = ?
/// OCCT BOPTools_AlgoTools3D.cxx L885-917
///
/// rcad  ? FaceSampleData  ?uv_domain =uv_centroid,= ?UV centroid
///  閿?=(OCCT =Hatcher =2D point-in-face, ?rcad =FaceSampleData
///  椤?椤? ?UV centroid = ? ?
// (point_in_face, classify_by_off_solid_edge removed  ?dead after ComputeState alignment)

/// = =3D  閿??u64 key,= 椤掑倵鍋??
pub(crate) fn quantize_pos(p: DVec3, tolerance: f64) -> u64 {
    let scale = 1.0 / tolerance;
    let x = (p.x * scale).round() as i64;
    let y = (p.y * scale).round() as i64;
    let z = (p.z * scale).round() as i64;
    //  = ?u64
    let xb = (x as u64) & 0x3FFFFF;
    let yb = (y as u64) & 0x3FFFFF;
    let zb = (z as u64) & 0x3FFFFF;
    (xb << 42) | (yb << 21) | zb
}

///  IsInternalFace  椤??(BOPTools_AlgoTools.cxx L791-872)
///
///  椤?
/// Level 1:  椤?= 閳?=solid  閿?椤? 1  = 閵??
/// 椤? 閳╁唭? ?solid = ?
/// Level 2: ComputeState == =  solid  閿?= 椤愩儺鍞??
/// = ?PointInFace =classify_point ?
///
///  閳?? Some(true) =  ?solid = ?(IN)
/// Some(false) =  ?solid = ?(OUT)
/// None =  
/// Check if a DS vertex lies on the boundary edge between sv/ev, and if so add it
/// to split_verts with its parametric position t.
///  FillImagesEdges checks pave blocks per edge (global scope).
pub(crate) fn check_and_add_split_vertex(
    ds: &DS,
    sv: usize,
    ev: usize,
    vi: usize,
    p_a: DVec3,
    ab: DVec3,
    ab_len2: f64,
    split_verts: &mut Vec<(usize, f64)>,
) {
    if vi == sv || vi == ev {
        return;
    }
    let p = ds.vertex_point(vi);
    let ap = p - p_a;
    let t = ap.dot(ab) / ab_len2;
    if t > 1e-8 && t < 1.0 - 1e-8 {
        let proj = p_a + ab * t;
        if (p - proj).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT {
            split_verts.push((vi, t));
        }
    }
}

pub(crate) fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}
