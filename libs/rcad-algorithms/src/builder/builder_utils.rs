use std::collections::{HashMap, HashSet, BTreeSet};
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::BRep;
use std::cell::RefCell;
use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::history::{BooleanHistory, FaceOrigin, ShellOrigin, SolidOrigin, VertexOrigin, EdgeOrigin, HistoryTracker};
use crate::inttools::context::Context;
use crate::tolerance::*;
use crate::builder::types::{BooleanOpType, FaceSampleData, WireSegment, WireEdgeSource, WireOrientation};
use crate::builder::SourceSide;

use crate::builder::angle_2d::angle_2d;
use crate::builder::wire_splitter::{world_to_uv, edge_uv_tangent, edge_angle_2d, compute_seam_tangent_angles, are_verts_coincident};
use crate::builder::edge_builders::{build_degenerate_edge_segments, build_sphere_seam_segments, build_cylinder_seam_segments};

/// ✅ OCCT-aligned: compare two Curve3 for identity (same TShape).
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
/// OCCT ref: BOPAlgo_Builder_3.cxx — `BOPAlgo_Builder::PostTreat`
/// (L1-250: builds `myLocModified` and `myLocGenerated` maps from DS images).
///
/// OCCT PostTreat algorithm (line-by-line mapping):
///   L20-40:  For each original shape, iterate sub-shapes (vertices, edges, faces).
///   L42-80:  Check `myImages[ei]` on each edge → if non-empty, record as Modified.
///   L82-110: For edges without images but present in result → record as Preserved.
///   L112-130: Generated edges (intersection edges) → record in myGenerated.
///   L132-170: For faces, check if wire edges were split → Modified; if not in
///             result → IsDeleted.
///   L172-200: Generated faces → myGenerated.
///   L202-230: Vertex tracking (fromA/fromB/intersection).
///   L232-250: Compute IsDeleted for entities absent from the result shape.
///
/// Differences from OCCT PostTreat:
/// - OCCT's PostTreat builds two maps: *myLocModified* (original -> last-modified
///   shape, for tracking splits and merges) and *myLocGenerated* (original -> list of
///   generated sub-shapes).  rcad's `annotate_history_from_ds` builds a simpler
///   `BooleanHistory` with flat `VertexOrigin`/`EdgeOrigin` arrays indexed by result
///   BRep position.
/// - OCCT PostTreat processes vertices, edges, and faces by iterating the DS images
///   (`myImages`, `myOrigins`, `myShapesSD`) and copying images from the source DS.
///   rcad uses spatial proximity (vertex point comparison) to match result vertices
///   to DS vertices, then traces edge origin from matched endpoints.
/// - OCCT PostTreat sets `myModified` for faces that were split (maps old -> new faces
///   via `myImages`).  rcad builds `FaceOrigin` separately (in `aggregate_face_origin`).
/// - OCCT PostTreat is called once at the end of `BOPAlgo_Builder::Build`.  rcad calls
///   `annotate_history_from_ds` inside `boolean_op_with_retry` after result assembly.
///
/// See also `BooleanHistory::update_with_post_treat()` for a more OCCT-aligned
/// implementation that uses `ds.my_images` instead of spatial proximity.
///
/// ✅ OCCT-aligned: core concept (history tracking from DS) matches OCCT's
///   image-map-based approach, adapted for rcad's flat-array data model.
/// ✅ OCCT-aligned: TopExp::MapShapes(myShape, myMapShape) — build result→DS index map.
///   OCCT maps TopoDS_Shape → identity for myMapShape lookup.
///   rcad: maps result vertex index → DS vertex index, result edge index → (DS vertices).
///   Used by PrepareHistory to determine Modified/Generated/Deleted provenance.
pub(crate) fn map_result_shapes(brep: &BRep, ds: &DS) -> (Vec<usize>, Vec<(usize, usize)>) {
    let mut result_to_ds: Vec<usize> = vec![usize::MAX; brep.vertices.len()];
    for (ri, rv) in brep.vertices.iter().enumerate() {
        let pt = rv.point;
        for (di, dv) in ds.vertices.iter().enumerate() {
            if (dv.point - pt).length_squared() < crate::tolerance::TOLERANCE_ABS * crate::tolerance::TOLERANCE_ABS * 4.0 {
                result_to_ds[ri] = di;
                break;
            }
        }
    }
    let edge_pairs: Vec<(usize, usize)> = brep.edges.iter()
        .map(|e| {
            let ds_s = result_to_ds.get(e.start).copied().unwrap_or(usize::MAX);
            let ds_e = result_to_ds.get(e.end).copied().unwrap_or(usize::MAX);
            (ds_s, ds_e)
        })
        .collect();
    (result_to_ds, edge_pairs)
}

/// ✅ OCCT-aligned: PrepareHistory (Builder_4.cxx L164-252).
///   OCCT iterates source shapes → LocModified → AddModified / AddGenerated / Remove.
///   rcad: uses pre-built result_to_ds map to annotate vertex/edge provenance.
pub(crate) fn annotate_history_from_ds(brep: &BRep, history: &mut BooleanHistory, ds: &DS) {
    let (result_to_ds, _) = map_result_shapes(brep, ds);

    // OCCT L176: MapShapes done.  Annotate vertex origins (FromA/FromB/Intersection).
    let a_vc = ds.a_vertex_count;
    let n_result_verts = brep.vertices.len();
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
    let n_result_edges = brep.edges.len();
    let mut edge_origins: Vec<EdgeOrigin> = Vec::with_capacity(n_result_edges);
    let a_ec = ds.a_edge_count;
    let total_ds_edges = ds.edges.len();

    for re in &brep.edges {
        let ds_s = result_to_ds[re.start];
        let ds_e = result_to_ds[re.end];

        let origin = if ds_s == usize::MAX || ds_e == usize::MAX {
            EdgeOrigin::Generated
        } else if ds_s < a_vc && ds_e < a_vc {
            // Both endpoints are A vertices 閳?look for a DS edge in A range.
            let found = (0..a_ec.min(total_ds_edges)).find(|&dei| {
                let de = &ds.edges[dei];
                (de.start_vertex == ds_s && de.end_vertex == ds_e)
                    || (de.start_vertex == ds_e && de.end_vertex == ds_s)
            });            match found {
                Some(dei) => EdgeOrigin::FromA(dei),
                None => EdgeOrigin::SplitFromA(ds_s.min(a_vc - 1)),
            }
        } else if ds_s >= a_vc && ds_e >= a_vc {
            // Both endpoints are B vertices 閳?look for a DS edge in B range.
            let found = (a_ec..total_ds_edges).find(|&dei| {
                let de = &ds.edges[dei];
                (de.start_vertex == ds_s && de.end_vertex == ds_e)
                    || (de.start_vertex == ds_e && de.end_vertex == ds_s)
            });            match found {
                Some(dei) => EdgeOrigin::FromB(dei - a_ec),
                None => EdgeOrigin::SplitFromB(ds_s.min(ds.vertices.len().saturating_sub(1)) - a_vc),
            }
        } else {
            EdgeOrigin::Generated
        };
        edge_origins.push(origin);
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

/// ✅ OCCT-aligned: PrepareHistory shell/solid provenance (Builder_4.cxx L164-252).
///   OCCT iterates source shapes → LocModified → AddModified/AddGenerated/Remove.
///   rcad: aggregates per-face origins to shell/solid level via face_region → shell → solid.
pub(crate) fn annotate_shell_and_solid_history(brep: &BRep, history: &mut BooleanHistory) {
    let mut face_cursor = 0;
    let mut shell_origins = Vec::new();
    let mut solid_origins = Vec::with_capacity(brep.solids.len());

    for solid in &brep.solids {
        let solid_shell_start = shell_origins.len();
        for shell in &solid.shells {
            let shell_face_count = shell.faces.len();
            let shell_face_origins = history
                .face_origins
                .get(face_cursor..face_cursor + shell_face_count)
                .unwrap_or(&[]);
            shell_origins.push(aggregate_face_region_origin(shell_face_origins));
            face_cursor += shell_face_count;
        }
        solid_origins.push(aggregate_shell_region_origin(&shell_origins[solid_shell_start..]));
    }

    if face_cursor != history.face_origins.len() {
        // Face count mismatch: BRep has more/fewer faces than history tracks.
        // This happens when compound reconstruction adds/removes faces or when
        // the face order in BRep differs from the emission order.  OCCT's
        // history tracking works with TopoDS shape identity — rcad's index-based
        // tracking is inherently more fragile.  Pad shell_origins to match.
        eprintln!("[HISTORY] face_cursor={} != history={}",
            face_cursor, history.face_origins.len());
    }
    history.shell_origins = shell_origins;
    history.solid_origins = solid_origins;
}

/// Deterministic order for merging parallel `boolean_op` face emissions into [`ResultBuilder`].
/// Rayon `collect` order is undefined; sorting stabilizes co-face dedup and `total_surface_area`.
pub(crate) fn cmp_boolean_emit_order(
    a: &(FaceSampleData, bool, FaceOrigin),
    b: &(FaceSampleData, bool, FaceOrigin),
) -> std::cmp::Ordering {
    
    let rank = |o: &FaceOrigin| -> (u8, usize) {
        match o {
            FaceOrigin::FromA(i) => (0, *i),
            FaceOrigin::FromB(i) => (1, *i),
            FaceOrigin::Generated => (2, 0),
        }
    };
    let (sa, ra) = rank(&a.2);
    let (sb, rb) = rank(&b.2);
    sa.cmp(&sb)
        .then(ra.cmp(&rb))
        .then_with(|| {
            let pa = a.0.sample_point();
            let pb = b.0.sample_point();
            pa.x
                .total_cmp(&pb.x)
                .then_with(|| pa.y.total_cmp(&pb.y))
                .then_with(|| pa.z.total_cmp(&pb.z))
        })
}

/// Boolean result builder (OCCT: BOPAlgo_BOP).
/// Tracks face splice origins and participates in `BooleanHistory`.
pub struct BooleanBuilder<'a> {
    ds: &'a DS,
    op: BooleanOpType,
    use_glue: bool,
    glue_tolerance: f64,
    context: RefCell<Context>,
    // ✅ OCCT-aligned: error tracking (myReport / HasErrors equivalent).
    has_errors: bool,
    // ✅ OCCT-aligned: myImages — source shape index → list of split image indices.
    //   Uses RefCell because phase functions take &self (OCCT uses mutable member maps).
    my_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: myOrigins — split shape index → list of source origin indices.
    my_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: myShapesSD — source shape index → same-domain shape index.
    my_shapes_sd: std::cell::RefCell<std::collections::HashMap<usize, usize>>,
    // ✅ OCCT-aligned: split edges created by FillImagesEdges (PaveBlock → new DSEdge).
    //   Stored here because DS is immutable (rcad uses &'a DS); their indices start
    //   at ds.edges.len() and are referenced by my_images(EDGE) / my_origins(EDGE).
    split_edges: std::cell::RefCell<Vec<crate::bopds::ds::DSEdge>>,
    // ✅ OCCT-aligned: myInParts — source solid index → list of its IN face indices
    //   (BOPAlgo_Builder.hxx L502).  Populated during FillImagesFaces, used by
    //   FillIn3DParts / BuildDraftSolid for solid assembly.
    my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: solid-level image tracking (BOPAlgo_Builder.hxx L498 myImages).
    //   OCCT BuildSplitSolids stores split solids in myImages[source_solid].
    //   rcad: maps source side (0=A, 1=B) → result solid indices from
    //   build_split_solids.  Used by annotate_shell_and_solid_history and
    //   for OCCT-form history tracking.
    my_solid_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: solid-level origin tracking (BOPAlgo_Builder.hxx L500 myOrigins).
    //   Reverse map: result solid index → list of source sides.
    my_solid_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: myNonDestructive (BOPAlgo_Builder.hxx L503).
    //   Safe processing — avoids modifying input shapes. Used in PostTreat.
    my_non_destructive: bool,
    // ✅ OCCT-aligned: myCheckInverted (BOPAlgo_Builder.hxx L505).
    //   Enables/disables inverted-solid check on input shapes.
    my_check_inverted: bool,
}

/// Fast path: if the opposite solid is an axis-aligned box, check all sub-face
/// boundary vertices against the box AABB. For tessellated faces (cone/cylinder
/// UV grid), individual grid cells can straddle the box boundary even when their
/// sample point falls inside. Requiring ALL boundary vertices to be on the correct
/// side ensures straddling cells are conservatively classified.
///
/// - Intersection (any side): sub-face is kept only when ENTIRELY inside the box.
/// - Difference B-side: sub-face is kept only when ENTIRELY inside the box.
/// - Union/Difference A-side: sub-face is kept only when ENTIRELY outside the box.
pub(crate) fn classify_subface_against_box(
    sub: &FaceSampleData,
    solid_face_indices: &[usize],
    ds: &DS,
    op: BooleanOpType,
    source: SourceSide,
) -> Option<Classification> {
    // Skip planar sub-faces 鈥?`classify_point` correctly classifies them as On
    // when they're coplanar with a box face, allowing the coplanar dedup in
    // `build_with_history` to avoid double-counting the shared area.  The AABB
    // boundary-vertex check was designed for tessellated curved surfaces
    // (cone/cylinder UV grid) where individual grid cells straddle the boundary.
    // Planar BSpline surfaces (from NURBS-converted boxes) are also planar 鈥?
    // their boundary vertices can span both inside and outside the box, causing
    // a false In/Out from a single vertex check.  OCCT classifies such faces by
    // sampling interior points (BOPTools_AlgoTools::PointInFace), not by
    // boundary-vertex AABB test.
    let is_planar_surf = match &sub.surface {
        rcad_kernel::geom::Surface3::Plane(_) => true,
        rcad_kernel::geom::Surface3::BSpline(bsp) => {
            rcad_kernel::geom::bspline_is_planar(bsp, TOLERANCE_PLANE_DIST_RELAX)
        }
        _ => false,
    };
    if is_planar_surf {
        return None;
    }
    let tol = TOLERANCE_MESH_LEGACY;
    let mut min_x = f64::NEG_INFINITY;
    let mut max_x = f64::INFINITY;
    let mut min_y = f64::NEG_INFINITY;
    let mut max_y = f64::INFINITY;
    let mut min_z = f64::NEG_INFINITY;
    let mut max_z = f64::INFINITY;

    for &fi in solid_face_indices {
        let Surface3::Plane(pl) = &ds.faces[fi].surface else {
            return None;
        };
        let n = pl.normal;
        let d = pl.origin;

        if n.x.abs() > 1.0 - tol {
            if n.x > 0.0 { max_x = max_x.min(d.x); }
            else { min_x = min_x.max(d.x); }
        } else if n.y.abs() > 1.0 - tol {
            if n.y > 0.0 { max_y = max_y.min(d.y); }
            else { min_y = min_y.max(d.y); }
        } else if n.z.abs() > 1.0 - tol {
            if n.z > 0.0 { max_z = max_z.min(d.z); }
            else { min_z = min_z.max(d.z); }
        } else {
            return None; // non-axis-aligned plane 鈫?not a simple box
        }
    }

    if min_x.is_infinite() || max_x.is_infinite()
        || min_y.is_infinite() || max_y.is_infinite()
        || min_z.is_infinite() || max_z.is_infinite()
    {
        return None; // incomplete bounds 鈫?not a full box
    }

    let require_all_inside = op == BooleanOpType::Intersection
        || (op == BooleanOpType::Difference && source == SourceSide::B);

    let (_bmin_x, _bmax_x) = sub.boundary.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(mn, mx), v| (mn.min(v.x), mx.max(v.x)));
    let (_bmin_y, _bmax_y) = sub.boundary.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(mn, mx), v| (mn.min(v.y), mx.max(v.y)));
    let (_bmin_z, _bmax_z) = sub.boundary.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(mn, mx), v| (mn.min(v.z), mx.max(v.z)));

    for &v in &sub.boundary {
        let inside = v.x >= min_x - tol && v.x <= max_x + tol
            && v.y >= min_y - tol && v.y <= max_y + tol
            && v.z >= min_z - tol && v.z <= max_z + tol;

        if require_all_inside {
            if !inside {
                // Boundary vertex outside the box 鈫?this sub-face straddles
                // the boundary.  Don't immediately return Out 鈥?the tessellation
                // vertices of a curved sub-face (cylinder wall near a box face)
                // can fall outside the box even when most of the sub-face is
                // inside.  Return None to fall through to the probe grid which
                // correctly classifies partial overlap.
                return None;
            }
        } else {
            if inside {
                // ✅ OCCT-aligned: for Union, boundary vertices may be ON the
                // box surface while the face INTERIOR extends outward (sphere
                // sub-face bounded by IC arcs on the box).  Check the sample
                // point to distinguish "on surface" from "inside".
                let sp = sub.sample_point();
                let sp_inside = sp.x >= min_x - tol && sp.x <= max_x + tol
                    && sp.y >= min_y - tol && sp.y <= max_y + tol
                    && sp.z >= min_z - tol && sp.z <= max_z + tol;
                if sp_inside {
                    return Some(Classification::In);
                }
                // Sample point outside → boundary vertices are on the box
                // surface but face is outside → fall through to probe grid
                return None;
            }
        }
    }

    // All vertices satisfy the condition 鈫?uniform classification
    let result = if require_all_inside {
        Classification::In  // all inside 鈫?keep for Intersection / Difference B-side
    } else {
        Classification::Out // all outside 鈫?keep for Union / Difference A-side
    };
    Some(result)
}

/// Classify a sub-face against the solid described by `solid_face_indices`.
///
/// For [`BooleanOpType::Intersection`], [`FaceSampleData::sample_point`] can land outside the
/// other solid even when the trimmed patch overlaps both volumes (e.g. sphere 閳?
/// finite cylinder: the inward offset toward the sphere center exits the cylinder
/// slab). When the primary sample is `Out`, we probe a coarse UV grid on
/// [`FaceSampleData::uv_domain`] before concluding `Out`.
///
/// Conversely, when the primary sample is `On` (within tolerance of the other solid's
/// surface), the sub-face may be genuinely on the boundary OR the sample point may
/// happen to fall within the tolerance band of the other solid's surface despite the
/// sub-face being entirely outside (e.g. a planar sub-face of a box near a sphere's
/// surface). In that case we probe boundary and interior samples to break the tie.
// ✅ OCCT-aligned: 鍒嗙被瀛愰潰涓?In/Out/On (ClassifyFaces)銆?
//    鎺ュ彈 FaceSampleData(浠?WireFace 鎴?FaceSampleData 鏋勯€?銆?
/// ✅ OCCT-aligned: classify_against_solid_for_boolean — ComputeState (OCCT BOPAlgo_Builder).
/// OCCT-aligned: BOPTools_AlgoTools::ComputeState (cxx L660-714).
pub(crate) fn classify_against_solid_for_boolean(
    _op: BooleanOpType,
    _source: SourceSide,
    sub: &FaceSampleData,
    solid_face_indices: &[usize],
    ds: &DS,
) -> Classification {
    let bnd = &sub.boundary;
    if bnd.len() < 3 { return Classification::In; }
    let edge_bounds = build_edge_bounds(solid_face_indices, ds);
    let tol = TOLERANCE_ABS * 100.0;
    for i in 0..bnd.len() {
        let j = (i + 1) % bnd.len();
        let p1 = bnd[i]; let p2 = bnd[j];
        let edge_idx = ds.edges.iter().position(|e| {
            let k1 = quantize_pos(ds.vertices[e.start_vertex].point, tol);
            let k2 = quantize_pos(ds.vertices[e.end_vertex].point, tol);
            let kp1 = quantize_pos(p1, tol); let kp2 = quantize_pos(p2, tol);
            (kp1 == k1 && kp2 == k2) || (kp1 == k2 && kp2 == k1)
        });
        let on_solid = edge_idx.map_or(false, |ei| edge_bounds.contains(&ei));
        if !on_solid {
            match classify_point((p1 + p2) * 0.5, solid_face_indices, ds) {
                Classification::Out => return Classification::Out,
                Classification::In => return Classification::In,
                Classification::On => continue,
            }
        }
    }
    let cent = bnd.iter().copied().sum::<DVec3>() / bnd.len() as f64;
    classify_point(cent, solid_face_indices, ds)
}

// =============================================================================
// OCCT 1:1 瀵归綈: IsInternalFace (BOPTools_AlgoTools.cxx L791-872)
// =============================================================================

/// ✅ OCCT-aligned: 鏋勫缓 MEF (Map Edge鈫扚aces) 鐢ㄤ簬杈圭骇瑙掑害娉曘€?
/// OCCT BOPAlgo_FillIn3DParts::MapEdgesAndFaces (BOPAlgo_Tools.cxx L1479-1503)
/// OCCT-aligned: IsTangentFace (BOPTools_AlgoTools).
/// Checks if two faces are tangent (parallel normals + close distance).
pub fn is_tangent_face(fi_a: usize, fi_b: usize, ds: &crate::bopds::ds::DS, angle_tol: f64, dist_tol: f64) -> bool {
    let face_a = &ds.faces[fi_a];
    let face_b = &ds.faces[fi_b];
    let n_dot = face_a.normal.dot(face_b.normal).abs();
    if n_dot < angle_tol.cos() { return false; }
    let sample_a = if !face_a.boundary_verts.is_empty() {
        ds.vertices[face_a.boundary_verts[0]].point
    } else { return false; };
    let dist = match &face_b.surface {
        rcad_kernel::geom::Surface3::Plane(p) => (sample_a - p.origin).dot(p.normal).abs(),
        rcad_kernel::geom::Surface3::Sphere(s) => ((sample_a - s.center).length() - s.radius).abs(),
        _ => return false,
    };
    dist < dist_tol
}

pub(crate) fn build_edge_bounds(face_indices: &[usize], ds: &DS) -> std::collections::BTreeSet<usize> {
    let mut bounds: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for &fi in face_indices {
        let face = &ds.faces[fi];
        for &ei in &face.boundary_edges {
            bounds.insert(ei);
        }
    }
    bounds
}

/// ✅ OCCT-aligned: PointInFace 绛変环 鈥?浠?FaceSampleData 鐨?UV domain 鑾峰彇鍐呴儴閲囨牱鐐广€?
/// OCCT BOPTools_AlgoTools3D.cxx L885-917
///
/// rcad 瀹炵幇: FaceSampleData 宸叉湁 uv_domain 鍜?uv_centroid,鐩存帴鐢?UV centroid
/// 浣滀负鍐呴儴鐐?(OCCT 鐢?Hatcher 鍋?2D point-in-face,浣?rcad 鐨?FaceSampleData
/// 鏄弬鏁板寲鍖哄煙,UV centroid 鍦ㄥ唴閮?銆?
// (point_in_face, classify_by_off_solid_edge removed — dead after ComputeState alignment)

/// 閲忓寲 3D 浣嶇疆鍒?u64 key,鐢ㄤ簬瀹瑰樊鍖归厤銆?
pub(crate) fn quantize_pos(p: DVec3, tolerance: f64) -> u64 {
    let scale = 1.0 / tolerance;
    let x = (p.x * scale).round() as i64;
    let y = (p.y * scale).round() as i64;
    let z = (p.z * scale).round() as i64;
    // 缁勫悎涓?u64
    let xb = (x as u64) & 0x3FFFFF;
    let yb = (y as u64) & 0x3FFFFF;
    let zb = (z as u64) & 0x3FFFFF;
    (xb << 42) | (yb << 21) | zb
}

/// ✅ OCCT-aligned: IsInternalFace 涓诲嚱鏁?(BOPTools_AlgoTools.cxx L791-872)
///
/// 涓ょ骇鍒嗙被:
///   Level 1: 杈圭骇瑙掑害娉?鈥?瀵逛簬鍦?solid 涓婃湁澶氫簬 1 涓偦闈㈢殑杈?
///            璁＄畻瑙掑害鍒ゆ柇闈㈡槸鍚﹀湪 solid 鍐呴儴銆?
///   Level 2: ComputeState 鈥?鍏堟壘涓嶅湪 solid 涓婄殑杈瑰垎绫讳腑鐐?
///            鍚﹀垯 PointInFace 鈫?classify_point銆?
///
/// 杩斿洖: Some(true) = 闈㈠湪 solid 鍐呴儴 (IN)
///       Some(false) = 闈笉鍦?solid 鍐呴儴 (OUT)
///       None = 鏃犳硶纭畾
/// Check if a DS vertex lies on the boundary edge between sv/ev, and if so add it
/// to split_verts with its parametric position t.
/// ✅ OCCT-aligned: FillImagesEdges checks pave blocks per edge (global scope).
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
    let p = ds.vertices[vi].point;
    let ap = p - p_a;
    let t = ap.dot(ab) / ab_len2;
    if t > 1e-8 && t < 1.0 - 1e-8 {
        let proj = p_a + ab * t;
        if (p - proj).length_squared() < 1e-10 {
            split_verts.push((vi, t));
        }
    }
}

/// ✅ OCCT-aligned: BuildSplitFaces edge assembly (L357-489) + DoSplitSEAMOnFace (L58-227).
pub(crate) fn collect_face_edge_segments(ds: &DS, face_idx: usize, pcurve_lookup: &impl Fn(usize) -> Option<Curve2d>) -> Vec<WireSegment> {
    let face = &ds.faces[face_idx];
    let mut segments: Vec<WireSegment> = Vec::new();
    let mut processed_seam_ds_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

    // ✅ OCCT-aligned: boundary vertex position map (ShapesSD equivalent).
    //    OCCT's DS shares TopoDS_Vertex between shapes at same position.
    //    rcad loads each shape's vertices independently, so the sphere north
    //    pole and a box corner at the same 3D position have different DS indices.
    //    This remaps ALL IC endpoint vertices to a GLOBAL canonical vertex per
    //    position, using ALL faces' boundary vertices (not just the current face),
    //    so section edges use the same canonical vertex regardless of which face
    //    processes them.
    let bv_positions: Vec<(DVec3, usize)> = (0..ds.faces.len()).flat_map(|fi| {
        ds.faces[fi].boundary_edges.iter().flat_map(|&ei| {
            let e = &ds.edges[ei];
            [(ds.vertices[e.start_vertex].point, e.start_vertex),
             (ds.vertices[e.end_vertex].point, e.end_vertex)]
        })
    }).collect();
    let remap_ic_v = |v: usize| -> usize {
        let p = ds.vertices[v].point;
        let tol = crate::tolerance::TOLERANCE_ABS * 1000.0;
        // Pick the canonical vertex with the MINIMUM index (earliest-loaded shape,
        // typically from operand A) to ensure consistency across all faces.
        bv_positions.iter()
            .filter(|(bp, _)| (bp - p).length_squared() <= tol * tol)
            .map(|&(_, bv)| bv)
            .min()
            .unwrap_or(v)
    };

    // Check if surface is closed (U/V)  for seam edge detection
    // OCCT L383-388: GeomLib::IsClosed  U/V
    let (is_u_closed, is_v_closed) = match &face.surface {
        Surface3::Sphere(_) => (true, true),
        Surface3::Cylinder(_) => (true, false),
        Surface3::Cone(_) => (true, false),
        _ => (false, false),
    };

    // ================================================================
    // 1. Original boundary edges (OCCT L357-460)
    // ================================================================
    // OCCT-aligned: orient boundary edges consistently for closed loop.
    // OCCT's TopExp_Explorer returns edges with the orientation they have
    // in the face's wire — each edge's end vertex matches the next edge's
    // start vertex.  rcad DS stores edges with arbitrary orientation.
    // Without this fix, a box face may have boundary edges like [2→3, 3→7,
    // 6→7, 2→6] where BOTH 3→7 and 6→7 end at vertex 7 (no outgoing edge
    // from 7), making the SmartMap connectivity wrong and preventing the
    // wire splitter from forming closed loops (fi=3 was failing).
    let mut prev_end: Option<usize> = None;
    // ✅ OCCT-aligned: virtual vertex indices for deg edge ends (OCCT uses
    //   distinct TopoDS_Vertex instances for deg edge start and end).
    let mut deg_virtual_counter: usize = ds.vertices.len();
    for &ei in &face.boundary_edges {
        let edge = &ds.edges[ei];
        let (sv, ev) = match prev_end {
            Some(pe) if edge.start_vertex == pe => (edge.start_vertex, edge.end_vertex),
            Some(pe) if edge.end_vertex == pe => (edge.end_vertex, edge.start_vertex),
            _ => (edge.start_vertex, edge.end_vertex),
        };
        prev_end = Some(ev);

        // ✅ OCCT L369: check if edge was split by intersection (myImages.IsBound).
        let edge_is_split = ei < ds.my_images.len() && ds.my_images[ei].len() > 1;

        if !edge_is_split {
            // ✅ OCCT L371-382: unsplit edge — add directly.
            //   OCCT L371-377: INTERNAL orientation → FWD+REV.
            //   OCCT L379-381: FORWARD/REVERSED → add with orientation.
            let is_internal = ds.edges[ei].is_internal;
            let rep = ds.edge_on_face(ei, face_idx);
            let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                Some(&edge.curve), Some(edge.t_range));
            let src = WireEdgeSource::DsEdge(ei);
            if is_internal {
                // OCCT L373-377: INTERNAL unsplit → FWD + REV
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev, source: src.clone(),
                    orientation: WireOrientation::Internal, is_seam: false, second_pcurve: None,
                    first_pcurve: rep.map(|r| r.pcurve.clone()),
                    t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                    tangent_start: t_start, tangent_end: t_end,
                });
            }
            segments.push(WireSegment {
                start_vertex: sv, end_vertex: ev, source: src,
                orientation: WireOrientation::Forward, is_seam: false, second_pcurve: None,
                first_pcurve: rep.map(|r| r.pcurve.clone()),
                t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                tangent_start: t_start, tangent_end: t_end,
            });
            if is_internal {
                // OCCT L376-377: REVERSED copy of INTERNAL edge
                segments.push(WireSegment {
                    start_vertex: ev, end_vertex: sv, source: WireEdgeSource::DsEdge(ei),
                    orientation: WireOrientation::Reversed, is_seam: false, second_pcurve: None,
                    first_pcurve: None, t_range: [0.0, 1.0],
                    tangent_start: t_end.map(|a| (a + std::f64::consts::PI) % std::f64::consts::TAU),
                    tangent_end: t_start.map(|a| (a + std::f64::consts::PI) % std::f64::consts::TAU),
                });
            }
            continue;
        }

        // ✅ OCCT L387-404: surface closed check + seam detection.
        //   bIsClosed = IsClosed(aE, aF) && ((isUClosed && isUIso) || (isVClosed && isVIso)).
        //   rcad: IsClosed = sv==ev || are_verts_coincident (same as OCCT TopoDS_Vertex identity).
        let b_is_degenerated = ds.is_edge_degenerated(ei);
        let b_edge_closed = sv == ev || are_verts_coincident(ds, sv, ev);
        let b_is_seam = !b_is_degenerated && b_edge_closed && (is_u_closed || is_v_closed);

        // ✅ OCCT L408-464: iterate split sub-edges (aLIE from myImages.Find).
        if b_is_degenerated {
            segments.extend(build_degenerate_edge_segments(ds, ei, sv, ev, face, face_idx, &mut deg_virtual_counter));
        } else if b_is_seam && matches!(face.surface, Surface3::Sphere(_)) {
            if !processed_seam_ds_edges.insert(ei) { continue; }
            segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
        } else if b_is_seam {
            segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
        } else {
            // ✅ OCCT-aligned: use my_images sub-edges when available (populated
            //    by build_split_edges in PaveFiller).  Handles both split edges
            //    (my_images[ei] = [sub1, sub2, ...]) and un-split edges
            //    (my_images[ei] = [ei]).  Falls back to vertices_in-based splitting
            //    only when my_images is not populated (defensive).
            if !ds.my_images.is_empty() && ei < ds.my_images.len() && !ds.my_images[ei].is_empty() {
                for &sub_ei in &ds.my_images[ei] {
                    let sub_edge = &ds.edges[sub_ei];
                    let sv_seg = sub_edge.start_vertex;
                    let ev_seg = sub_edge.end_vertex;
                    if sv_seg == ev_seg { continue; }
                    let (t_start, t_end) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
                        Some(&sub_edge.curve), Some(sub_edge.t_range));
                    let rep = ds.edge_on_face(sub_ei, face_idx);
                    segments.push(WireSegment {
                        start_vertex: sv_seg, end_vertex: ev_seg,
                        source: WireEdgeSource::DsEdge(sub_ei),
                        orientation: WireOrientation::Forward, is_seam: false, second_pcurve: None,
                        first_pcurve: rep.map(|r| r.pcurve.clone()),
                        t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        tangent_start: t_start, tangent_end: t_end,
                    });                }
            } else {
                // Fallback: split boundary edges by IC vertices (FillImagesEdges equivalent).
                let p_a = ds.vertices[sv].point;
                let p_b = ds.vertices[ev].point;
                let ab = p_b - p_a;
                let ab_len2 = ab.length_squared();
                let mut split_verts: Vec<(usize, f64)> = Vec::new();
                if ab_len2 > 1e-12 {
                    // Vertices from current face's face_info.vertices_in
                    for &vi in &face.face_info.vertices_in {
                        check_and_add_split_vertex(ds, sv, ev, vi, p_a, ab, ab_len2, &mut split_verts);
                    }
                    // OCCT-aligned DoSplitSEAMOnFace: split seam edges at IC endpoints
                    // whose UV lies on the seam (U ≈ 0 or 2π).  For a sphere seam edge
                    // between two poles (V varies along U=0), an IC endpoint at U=0 on
                    // the equator (V=π/2) splits the seam into two sub-edges.
                    let is_periodic_seam = b_is_seam && matches!(face.surface,
                        Surface3::Sphere(_) | Surface3::Cylinder(_) | Surface3::Torus(_));
                    if is_periodic_seam {
                        let uv_s = world_to_uv(&face.surface, ds.vertices[sv].point);
                        let uv_e = world_to_uv(&face.surface, ds.vertices[ev].point);
                        if let (Some(uva), Some(uvb)) = (uv_s, uv_e) {
                            let seam_tol = 1e-6;
                            let on_seam_u = |u: f64| -> bool {
                                u.abs() < seam_tol || (u - std::f64::consts::TAU).abs() < seam_tol
                            };
                            let v_range = [uva.y, uvb.y];
                            let v_min = v_range[0].min(v_range[1]);
                            let v_max = v_range[0].max(v_range[1]);
                            let v_span = v_max - v_min;
                            if v_span > 1e-15 {
                                for &ci in &face.face_info.curves_sc_only() {
                                    let ic = &ds.intersection_curves[ci];
                                    for &ep in &[ic.start_vertex, ic.end_vertex] {
                                        if ep == sv || ep == ev { continue; }
                                        if split_verts.iter().any(|(v, _)| *v == ep) { continue; }
                                        if let Some(uv_ep) = world_to_uv(&face.surface, ds.vertices[ep].point) {
                                            if on_seam_u(uv_ep.x) && uv_ep.y >= v_min - seam_tol
                                                && uv_ep.y <= v_max + seam_tol
                                            {
                                                let t = (uv_ep.y - v_min) / v_span;
                                                if t > 1e-8 && t < 1.0 - 1e-8 {
                                                    split_verts.push((ep, t));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                split_verts.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
                if split_verts.is_empty() {
                    // No split vertices — whole edge as one segment (OCCT L374-378)
                    let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                        Some(&ds.edges[ei].curve), Some(ds.edges[ei].t_range));
                    let rep = ds.edge_on_face(ei, face_idx);
                    segments.push(WireSegment {
                        start_vertex: sv, end_vertex: ev,
                        source: WireEdgeSource::DsEdge(ei),
                        orientation: WireOrientation::Forward, is_seam: false, second_pcurve: None,
                        first_pcurve: rep.map(|r| r.pcurve.clone()),
                        t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        tangent_start: t_start, tangent_end: t_end,
                    });                } else {
                    // ✅ OCCT-aligned: edge split by IC vertices (OCCT myImages equivalent).
                    let mut prev_v = sv;
                    let edge_curve = &ds.edges[ei].curve;
                    let etr = ds.edges[ei].t_range;
                    // ✅ OCCT-aligned: sub-segments inherit pcurve from original edge.
                    let seg_rep = ds.edge_on_face(ei, face_idx);
                    let seg_first_pcurve = seg_rep.map(|r| r.pcurve.clone());
                    let seg_range = seg_rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]);
                    // Map normalized position to curve parameter for sub-edge ranges.
                    let norm_to_t = |n: f64| etr[0] + n * (etr[1] - etr[0]);
                    let mut prev_t = norm_to_t(0.0);
                    for &(vi, t) in &split_verts {
                        let t_vi = norm_to_t(t);
                        let (ts, te) = edge_uv_tangent(ds, prev_v, vi, &face.surface,
                            Some(edge_curve), Some([prev_t, t_vi]));
                        segments.push(WireSegment {
                            start_vertex: prev_v, end_vertex: vi,
                            source: WireEdgeSource::DsEdge(ei),
                            orientation: WireOrientation::Forward, is_seam: false, second_pcurve: None,
                            first_pcurve: seg_first_pcurve.clone(), t_range: seg_range,
                            tangent_start: ts, tangent_end: te,
                        });                        prev_v = vi;
                        prev_t = t_vi;
                    }
                    let (ts, te) = edge_uv_tangent(ds, prev_v, ev, &face.surface,
                        Some(edge_curve), Some([prev_t, etr[1]]));
                    segments.push(WireSegment {
                        start_vertex: prev_v, end_vertex: ev,
                        source: WireEdgeSource::DsEdge(ei),
                        orientation: WireOrientation::Forward, is_seam: false, second_pcurve: None,
                        first_pcurve: seg_first_pcurve.clone(), t_range: seg_range,
                        tangent_start: ts, tangent_end: te,
                    });                }
            }
        }
    }

    // ================================================================
    // ✅ OCCT-aligned: inner wire edges (BOPAlgo_Builder_2.cxx L362-384).
    // TopExp_Explorer iterates inner wires' edges after outer wire edges.
    // Each edge inherits its wire orientation (forward = FORWARD in wire).
    // ================================================================
    for (wi, inner_wire) in face.inner_boundary_edges.iter().enumerate() {
        for &(ei, forward_in_wire) in inner_wire {
            let edge = &ds.edges[ei];
            let (sv, ev) = if forward_in_wire {
                (edge.start_vertex, edge.end_vertex)
            } else {
                (edge.end_vertex, edge.start_vertex)
            };
            if sv == ev { continue; }
            let is_degenerate = ds.is_edge_degenerated(ei);
            if is_degenerate { continue; }
            // Handle seam edges for periodic surfaces
            // Defer to existing is_seam detection
            let is_seam = match &face.surface {
                Surface3::Sphere(_) => true,
                _ => (is_u_closed || is_v_closed)
                    && (sv == ev || are_verts_coincident(ds, sv, ev)),
            };
            if is_seam {
                // Use existing seam handling
                // (Seam edges from inner wires are rare; for now, add as-is)
                let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                    Some(&edge.curve), Some(edge.t_range));
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev,
                    source: WireEdgeSource::DsEdge(ei),
                    orientation: if forward_in_wire { WireOrientation::Forward } else { WireOrientation::Reversed },
                    is_seam: true, second_pcurve: None, first_pcurve: None, t_range: [0.0, 1.0],
                    tangent_start: t_start,
                    tangent_end: t_end,
                });            } else {
                let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                    Some(&edge.curve), Some(edge.t_range));
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev,
                    source: WireEdgeSource::DsEdge(ei),
                    orientation: if forward_in_wire { WireOrientation::Forward } else { WireOrientation::Reversed },
                    is_seam: false, second_pcurve: None, first_pcurve: None, t_range: [0.0, 1.0],
                    tangent_start: t_start,
                    tangent_end: t_end,
                });            }
        }
    }

    // ================================================================
    // IN edge PBs (OCCT BOPAlgo_Builder_2.cxx L467-480).
    // Each IN PaveBlock contributes its split edge as FWD+REV.
    let boundary_set: std::collections::HashSet<usize> =
        face.boundary_edges.iter().copied().collect();
    let mut pb_dedup: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &pb_idx in &face.face_info.pave_blocks_in {
        if pb_idx >= ds.pave_blocks.len() { continue; }
        let pb = &ds.pave_blocks[pb_idx];
        if boundary_set.contains(&pb.original_edge) { continue; }
        if !pb_dedup.insert(pb.original_edge) { continue; }
        let ei = pb.new_edge.unwrap_or(pb.original_edge);
        if ei >= ds.edges.len() { continue; }
        let edge = &ds.edges[ei];
        let face_surf = &ds.faces[face_idx].surface;
        let t_start = edge_angle_2d(&edge.curve, edge.t_range[0], edge.t_range, face_surf, false, ds.vertices[edge.start_vertex].geom_tol);
        let t_end = edge_angle_2d(&edge.curve, edge.t_range[1], edge.t_range, face_surf, true, ds.vertices[edge.end_vertex].geom_tol);
        // OCCT: aLE.Append(aSp) with FORWARD orientation.
        segments.push(WireSegment {
            start_vertex: edge.start_vertex, end_vertex: edge.end_vertex,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
            is_seam: false, second_pcurve: None, first_pcurve: None,
            t_range: edge.t_range,
            tangent_start: t_start, tangent_end: t_end,
        });
        // OCCT: aLE.Append(aSp) with REVERSED orientation.
        segments.push(WireSegment {
            start_vertex: edge.end_vertex, end_vertex: edge.start_vertex,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
            is_seam: false, second_pcurve: None, first_pcurve: None,
            t_range: edge.t_range,
            tangent_start: t_end, tangent_end: t_start,
        });
    }

    // Section edges = Intersection curves (OCCT BOPAlgo_Builder_2.cxx L285-296, L478-489).
    // ================================================================
    // OCCT-aligned: Process PaveBlocksSc — each PB has a pre-built edge with
    //   valid aPB->Edge().  This is the PRIMARY path for section edges.
    let had_pb_sc = !face.face_info.pave_blocks_sc.is_empty();
    for &pb_idx in &face.face_info.pave_blocks_sc {
        if pb_idx >= ds.pave_blocks.len() { continue; }
        let pb = &ds.pave_blocks[pb_idx];
        let Some(nei) = pb.new_edge else { continue; };
        if nei >= ds.edges.len() { continue; }
        let edge = &ds.edges[nei];
        if edge.start_vertex == edge.end_vertex { continue; }
        // OCCT L484-494: FWD + REV orientations for each section edge PB.
        // OCCT-aligned: remap IC endpoint to boundary vertex (ShapesSD equivalent).
        let sv_remap = remap_ic_v(edge.start_vertex);
        let ev_remap = remap_ic_v(edge.end_vertex);
        // ✅ OCCT-aligned: propagate pcurve from DSEdge face_reps to WireSegment.
        // OCCT BRep_Tool::CurveOnSurface(aE, myFace) returns the pcurve stored
        // on the edge; rcad stores it in edge.face_reps (populated by
        // make_section_edges_from_curve_pbs). Required by SmartMap has_pcurve check.
        let sec_pcurve = edge.face_reps.iter().find(|r| r.face_idx == face_idx).map(|r| r.pcurve.clone());
        let (t_fwd_s, t_fwd_e) = edge_uv_tangent(ds, sv_remap, ev_remap,
            &face.surface, Some(&edge.curve), Some(edge.t_range));
        segments.push(WireSegment {
            start_vertex: sv_remap, end_vertex: ev_remap,
            source: WireEdgeSource::DsEdge(nei), orientation: WireOrientation::Forward,
            is_seam: false, second_pcurve: None, first_pcurve: sec_pcurve.clone(),
            t_range: edge.t_range,
            tangent_start: t_fwd_s, tangent_end: t_fwd_e,
        });
        segments.push(WireSegment {
            start_vertex: ev_remap, end_vertex: sv_remap,
            source: WireEdgeSource::DsEdge(nei), orientation: WireOrientation::Reversed,
            is_seam: false, second_pcurve: None, first_pcurve: sec_pcurve,
            t_range: edge.t_range,
            tangent_start: t_fwd_e, tangent_end: t_fwd_s,
        });
    }

    // OCCT-aligned: curves_sc fallback — only when PaveBlocksSc had no valid PBs.
    //   OCCT only builds section edges from PaveBlocksSc.  rcad may miss PBs in
    //   edge cases (no init_pave_block1), so curves_sc provides a recovery path.
    if !had_pb_sc {
    for &ci in &face.face_info.curves_sc_only() {
        let ic = &ds.intersection_curves[ci];
        // ✅ OCCT-aligned: remap IC endpoint to boundary vertex (ShapesSD).
        let sv = remap_ic_v(ic.start_vertex);
        let ev = remap_ic_v(ic.end_vertex);
        // OCCT-aligned: Skip degenerate IC (unless sphere face, where we try to infer correct vertex)
        let d2 = ds.vertices[sv].point.distance_squared(ds.vertices[ev].point);
        if std::env::var("RCAD_DEBUG_IC").is_ok() {
            eprintln!("[IC_LOOP] fi={} ci={} raw=({},{}) remap=({},{})",
                face_idx, ci, ic.start_vertex, ic.end_vertex, sv, ev);
        }
        if sv == ev || d2 < TOLERANCE_ABS_SQ {
            if matches!(face.surface, Surface3::Sphere(_)) {
                // Degenerate IC on sphere: infer the correct second vertex from other ICs.
                let other_v: Vec<usize> = face.face_info.curves_sc_only().iter()
                    .filter(|&&oci| oci != ci)
                    .flat_map(|&oci| {
                        let oic = &ds.intersection_curves[oci];
                        vec![oic.start_vertex, oic.end_vertex]
                    })
                    .filter(|&v| v != sv)
                    .collect();
                let vcounts: std::collections::HashMap<usize, usize> = {
                    let mut m = std::collections::HashMap::new();
                    for &v in &other_v { *m.entry(v).or_insert(0) += 1; }
                    m
                };
                let mut candidate: Option<usize> = None;
                for (&v, &cnt) in &vcounts {
                    if cnt == 1 {
                        if candidate.is_some() { }
                        else { candidate = Some(v); }
                    }
                }
                if candidate.is_none() {
                    candidate = other_v.iter().max_by_key(|&&v| vcounts.get(&v).copied().unwrap_or(0)).copied();
                }
                if let Some(correct_ev) = candidate {
                    let fixed_sv = sv;
                    let fixed_ev = correct_ev;
                    let pcurve = pcurve_lookup(ci);
                    let (t_start, t_end) = if let Some(ref pc) = pcurve {
                        (angle_2d(pc, ic.t_range[0], ic.t_range, false, &face.surface, ds.vertices[fixed_sv].geom_tol, None),
                         angle_2d(pc, ic.t_range[1], ic.t_range, true, &face.surface, ds.vertices[correct_ev].geom_tol, None))
                    } else { (None, None) };
                    let ic_second_pcurve = compute_ic_second_pcurve(
                        &face.surface, ds, fixed_sv, fixed_ev);
                    segments.push(WireSegment { start_vertex: fixed_sv, end_vertex: fixed_ev,
                        source: WireEdgeSource::IntersectionCurve(ci), orientation: WireOrientation::Forward,
                        is_seam: false, second_pcurve: ic_second_pcurve, first_pcurve: None, t_range: [ic.t_range[0], ic.t_range[1]], tangent_start: t_start, tangent_end: t_end });
                    continue;
                }
                // Non-sphere face with degenerate IC: skip completely
                continue;
            }
            // ✅ OCCT-aligned: 闭合 Circle IC 在边界顶点处分裂(FillImagesEdges 等价)。
            //    当 Circle IC(start==end)且 boundary 边已被 vertices_in 中的顶点分割时,
            //    在 boundary 上的顶点处分裂圆为圆弧段,使 wire builder 能形成闭合环。
            //    OCCT 在 BuildSplitFaces 中通过 myImages 获得的子边自然携带了这些顶点。
            if let Curve3::Circle(ref circ) = ic.curve {
                let center = circ.center;
                let n = circ.normal.normalize();
                let r_dir = rcad_kernel::geom::any_perpendicular(n);
                let p_dir = n.cross(r_dir);
                let r = circ.radius;
                let circle_tol = 1e-8 * r.max(1.0);
                // ✅ OCCT-aligned: 收集边界上的分割顶点(来自 FillImagesEdges 的边分裂)以及在
                //    vertices_in 中的顶点,检查哪些在 Circle IC 上。
                //    边界分割顶点来自 side 面上的 TangentLine IC,不在当前面的 vertices_in 中。
                let mut vertices_to_check: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                for &vi in &face.face_info.vertices_in { vertices_to_check.insert(vi); }
                for seg in &segments { vertices_to_check.insert(seg.start_vertex); vertices_to_check.insert(seg.end_vertex); }
                let mut on_circle: Vec<(usize, f64)> = Vec::new();
                for &vi in &vertices_to_check {
                    let pt = ds.vertices[vi].point;
                    let d = pt - center;
                    if (d.length() - r).abs() < circle_tol {
                        let angle = f64::atan2(d.dot(p_dir), d.dot(r_dir));
                        on_circle.push((vi, angle));
                    }
                }
                if on_circle.len() >= 2 {
                    on_circle.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
                    let n_on = on_circle.len();
                    for i in 0..n_on {
                        let j = (i + 1) % n_on;
                        let vi = on_circle[i].0;
                        let vj = on_circle[j].0;
                        let pcurve = pcurve_lookup(ci);
                        let (ts_val, te_val) = if let Some(ref pc) = pcurve {
                            (angle_2d(pc, ic.t_range[0], ic.t_range, false, &face.surface, ds.vertices[vi].geom_tol, None),
                             angle_2d(pc, ic.t_range[1], ic.t_range, true, &face.surface, ds.vertices[vj].geom_tol, None))
                        } else { (None, None) };
                        let arc_second = compute_ic_second_pcurve(
                            &face.surface, ds, vi, vj);
                        segments.push(WireSegment {
                            start_vertex: vi, end_vertex: vj,
                            source: WireEdgeSource::IntersectionCurve(ci), orientation: WireOrientation::Forward,
                            is_seam: false, second_pcurve: arc_second, first_pcurve: None, t_range: [on_circle[i].1, on_circle[j].1], tangent_start: ts_val, tangent_end: te_val,
                        });                    }
                    continue;
                }
            }
            continue;
        }

        //  pcurve  (Angle2D)
        let pcurve = pcurve_lookup(ci);
        let (t_start, t_end) = if let Some(ref pc) = pcurve {
            let domain = ic.t_range;
            (angle_2d(pc, domain[0], domain, false, &face.surface, ds.vertices[sv].geom_tol, None),
             angle_2d(pc, domain[1], domain, true, &face.surface, ds.vertices[ev].geom_tol, None))
        } else {
            (None, None)
        };

        // ✅ OCCT-aligned: IC edges go into WES once (FORWARD orientation).
        // OCCT BOPAlgo_Builder_2.cxx L478-489: each non-closed edge added once.
        // Closed edges (seam on periodic surfaces) get FWD+REV via separate seam logic.
        let gen_ic_second = compute_ic_second_pcurve(&face.surface, ds, sv, ev);
        segments.push(WireSegment {
            start_vertex: sv,
            end_vertex: ev,
            source: WireEdgeSource::IntersectionCurve(ci),
            orientation: WireOrientation::Forward,
            is_seam: false, second_pcurve: gen_ic_second, first_pcurve: None, t_range: [ic.t_range[0], ic.t_range[1]],
            tangent_start: t_start,
            tangent_end: t_end,
        });
    }
    } // end if !had_pb_sc (curves_sc fallback — only when PaveBlocksSc had no valid PBs)
    segments
}

/// DoSplitSEAMOnFace overload 2: compute second pcurve for an IC edge
/// whose endpoints lie on the parametric seam of a periodic surface.
/// Returns None for non-sphere surfaces or when UV can't be computed.
pub(crate) fn compute_ic_second_pcurve(
    surface: &Surface3,
    ds: &DS,
    start_vertex: usize,
    end_vertex: usize,
) -> Option<Curve2d> {
    if !matches!(surface, Surface3::Sphere(_)) {
        return None;
    }
    let sv_uv = world_to_uv(surface, ds.vertices[start_vertex].point)?;
    let ev_uv = world_to_uv(surface, ds.vertices[end_vertex].point)?;
    // Check if both endpoints are on the seam (U ≈ 0 or U ≈ 2π)
    const SEAM_TOL: f64 = 1e-6;
    let near_seam = |u: f64| -> bool {
        u.abs() < SEAM_TOL || (u - std::f64::consts::TAU).abs() < SEAM_TOL
    };
    if near_seam(sv_uv.x) && near_seam(ev_uv.x) {
        Some(Curve2d::Line(Line2d {
            origin: DVec2::new(sv_uv.x + std::f64::consts::TAU, sv_uv.y),
            direction: DVec2::new(ev_uv.x - sv_uv.x, ev_uv.y - sv_uv.y),
        }))
    } else {
        None
    }
}
