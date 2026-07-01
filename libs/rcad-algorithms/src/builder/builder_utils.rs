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
use crate::builder::wire_splitter::{world_to_uv, edge_uv_tangent, edge_angle_2d, are_verts_coincident, is_edge_isoline};
use crate::builder::edge_builders::{build_sphere_seam_segments, build_cylinder_seam_segments, is_split_to_reverse};

/// 鉁?OCCT-aligned: compare two Curve3 for identity (same TShape).
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
/// OCCT ref: BOPAlgo_Builder_3.cxx 鈥?`BOPAlgo_Builder::PostTreat`
/// (L1-250: builds `myLocModified` and `myLocGenerated` maps from DS images).
///
/// OCCT PostTreat algorithm (line-by-line mapping):
///   L20-40:  For each original shape, iterate sub-shapes (vertices, edges, faces).
///   L42-80:  Check `myImages[ei]` on each edge 鈫?if non-empty, record as Modified.
///   L82-110: For edges without images but present in result 鈫?record as Preserved.
///   L112-130: Generated edges (intersection edges) 鈫?record in myGenerated.
///   L132-170: For faces, check if wire edges were split 鈫?Modified; if not in
///             result 鈫?IsDeleted.
///   L172-200: Generated faces 鈫?myGenerated.
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
/// 鉁?OCCT-aligned: core concept (history tracking from DS) matches OCCT's
///   image-map-based approach, adapted for rcad's flat-array data model.
/// 鉁?OCCT-aligned: TopExp::MapShapes(myShape, myMapShape) 鈥?build result鈫扗S index map.
///   OCCT maps TopoDS_Shape 鈫?identity for myMapShape lookup.
///   rcad: maps result vertex index 鈫?DS vertex index, result edge index 鈫?(DS vertices).
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

/// 鉁?OCCT-aligned: PrepareHistory (Builder_4.cxx L164-252).
///   OCCT iterates source shapes 鈫?LocModified 鈫?AddModified / AddGenerated / Remove.
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
            // Both endpoints are A vertices 闁?look for a DS edge in A range.
            let found = (0..a_ec.min(total_ds_edges)).find(|&dei| {
                let de = &ds.edges[dei];
                (de.start_vertex == ds_s && de.end_vertex == ds_e)
                    || (de.start_vertex == ds_e && de.end_vertex == ds_s)
            });            match found {
                Some(dei) => EdgeOrigin::FromA(dei),
                None => EdgeOrigin::SplitFromA(ds_s.min(a_vc - 1)),
            }
        } else if ds_s >= a_vc && ds_e >= a_vc {
            // Both endpoints are B vertices 闁?look for a DS edge in B range.
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

/// 鉁?OCCT-aligned: PrepareHistory shell/solid provenance (Builder_4.cxx L164-252).
///   OCCT iterates source shapes 鈫?LocModified 鈫?AddModified/AddGenerated/Remove.
///   rcad: aggregates per-face origins to shell/solid level via face_region 鈫?shell 鈫?solid.
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
        // history tracking works with TopoDS shape identity 鈥?rcad's index-based
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
    // 鉁?OCCT-aligned: error tracking (myReport / HasErrors equivalent).
    has_errors: bool,
    // 鉁?OCCT-aligned: myImages 鈥?source shape index 鈫?list of split image indices.
    //   Uses RefCell because phase functions take &self (OCCT uses mutable member maps).
    my_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // 鉁?OCCT-aligned: myOrigins 鈥?split shape index 鈫?list of source origin indices.
    my_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // 鉁?OCCT-aligned: myShapesSD 鈥?source shape index 鈫?same-domain shape index.
    my_shapes_sd: std::cell::RefCell<std::collections::HashMap<usize, usize>>,
    // 鉁?OCCT-aligned: split edges created by FillImagesEdges (PaveBlock 鈫?new DSEdge).
    //   Stored here because DS is immutable (rcad uses &'a DS); their indices start
    //   at ds.edges.len() and are referenced by my_images(EDGE) / my_origins(EDGE).
    split_edges: std::cell::RefCell<Vec<crate::bopds::ds::DSEdge>>,
    // 鉁?OCCT-aligned: myInParts 鈥?source solid index 鈫?list of its IN face indices
    //   (BOPAlgo_Builder.hxx L502).  Populated during FillImagesFaces, used by
    //   FillIn3DParts / BuildDraftSolid for solid assembly.
    my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // 鉁?OCCT-aligned: solid-level image tracking (BOPAlgo_Builder.hxx L498 myImages).
    //   OCCT BuildSplitSolids stores split solids in myImages[source_solid].
    //   rcad: maps source side (0=A, 1=B) 鈫?result solid indices from
    //   build_split_solids.  Used by annotate_shell_and_solid_history and
    //   for OCCT-form history tracking.
    my_solid_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // 鉁?OCCT-aligned: solid-level origin tracking (BOPAlgo_Builder.hxx L500 myOrigins).
    //   Reverse map: result solid index 鈫?list of source sides.
    my_solid_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // 鉁?OCCT-aligned: myNonDestructive (BOPAlgo_Builder.hxx L503).
    //   Safe processing 鈥?avoids modifying input shapes. Used in PostTreat.
    my_non_destructive: bool,
    // 鉁?OCCT-aligned: myCheckInverted (BOPAlgo_Builder.hxx L505).
    //   Enables/disables inverted-solid check on input shapes.
    my_check_inverted: bool,
}

/// Fast path: if the opposite solid is an axis-aligned box, check all sub-face
/// boundary vertices against the box AABB. For tessellated faces (cone/cylinder
/// UV grid), individual grid cells can straddle the box boundary even when their
/// sample point falls inside. Requiring ALL boundary vertices to be on the correct
/// side ensures straddling cells are conservatively classified.
///
/// - Intersection (any side): face is kept only when ENTIRELY inside the box.
/// - Difference B-side: face is kept only when ENTIRELY inside the box.
/// - Union/Difference A-side: face is kept only when ENTIRELY outside the box.
pub(crate) fn classify_face_against_box(
    sub: &FaceSampleData,
    solid_face_indices: &[usize],
    ds: &DS,
    op: BooleanOpType,
    source: SourceSide,
) -> Option<Classification> {
    // Skip planar sub-faces 閳?`classify_point` correctly classifies them as On
    // when they're coplanar with a box face, allowing the coplanar dedup in
    // `build_with_history` to avoid double-counting the shared area.  The AABB
    // boundary-vertex check was designed for tessellated curved surfaces
    // (cone/cylinder UV grid) where individual grid cells straddle the boundary.
    // Planar BSpline surfaces (from NURBS-converted boxes) are also planar 閳?
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
            return None; // non-axis-aligned plane 閳?not a simple box
        }
    }

    if min_x.is_infinite() || max_x.is_infinite()
        || min_y.is_infinite() || max_y.is_infinite()
        || min_z.is_infinite() || max_z.is_infinite()
    {
        return None; // incomplete bounds 閳?not a full box
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
                // Boundary vertex outside the box 閳?this sub-face straddles
                // the boundary.  Don't immediately return Out 閳?the tessellation
                // vertices of a curved sub-face (cylinder wall near a box face)
                // can fall outside the box even when most of the sub-face is
                // inside.  Return None to fall through to the probe grid which
                // correctly classifies partial overlap.
                return None;
            }
        } else {
            if inside {
                // 鉁?OCCT-aligned: for Union, boundary vertices may be ON the
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
                // Sample point outside 鈫?boundary vertices are on the box
                // surface but face is outside 鈫?fall through to probe grid
                return None;
            }
        }
    }

    // All vertices satisfy the condition 閳?uniform classification
    let result = if require_all_inside {
        Classification::In  // all inside 閳?keep for Intersection / Difference B-side
    } else {
        Classification::Out // all outside 閳?keep for Union / Difference A-side
    };
    Some(result)
}

/// Classify a sub-face against the solid described by `solid_face_indices`.
///
/// For [`BooleanOpType::Intersection`], [`FaceSampleData::sample_point`] can land outside the
/// other solid even when the trimmed patch overlaps both volumes (e.g. sphere 闁?
/// finite cylinder: the inward offset toward the sphere center exits the cylinder
/// slab). When the primary sample is `Out`, we probe a coarse UV grid on
/// [`FaceSampleData::uv_domain`] before concluding `Out`.
///
/// Conversely, when the primary sample is `On` (within tolerance of the other solid's
/// surface), the sub-face may be genuinely on the boundary OR the sample point may
/// happen to fall within the tolerance band of the other solid's surface despite the
/// sub-face being entirely outside (e.g. a planar sub-face of a box near a sphere's
/// surface). In that case we probe boundary and interior samples to break the tie.
// 鉁?OCCT-aligned: 閸掑棛琚€涙劙娼版稉?In/Out/On (ClassifyFaces)閵?
//    閹恒儱褰?FaceSampleData(娴?WireFace 閹?FaceSampleData 閺嬪嫰鈧?閵?
/// 鉁?OCCT-aligned: classify_against_solid_for_boolean 鈥?ComputeState (OCCT BOPAlgo_Builder).
/// OCCT-aligned: BOPTools_AlgoTools::ComputeState (cxx L660-714).
pub(crate) fn compute_state(
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
// OCCT 1:1 鐎靛綊缍? IsInternalFace (BOPTools_AlgoTools.cxx L791-872)
// =============================================================================

/// 鉁?OCCT-aligned: 閺嬪嫬缂?MEF (Map Edge閳墯aces) 閻劋绨潏鍦獓鐟欐帒瀹冲▔鏇樷偓?
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

/// 鉁?OCCT-aligned: PointInFace 缁涘鐜?閳?娴?FaceSampleData 閻?UV domain 閼惧嘲褰囬崘鍛村劥闁插洦鐗遍悙骞库偓?
/// OCCT BOPTools_AlgoTools3D.cxx L885-917
///
/// rcad 鐎圭偟骞? FaceSampleData 瀹稿弶婀?uv_domain 閸?uv_centroid,閻╁瓨甯撮悽?UV centroid
/// 娴ｆ粈璐熼崘鍛村劥閻?(OCCT 閻?Hatcher 閸?2D point-in-face,娴?rcad 閻?FaceSampleData
/// 閺勵垰寮弫鏉垮閸栧搫鐓?UV centroid 閸︺劌鍞撮柈?閵?
// (point_in_face, classify_by_off_solid_edge removed 鈥?dead after ComputeState alignment)

/// 闁插繐瀵?3D 娴ｅ秶鐤嗛崚?u64 key,閻劋绨€圭懓妯婇崠褰掑帳閵?
pub(crate) fn quantize_pos(p: DVec3, tolerance: f64) -> u64 {
    let scale = 1.0 / tolerance;
    let x = (p.x * scale).round() as i64;
    let y = (p.y * scale).round() as i64;
    let z = (p.z * scale).round() as i64;
    // 缂佸嫬鎮庢稉?u64
    let xb = (x as u64) & 0x3FFFFF;
    let yb = (y as u64) & 0x3FFFFF;
    let zb = (z as u64) & 0x3FFFFF;
    (xb << 42) | (yb << 21) | zb
}

/// 鉁?OCCT-aligned: IsInternalFace 娑撹鍤遍弫?(BOPTools_AlgoTools.cxx L791-872)
///
/// 娑撱倗楠囬崚鍡欒:
///   Level 1: 鏉堝湱楠囩憴鎺戝濞?閳?鐎甸€涚艾閸?solid 娑撳﹥婀佹径姘艾 1 娑擃亪鍋﹂棃銏㈡畱鏉?
///            鐠侊紕鐣荤憴鎺戝閸掋倖鏌囬棃銏℃Ц閸氾箑婀?solid 閸愬懘鍎撮妴?
///   Level 2: ComputeState 閳?閸忓牊澹樻稉宥呮躬 solid 娑撳﹦娈戞潏鐟板瀻缁鑵戦悙?
///            閸氾箑鍨?PointInFace 閳?classify_point閵?
///
/// 鏉╂柨娲? Some(true) = 闂堛垹婀?solid 閸愬懘鍎?(IN)
///       Some(false) = 闂堫澀绗夐崷?solid 閸愬懘鍎?(OUT)
///       None = 閺冪姵纭剁涵顔肩暰
/// Check if a DS vertex lies on the boundary edge between sv/ev, and if so add it
/// to split_verts with its parametric position t.
/// 鉁?OCCT-aligned: FillImagesEdges checks pave blocks per edge (global scope).
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

/// 鉁?OCCT-aligned: BuildSplitFaces edge assembly (L357-489) + DoSplitSEAMOnFace (L58-227).
pub(crate) fn collect_face_edge_segments(ds: &DS, face_idx: usize, pcurve_lookup: &impl Fn(usize) -> Option<Curve2d>) -> Vec<WireSegment> {
    let face = &ds.faces[face_idx];
    let mut segments: Vec<WireSegment> = Vec::new();
    let mut processed_seam_ds_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

    // 鉁?OCCT-aligned: boundary vertex position map (ShapesSD equivalent).
    //   OCCT's DS shares vertices via ShapesSD during PaveFiller.
    //   rcad: vertex remapping is done in make_section_edges_from_curve_pbs,
    //   so IC endpoints already reference canonical vertices by this point.

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
    // in the face's wire 鈥?each edge's end vertex matches the next edge's
    // start vertex.  rcad DS stores edges with arbitrary orientation.
    // Without this fix, a box face may have boundary edges like [2鈫?, 3鈫?,
    // 6鈫?, 2鈫?] where BOTH 3鈫? and 6鈫? end at vertex 7 (no outgoing edge
    // from 7), making the SmartMap connectivity wrong and preventing the
    // wire splitter from forming closed loops (fi=3 was failing).
    let mut prev_end: Option<usize> = None;
    for &ei in &face.boundary_edges {
        let edge = &ds.edges[ei];
        let (sv, ev) = match prev_end {
            Some(pe) if edge.start_vertex == pe => (edge.start_vertex, edge.end_vertex),
            Some(pe) if edge.end_vertex == pe => (edge.end_vertex, edge.start_vertex),
            _ => (edge.start_vertex, edge.end_vertex),
        };
        prev_end = Some(ev);

        // 鉁?OCCT L369: check if edge was split by intersection (myImages.IsBound).
        let edge_is_split = ei < ds.my_images.len() && ds.my_images[ei].len() > 1;

        if !edge_is_split {
            // 鉁?OCCT L395-404: seam detection for unsplit edges on periodic surfaces.
            //   OCCT iterates all wire edges uniformly (no split/unsplit distinction);
            //   rcad processes unsplit edges here 鈥?must detect seam before adding.
            let b_is_degenerated = ds.is_edge_degenerated(ei);
            let b_is_seam = !b_is_degenerated && (is_u_closed || is_v_closed)
                && ds.edge_on_face(ei, face_idx).map_or(false, |rep| {
                    let (is_uiso, is_v_iso) = is_edge_isoline(&rep.pcurve, rep.pcurve_range);
                    (is_u_closed && is_uiso) || (is_v_closed && is_v_iso)
                });
            if b_is_seam {
                if matches!(face.surface, Surface3::Sphere(_)) {
                    if !processed_seam_ds_edges.insert(ei) { continue; }
                    segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
                } else {
                    segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
                }
                continue;
            }
            // 鉁?OCCT L371-382: unsplit edge 鈥?add directly.
            //   OCCT L371-377: INTERNAL orientation 鈫?FWD+REV.
            //   OCCT L379-381: FORWARD/REVERSED 鈫?add with orientation.
            let is_internal = ds.edges[ei].is_internal;
            let rep = ds.edge_on_face(ei, face_idx);
            let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                Some(&edge.curve), Some(edge.t_range));
            let src = WireEdgeSource::DsEdge(ei);
            if is_internal {
                // OCCT L373-377: INTERNAL unsplit 鈫?FWD + REV
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev, source: src.clone(),
                    orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                    first_pcurve: rep.map(|r| r.pcurve.clone()),
                    t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                });
                segments.push(WireSegment {
                    start_vertex: ev, end_vertex: sv, source: src,
                    orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
                    first_pcurve: None, t_range: [0.0, 1.0],
                });
            } else {
                // OCCT L379-381: non-INTERNAL 鈫?add with orientation
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev, source: src,
                    orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                    first_pcurve: rep.map(|r| r.pcurve.clone()),
                    t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                });
            }
            continue;
        }

        // 鉁?OCCT L395-404: bIsClosed via IsClosed + IsEdgeIsoline.
        // On U-closed periodic surfaces (Sphere, Cylinder, Cone), seam edges
        // are U-isolines. Vertex coincidence NOT required (sphere pole-to-pole).
        let b_is_degenerated = ds.is_edge_degenerated(ei);
        let b_is_seam = if !b_is_degenerated && (is_u_closed || is_v_closed) {
            if let Some(rep) = ds.edge_on_face(ei, face_idx) {
                let (is_uiso, is_v_iso) = is_edge_isoline(&rep.pcurve, rep.pcurve_range);
                (is_u_closed && is_uiso) || (is_v_closed && is_v_iso)
            } else {
                false
            }
        } else {
            false
        };

        // 鉁?OCCT L408-464: iterate split sub-edges (aLIE from myImages.Find).
        if b_is_degenerated {
            // OCCT L413-417: iterate sub-edges, set orientation, append
            for &sub_ei in &ds.my_images[ei] {
                let sub_edge = &ds.edges[sub_ei];
                let sv_seg = sub_edge.start_vertex;
                let ev_seg = sub_edge.end_vertex;
                if sv_seg == ev_seg { continue; }
                segments.push(WireSegment {
                    start_vertex: sv_seg, end_vertex: ev_seg,
                    source: WireEdgeSource::DsEdge(sub_ei),
                    orientation: WireOrientation::Forward,
                    is_closed_on_face: true, second_pcurve: None, first_pcurve: None,
                    t_range: [0.0, 1.0],
                });
            }
        } else if b_is_seam && matches!(face.surface, Surface3::Sphere(_)) {
            if !processed_seam_ds_edges.insert(ei) { continue; }
            segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
        } else if b_is_seam {
            segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
        } else {
            // 鉁?OCCT-aligned L408-464: three-branch split edge processing.
            //   For each sub-edge from my_images, after degenerated (handled above):
            //   1. INTERNAL original (L420-426) -> FWD+REV
            //   2. Seam bIsClosed (L429-455)   -> FWD+REV with fence (handled above)
            //   3. Normal (L457-462)            -> orientation + IsSplitToReverseWithWarn
            if !ds.my_images.is_empty() && ei < ds.my_images.len() && !ds.my_images[ei].is_empty() {
                let is_original_internal = ds.edges[ei].is_internal;
                for &sub_ei in &ds.my_images[ei] {
                    let sub_edge = &ds.edges[sub_ei];
                    let sv_seg = sub_edge.start_vertex;
                    let ev_seg = sub_edge.end_vertex;
                    if sv_seg == ev_seg { continue; }
                    // OCCT L420-426: INTERNAL original -> each sub-edge FWD+REV
                    if is_original_internal {
                        let (t_start, t_end) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
                            Some(&sub_edge.curve), Some(sub_edge.t_range));
                        let rep = ds.edge_on_face(sub_ei, face_idx)
                            .or_else(|| ds.edge_on_face(ei, face_idx));
                        segments.push(WireSegment {
                            start_vertex: sv_seg, end_vertex: ev_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: rep.map(|r| r.pcurve.clone()),
                            t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        });
                        segments.push(WireSegment {
                            start_vertex: ev_seg, end_vertex: sv_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: None, t_range: [0.0, 1.0],
                        });
                        continue;
                    }
                    // OCCT L457-462: normal split -> orientation + IsSplitToReverseWithWarn
                    let needs_reverse = is_split_to_reverse(ds, sub_ei, ei);
                    let (t_fwd, t_rev) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
                        Some(&sub_edge.curve), Some(sub_edge.t_range));
                    let rep = ds.edge_on_face(sub_ei, face_idx)
                        .or_else(|| ds.edge_on_face(ei, face_idx));
                    if needs_reverse {
                        segments.push(WireSegment {
                            start_vertex: ev_seg, end_vertex: sv_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: rep.map(|r| r.pcurve.clone()),
                            t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        });
                    } else {
                        segments.push(WireSegment {
                            start_vertex: sv_seg, end_vertex: ev_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: rep.map(|r| r.pcurve.clone()),
                            t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        });
                    }                }
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
                    // whose UV lies on the seam (U 鈮?0 or 2蟺).  For a sphere seam edge
                    // between two poles (V varies along U=0), an IC endpoint at U=0 on
                    // the equator (V=蟺/2) splits the seam into two sub-edges.
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
                    // No split vertices 鈥?whole edge as one segment (OCCT L374-378)
                    let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                        Some(&ds.edges[ei].curve), Some(ds.edges[ei].t_range));
                    let rep = ds.edge_on_face(ei, face_idx);
                    segments.push(WireSegment {
                        start_vertex: sv, end_vertex: ev,
                        source: WireEdgeSource::DsEdge(ei),
                        orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                        first_pcurve: rep.map(|r| r.pcurve.clone()),
                        t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                    });                } else {
                    // 鉁?OCCT-aligned: edge split by IC vertices (OCCT myImages equivalent).
                    let mut prev_v = sv;
                    let edge_curve = &ds.edges[ei].curve;
                    let etr = ds.edges[ei].t_range;
                    // 鉁?OCCT-aligned: sub-segments inherit pcurve from original edge.
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
                            orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: seg_first_pcurve.clone(), t_range: seg_range,
                        });                        prev_v = vi;
                        prev_t = t_vi;
                    }
                    let (ts, te) = edge_uv_tangent(ds, prev_v, ev, &face.surface,
                        Some(edge_curve), Some([prev_t, etr[1]]));
                    segments.push(WireSegment {
                        start_vertex: prev_v, end_vertex: ev,
                        source: WireEdgeSource::DsEdge(ei),
                        orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                        first_pcurve: seg_first_pcurve.clone(), t_range: seg_range,
                    });                }
            }
        }
    }

    // ================================================================
    // 鉁?OCCT-aligned: inner wire edges 鈥?same processing as outer boundary.
    // OCCT TopExp_Explorer iterates all wires' edges in one loop.
    // rcad stores them separately, so we apply identical logic here.
    // ================================================================
    for inner_wire in &face.inner_boundary_edges {
        for &(ei, forward_in_wire) in inner_wire {
            let edge = &ds.edges[ei];
            let (sv, ev) = if forward_in_wire {
                (edge.start_vertex, edge.end_vertex)
            } else {
                (edge.end_vertex, edge.start_vertex)
            };
            if sv == ev { continue; }

        let edge_is_split = ei < ds.my_images.len() && ds.my_images[ei].len() > 1;

        if !edge_is_split {
            let is_internal = ds.edges[ei].is_internal;
            let rep = ds.edge_on_face(ei, face_idx);
            let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                Some(&edge.curve), Some(edge.t_range));
            let src = WireEdgeSource::DsEdge(ei);
            if is_internal {
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev, source: src.clone(),
                    orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                    first_pcurve: rep.map(|r| r.pcurve.clone()),
                    t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                });
                segments.push(WireSegment {
                    start_vertex: ev, end_vertex: sv, source: src,
                    orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
                    first_pcurve: None, t_range: [0.0, 1.0],
                });
            } else {
                segments.push(WireSegment {
                    start_vertex: sv, end_vertex: ev, source: src,
                    orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                    first_pcurve: rep.map(|r| r.pcurve.clone()),
                    t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
                });
            }
            continue;
        }

        let b_is_degenerated = ds.is_edge_degenerated(ei);
        let b_is_seam = if !b_is_degenerated && (is_u_closed || is_v_closed) {
            if let Some(rep) = ds.edge_on_face(ei, face_idx) {
                let (is_uiso, is_v_iso) = is_edge_isoline(&rep.pcurve, rep.pcurve_range);
                (is_u_closed && is_uiso) || (is_v_closed && is_v_iso)
            } else {
                false
            }
        } else {
            false
        };

        if b_is_degenerated {
            // OCCT L413-417: iterate sub-edges, set orientation, append
            for &sub_ei in &ds.my_images[ei] {
                let sub_edge = &ds.edges[sub_ei];
                let sv_seg = sub_edge.start_vertex;
                let ev_seg = sub_edge.end_vertex;
                if sv_seg == ev_seg { continue; }
                segments.push(WireSegment {
                    start_vertex: sv_seg, end_vertex: ev_seg,
                    source: WireEdgeSource::DsEdge(sub_ei),
                    orientation: WireOrientation::Forward,
                    is_closed_on_face: true, second_pcurve: None, first_pcurve: None,
                    t_range: [0.0, 1.0],
                });
            }
        } else if b_is_seam && matches!(face.surface, Surface3::Sphere(_)) {
            if !processed_seam_ds_edges.insert(ei) { continue; }
            segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
        } else if b_is_seam {
            segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
        } else {
            if !ds.my_images.is_empty() && ei < ds.my_images.len() && !ds.my_images[ei].is_empty() {
                let is_original_internal = ds.edges[ei].is_internal;
                for &sub_ei in &ds.my_images[ei] {
                    let sub_edge = &ds.edges[sub_ei];
                    let sv_seg = sub_edge.start_vertex;
                    let ev_seg = sub_edge.end_vertex;
                    if sv_seg == ev_seg { continue; }
                    if is_original_internal {
                        let (t_start, t_end) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
                            Some(&sub_edge.curve), Some(sub_edge.t_range));
                        let rep = ds.edge_on_face(sub_ei, face_idx)
                            .or_else(|| ds.edge_on_face(ei, face_idx));
                        segments.push(WireSegment {
                            start_vertex: sv_seg, end_vertex: ev_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: rep.map(|r| r.pcurve.clone()),
                            t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        });
                        segments.push(WireSegment {
                            start_vertex: ev_seg, end_vertex: sv_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: None, t_range: [0.0, 1.0],
                        });
                        continue;
                    }
                    let needs_reverse = is_split_to_reverse(ds, sub_ei, ei);
                    let (t_fwd, t_rev) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
                        Some(&sub_edge.curve), Some(sub_edge.t_range));
                    let rep = ds.edge_on_face(sub_ei, face_idx)
                        .or_else(|| ds.edge_on_face(ei, face_idx));
                    if needs_reverse {
                        segments.push(WireSegment {
                            start_vertex: ev_seg, end_vertex: sv_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: rep.map(|r| r.pcurve.clone()),
                            t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        });
                    } else {
                        segments.push(WireSegment {
                            start_vertex: sv_seg, end_vertex: ev_seg,
                            source: WireEdgeSource::DsEdge(sub_ei),
                            orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: rep.map(|r| r.pcurve.clone()),
                            t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                        });
                    }                }
            } else {
                let p_a = ds.vertices[sv].point;
                let p_b = ds.vertices[ev].point;
                let ab = p_b - p_a;
                let ab_len2 = ab.length_squared();
                let mut split_verts: Vec<(usize, f64)> = Vec::new();
                if ab_len2 > 1e-12 {
                    for &vi in &face.face_info.vertices_in {
                        check_and_add_split_vertex(ds, sv, ev, vi, p_a, ab, ab_len2, &mut split_verts);
                    }
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
                    let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
                        Some(&ds.edges[ei].curve), Some(ds.edges[ei].t_range));
                    let rep = ds.edge_on_face(ei, face_idx);
                    segments.push(WireSegment {
                        start_vertex: sv, end_vertex: ev,
                        source: WireEdgeSource::DsEdge(ei),
                        orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                        first_pcurve: rep.map(|r| r.pcurve.clone()),
                        t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
                    });                } else {
                    let mut prev_v = sv;
                    let edge_curve = &ds.edges[ei].curve;
                    let etr = ds.edges[ei].t_range;
                    let seg_rep = ds.edge_on_face(ei, face_idx);
                    let seg_first_pcurve = seg_rep.map(|r| r.pcurve.clone());
                    let seg_range = seg_rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]);
                    let norm_to_t = |n: f64| etr[0] + n * (etr[1] - etr[0]);
                    let mut prev_t = norm_to_t(0.0);
                    for &(vi, t) in &split_verts {
                        let t_vi = norm_to_t(t);
                        let (ts, te) = edge_uv_tangent(ds, prev_v, vi, &face.surface,
                            Some(edge_curve), Some([prev_t, t_vi]));
                        segments.push(WireSegment {
                            start_vertex: prev_v, end_vertex: vi,
                            source: WireEdgeSource::DsEdge(ei),
                            orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                            first_pcurve: seg_first_pcurve.clone(), t_range: seg_range,
                        });                        prev_v = vi;
                        prev_t = t_vi;
                    }
                    let (ts, te) = edge_uv_tangent(ds, prev_v, ev, &face.surface,
                        Some(edge_curve), Some([prev_t, etr[1]]));
                    segments.push(WireSegment {
                        start_vertex: prev_v, end_vertex: ev,
                        source: WireEdgeSource::DsEdge(ei),
                        orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
                        first_pcurve: seg_first_pcurve.clone(), t_range: seg_range,
                    });                }
            }
        }
        }  // end inner_wire edge loop
    }  // end inner_wire loop

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
            is_closed_on_face: false, second_pcurve: None, first_pcurve: None,
            t_range: edge.t_range,
        });
        // OCCT: aLE.Append(aSp) with REVERSED orientation.
        segments.push(WireSegment {
            start_vertex: edge.end_vertex, end_vertex: edge.start_vertex,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
            is_closed_on_face: false, second_pcurve: None, first_pcurve: None,
            t_range: edge.t_range,
        });
    }

    // OCCT-aligned: Process PaveBlocksSc — each PB has a pre-built edge with
    //   valid aPB->Edge().  This is the PRIMARY path for section edges.
    let mut sc_dedup: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &pb_idx in &face.face_info.pave_blocks_sc {
        if pb_idx >= ds.pave_blocks.len() { continue; }
        let pb = &ds.pave_blocks[pb_idx];
        let ei = pb.new_edge.unwrap_or(pb.original_edge);
        if ei >= ds.edges.len() { continue; }
        // Section edges in boundary_set are handled by the boundary_edges loop above.
        if boundary_set.contains(&ei) { continue; }
        if ds.is_edge_degenerated(ei) { continue; }
        // Dedup: each physical edge appears at most once (one FWD + one REV).
        if !sc_dedup.insert(ei) { continue; }
        let edge = &ds.edges[ei];
        // OCCT-aligned: section edges contribute FWD+REV pair to the wire,
        // matching how OCCT adds a single TopoDS_Edge with bidirectional
        // traversal in BOPAlgo_BuilderFace::PerformLoops → WireSplitter.
        let sv = edge.start_vertex;
        let ev = edge.end_vertex;
        let is_deg = (ds.vertices[sv].point - ds.vertices[ev].point).length_squared() < TOLERANCE_ABS_SQ;
        if is_deg { continue; }
        // ✅ OCCT-aligned: propagate pcurve from DSEdge face_reps to WireSegment.
        let sec_pcurve = edge.face_reps.iter().find(|r| r.face_idx == face_idx).map(|r| r.pcurve.clone());
        segments.push(WireSegment {
            start_vertex: sv, end_vertex: ev,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
            is_closed_on_face: false, second_pcurve: None, first_pcurve: sec_pcurve.clone(),
            t_range: edge.t_range,
        });
        segments.push(WireSegment {
            start_vertex: ev, end_vertex: sv,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
            is_closed_on_face: false, second_pcurve: None, first_pcurve: sec_pcurve,
            t_range: edge.t_range,
        });
    }
    segments
}
