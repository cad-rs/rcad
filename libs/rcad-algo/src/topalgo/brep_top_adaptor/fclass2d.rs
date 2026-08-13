// OCCT IntTools_FClass2d (IntTools_FClass2d.hxx / .cxx)
// 2D face classifier used by IntTools_Context::IsPointInFace /
// IsPointInOnFace / IsHole.
//
// OCCT IntTools_FClass2d.cxx L77-621 (Init), L625-633 (PerformInfinitePoint),
// L637-804 (Perform), L808-943 (TestOnRestriction).
//
// Builds a per-wire polygon of the face's UV boundary (adaptive sampling via
// Geom2dInt_Geom2dCurveTool::NbSamples), records wire orientation (TabOrien),
// and classifies a 2D point with CSLib_Class2d (ray casting + ON detection).
// Falls back to BRepClass_FaceClassifier when the point is ON the boundary or
// a wire is bad.
//
// rcad data-model notes:
// - OCCT TopoDS_Face ↔ rcad DS face index + ds.face_surface(fi).
// - OCCT always has pcurves on face-boundary edges; rcad's DS builds them
//   incrementally (MakePCurves runs after VF/EF). When a boundary edge has no
//   pcurve yet, its 3D curve samples are projected onto the face surface
//   (rcad architecture difference; OCCT IntTools_FClass2d.cxx L159-160 does
//   `if (aC2D.IsNull()) return;`).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topods::{Orientation, ShapeType, TShape, TWireData, tshape_flags};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::{CONFUSION, SQUARE_CONFUSION};

use crate::topalgo::shape_source::{edge_pcurve_on_face, ShapeSource};
use crate::topalgo::brep_class::face_classifier::FClassifier;
use crate::topalgo::brep_top_adaptor::class2d::{Class2d, Class2dResult};
use crate::topalgo::gcpnts::QuasiUniformDeflection;

/// OCCT TopAbs_State — classification result (TKBRep/TopAbs layer).
pub use rcad_kernel::topods::State;

// ---------------------------------------------------------------------------
// Small helpers (OCCT support functions)
// ---------------------------------------------------------------------------

/// OCCT Geom2dInt_Geom2dCurveTool::NbSamples (Geom2dInt_Geom2dCurveTool.cxx
/// L73-91) — number of sampling points for a pcurve.
fn nb_samples(c: &Curve2d, u0: f64, u1: f64) -> i32 {
    // OCCT Adaptor2d_Curve2d::NbSamples() default (Adaptor2d_Curve2d.cxx L258-261).
    let mut nbs = 20;
    if let Curve2d::Circle(a_circ) = c {
        // Try to reach deflection = eps*R, eps = 0.01.
        if a_circ.radius > 1.0 {
            let angl = 0.283079; // 2.*std::acos(1. - eps);
            let n = ((u1 - u0).abs() / angl) as i32;
            nbs = nbs.max(n);
        }
    }
    nbs
}

/// OCCT GeomInt::AdjustPeriodic (GeomInt.cxx L21-48) — translate a parameter
/// by whole periods so it lands inside [theParMin, theParMax]. Returns
/// `(new_par, offset)`.
fn adjust_periodic(
    the_par: f64,
    the_par_min: f64,
    the_par_max: f64,
    the_period: f64,
) -> (f64, f64) {
    let the_offset;
    let mut the_new_par = the_par;
    let b_min = the_par_min - the_par > 0.0;
    let b_max = the_par - the_par_max > 0.0;
    if b_min || b_max {
        let dp = if b_min {
            the_par_max - the_par
        } else {
            the_par_min - the_par
        };
        let a_nb_per = (dp / the_period).trunc(); // modf() integer part
        the_offset = a_nb_per * the_period;
        the_new_par += the_offset;
    } else {
        the_offset = 0.0;
    }
    (the_new_par, the_offset)
}

/// OCCT Poly::PolygonProperties (Poly.hxx L165-196) — signed area and
/// perimeter of a 2D polygon. Area is negative when bypassed clockwise.
fn polygon_properties(pts: &[DVec2]) -> (f64, f64) {
    let n = pts.len();
    if n < 2 {
        return (0.0, 0.0);
    }
    let a_ref_pnt = pts[0];
    let mut a_prev_pt = pts[1] - a_ref_pnt;
    let mut the_area = 0.0;
    let mut the_perimeter = a_prev_pt.length();
    for i in 2..n {
        let a_curr_pt = pts[i] - a_ref_pnt;
        let a_delta = a_prev_pt.x * a_curr_pt.y - a_prev_pt.y * a_curr_pt.x; // Crossed
        the_area += a_delta;
        the_perimeter += (a_prev_pt - a_curr_pt).length();
        a_prev_pt = a_curr_pt;
    }
    the_perimeter += a_prev_pt.length();
    the_area *= 0.5;
    (the_area, the_perimeter)
}

/// OCCT BRepTools_WireExplorer (BRepTools_WireExplorer.cxx L121-390) —
/// reorder a wire's edges into a continuous loop following the effective
/// edge orientations (V1 -> V2). For a closed wire the explorer starts at the
/// first stored edge's V1 and walks the vertex-adjacency chain. Returns the
/// edges in traversal order (only FORWARD/REVERSED edges).
fn order_wire_edges(
    ds: &dyn ShapeSource,
    edges: &[(usize, Orientation)],
) -> Vec<(usize, Orientation)> {
    use std::collections::HashMap;
    let n = edges.len();
    let mut by_start: HashMap<(u64, u32), Vec<usize>> = HashMap::new();
    let mut edge_v2: Vec<(u64, u32)> = Vec::with_capacity(n);
    for (i, &(ei, ori)) in edges.iter().enumerate() {
        let (vf, vl) = match ds.shape_at(ei) {
            Shape { data, .. } => match &*data {
                TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
                _ => continue,
            },
        };
        // OCCT TopoDS_Iterator(E, cumLoc) composes the edge Location into the
        // vertices (TopoDS_Iterator.cxx L76-78): vertex eff location =
        // edge.Location * vertex.Location. Identity fast path in
        // compose_edge_vertex_location keeps ordinary edges untouched.
        let edge_loc = ds.shape_at(ei).location;
        let locs = ds.locations();
        let vf_loc = crate::bop::algo::compose_edge_vertex_location(edge_loc, vf.location, locs);
        let vl_loc = crate::bop::algo::compose_edge_vertex_location(edge_loc, vl.location, locs);
        // OCCT TopExp::Vertices(E, V1, V2, CumOri=true) (TopExp.cxx L214-252):
        // iterate the edge's vertex sub-shapes with the edge orientation
        // compounded (TopoDS_Iterator(E, CumOri)); V1 = the first vertex whose
        // compound orientation is FORWARD, V2 = the first whose compound
        // orientation is REVERSED. This is NOT the storage order: GWedge edges
        // store [REV(high-param), FWD(low-param)], so a FORWARD edge has
        // V1=low, V2=high and a REVERSED edge V1=high, V2=low.
        let compound = |o: Orientation| match (o, ori) {
            (Orientation::Forward, Orientation::Forward)
            | (Orientation::Reversed, Orientation::Reversed) => Orientation::Forward,
            _ => Orientation::Reversed,
        };
        let o0 = compound(vf.orientation);
        let o1 = compound(vl.orientation);
        // V1 = first vertex with compound orientation FORWARD; V2 = first with
        // compound orientation REVERSED (OCCT TopExp::Vertices L214-252).
        let v1 = if o0 == Orientation::Forward {
            (vf.ptr_id(), vf_loc)
        } else {
            (vl.ptr_id(), vl_loc)
        };
        let v2 = if o1 == Orientation::Reversed {
            (vl.ptr_id(), vl_loc)
        } else if o0 == Orientation::Reversed {
            (vf.ptr_id(), vf_loc)
        } else {
            // No REVERSED vertex (closed/degenerated edge, OCCT V2 null):
            // keep V2 = V1 so the closed-loop preference below applies.
            v1
        };
        edge_v2.push(v2);
        by_start.entry(v1).or_default().push(i);
    }
    if n == 0 {
        return Vec::new();
    }
    let mut used = vec![false; n];
    let mut result = Vec::with_capacity(n);
    // Closed wire: start from the first stored edge (V1 of that edge).
    let mut cur = 0usize;
    result.push(edges[cur]);
    used[cur] = true;
    let mut cur_v2 = edge_v2[cur];
    while result.len() < n {
        // Among the unused edges starting at the current vertex, prefer a
        // closed-loop edge (V2 == current vertex) — OCCT BRepTools_WireExplorer
        // traverses a loop (cylinder/cone lateral TopEdge/BottomEdge circle)
        // in place before leaving the vertex. Without this, the cylinder wire
        // [TopR, EndF, BottomF, StartR] walks [TopR, StartR, EndF] and skips
        // the bottom circle, producing a self-crossing UV polygon. Fall back
        // to the first unused edge otherwise.
        let next_opt = by_start.get(&cur_v2).and_then(|list| {
            list.iter()
                .copied()
                .filter(|&i| !used[i])
                .min_by_key(|&i| if edge_v2[i] == cur_v2 { 0 } else { 1 })
        });
        match next_opt {
            Some(i) => {
                result.push(edges[i]);
                used[i] = true;
                cur_v2 = edge_v2[i];
            }
            None => break,
        }
    }
    // If the walk did not visit every edge (e.g. a wire containing a closed
    // edge — cylinder/cone lateral seam — where the first stored edge is the
    // closed circle), the stored order already samples the boundary correctly.
    // Fall back to the stored order rather than dropping edges.
    if result.len() < n {
        return edges.to_vec();
    }
    result
}

/// OCCT ElCLib::Parameter / ElCLib::Value on a gp_Lin2d(origin, dir) with
/// unit direction.
fn elclib_parameter(line_origin: DVec2, line_dir: DVec2, p: DVec2) -> f64 {
    (p - line_origin).dot(line_dir)
}
fn elclib_value(t: f64, line_origin: DVec2, line_dir: DVec2) -> DVec2 {
    line_origin + line_dir * t
}

/// Surface-type periodic flags, matching OCCT GeomAbs_SurfaceType checks in
/// IntTools_FClass2d::Init (L588-619) and surf->IsUPeriodic/IsVPeriodic in
/// Perform (L655-658): U periodic for cone/cylinder/torus/sphere/surface of
/// revolution; V periodic only for torus.
fn surface_periodic(surf: &Surface3) -> (bool, bool) {
    match surf {
        Surface3::Sphere(_) | Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Revolution(_) => {
            (true, false)
        }
        Surface3::Torus(_) => (true, true),
        _ => (false, false),
    }
}

/// Per-edge UV/3D evaluator.
///
/// OCCT evaluates the pcurve `C` and its 3D curve `C3d` at sample parameters.
/// rcad adds the Projected variant for edges that lack a pcurve (rcad DS is
/// built incrementally; OCCT faces always carry pcurves).
enum EdgeEval {
    /// OCCT path: 2D from the pcurve, 3D from the edge's 3D curve.
    OnPcurve { c2d: Curve2d, c3d: Option<Curve3> },
    /// rcad path: 3D curve samples projected onto the face surface.
    Projected { c3d: Curve3, surf: Surface3, ori: Orientation },
}

impl EdgeEval {
    fn point(&self, u: f64) -> (DVec2, Option<DVec3>) {
        match self {
            EdgeEval::OnPcurve { c2d, c3d } => (
                c2d.point_at(u),
                c3d.as_ref().map(|c| c.point_at(u)),
            ),
            EdgeEval::Projected { c3d, surf, ori } => {
                let p3d = c3d.point_at(u);
                let proj = rcad_kernel::base::geom_api::project::closest_point_on_surface(surf, p3d, 64);
                let uv = glam::DVec2::new(proj.params.0, proj.params.1);
                // OCCT periodic seam (BRepPrim_OneAxis::LateralWire): the seam
                // edge's End instance sits at the periodic image u=2*PI (the
                // edge has a 2D Location), while Start sits at u=0. rcad's DS
                // has no 2D edge Locations, so the End (FORWARD) image is
                // reconstructed by shifting u by one period (the projected u of
                // the seam's 3D line is u=0; the FORWARD End instance maps to
                // u=2*PI, the REVERSED Start instance stays at u=0). A seam edge
                // is a line on a U-periodic surface of revolution, or the
                // meridian circle of a sphere (BRepPrim_Sphere::SetMeridian —
                // the arc in the XZ plane).
                let is_u_per = surface_periodic(surf).0;
                let is_seam = matches!(c3d, Curve3::Line(_))
                    || (matches!(c3d, Curve3::Circle(_)) && matches!(surf, Surface3::Sphere(_)));
                let mut u_out = uv.x;
                if is_u_per && is_seam && *ori == Orientation::Forward {
                    u_out += std::f64::consts::TAU;
                }
                (DVec2::new(u_out, uv.y), Some(p3d))
            }
        }
    }

    /// OCCT Geom2dInt_Geom2dCurveTool::NbSamples on the pcurve; rcad default
    /// (20) when no pcurve is available.
    fn nb_samples(&self, u0: f64, u1: f64) -> i32 {
        match self {
            EdgeEval::OnPcurve { c2d, .. } => nb_samples(c2d, u0, u1),
            EdgeEval::Projected { .. } => 20,
        }
    }

    /// OCCT BRep_Tool::IsClosed(edge, Face) approximation — a closed pcurve
    /// (circle/ellipse) is treated like a seam edge in the sampling loop.
    fn is_closed_pcurve(&self) -> bool {
        match self {
            EdgeEval::OnPcurve { c2d, .. } => {
                matches!(c2d, Curve2d::Circle(_) | Curve2d::Ellipse(_))
            }
            EdgeEval::Projected { .. } => false,
        }
    }
}

/// OCCT IntTools_FClass2d — 2D point-in-face classifier.
pub struct FClass2d {
    /// OCCT TabClass — per-wire CSLib_Class2d classifiers.
    tab_class: Vec<Class2d>,
    /// OCCT TabOrien — per-wire orientation (1 / 0 / -1).
    tab_orien: Vec<i32>,
    /// OCCT Toluv — UV tolerance.
    toluv: f64,
    /// OCCT Face — DS face index.
    face: usize,
    /// OCCT U1/V1/U2/V2 — periodic recadre bounds.
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    /// OCCT Umin/Umax/Vmin/Vmax — UV bounds of the sampled boundary.
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
    /// OCCT myIsHole.
    my_is_hole: bool,
    /// rcad: true when every boundary edge carries a pcurve. The precise
    /// BRepClass_FaceClassifier fallback needs pcurves (it probes the boundary
    /// pcurves); rcad's primitive faces build pcurves incrementally, so for
    /// pcurve-less faces the fallback falls back to the polygon ON/domain
    /// result.
    has_pcurves: bool,
    /// rcad: raw per-wire UV boundary polygons (final orientation-normalized
    /// sampling of the boundary pcurves). Exposed for the hatcher-equivalent
    /// "point in face" (BOPTools_AlgoTools3D::PointInFace), which intersects a
    /// vertical 2D line with the boundary. Not an OCCT field — OCCT trims the
    /// exact pcurves via Geom2dHatch_Hatcher; rcad reuses this sampling.
    uv_polygons: Vec<Vec<DVec2>>,
}

impl FClass2d {
    /// OCCT constructor: IntTools_FClass2d(aFace, TolUV) → Init(Face, Toluv).
    pub fn new(ds: &dyn ShapeSource, face_idx: usize, tol_uv: f64) -> Self {
        let mut f = FClass2d::blank(face_idx, tol_uv);
        f.init(ds, face_idx, tol_uv, None);
        f
    }

    /// OCCT IntTools_Context::FClass2d(aFace) where aFace is a temporary face
    /// built from a single analyzed loop wire (BOPAlgo_BuilderFace.cxx
    /// L437-445): the loop edges form the face's only (outer) wire. The DS face
    /// index supplies the surface and the edge pcurves (BRep_Tool::CurveOnSurface
    /// matches by surface identity, so the temp face reuses the DS face's
    /// pcurves).
    pub fn new_for_loop(
        ds: &dyn ShapeSource,
        face_idx: usize,
        tol_uv: f64,
        loop_edges: &[Shape],
    ) -> Self {
        let mut f = FClass2d::blank(face_idx, tol_uv);
        f.init(ds, face_idx, tol_uv, Some(loop_edges));
        f
    }

    fn blank(face_idx: usize, tol_uv: f64) -> Self {
        FClass2d {
            tab_class: Vec::new(),
            tab_orien: Vec::new(),
            toluv: tol_uv,
            face: face_idx,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            u_min: f64::INFINITY,
            u_max: f64::NEG_INFINITY,
            v_min: f64::INFINITY,
            v_max: f64::NEG_INFINITY,
            my_is_hole: true,
            has_pcurves: true,
            uv_polygons: Vec::new(),
        }
    }

    /// OCCT IntTools_FClass2d::IsHole — the face's outer wire is a hole.
    pub fn is_hole(&self) -> bool {
        self.my_is_hole
    }

    /// rcad: the face's per-wire UV boundary polygons (the final
    /// orientation-normalized pcurve sampling, one per classifiable wire).
    /// Used by the PointInFace hatcher equivalent to intersect a vertical 2D
    /// line with the face boundary. Not an OCCT method — see `uv_polygons`.
    pub fn uv_polygons(&self) -> &[Vec<DVec2>] {
        &self.uv_polygons
    }

    /// OCCT IntTools_FClass2d::PerformInfinitePoint.
    pub fn perform_infinite_point(&self, ds: &dyn ShapeSource) -> State {
        if self.u_max == f64::NEG_INFINITY
            || self.v_max == f64::NEG_INFINITY
            || self.u_min == f64::INFINITY
            || self.v_min == f64::INFINITY
        {
            return State::In;
        }
        let p = DVec2::new(
            self.u_min - (self.u_max - self.u_min),
            self.v_min - (self.v_max - self.v_min),
        );
        self.perform(ds, p, false)
    }

    /// OCCT IntTools_FClass2d::Init (IntTools_FClass2d.cxx L77-621). When
    /// `loop_edges` is Some, the classifier is built on a temporary face whose
    /// only wire is the analyzed loop (BOPAlgo_BuilderFace.cxx L437-445).
    pub fn init(
        &mut self,
        ds: &dyn ShapeSource,
        a_face: usize,
        tol_uv: f64,
        loop_edges: Option<&[Shape]>,
    ) {
        self.toluv = tol_uv;
        self.face = a_face;
        self.my_is_hole = true;
        self.has_pcurves = true;
        self.tab_class.clear();
        self.tab_orien.clear();
        self.uv_polygons.clear();

        let a_pr_cf = CONFUSION;
        let a_pr_cf2 = a_pr_cf * a_pr_cf;

        self.u_min = f64::INFINITY; // RealLast
        self.v_min = f64::INFINITY;
        self.u_max = f64::NEG_INFINITY;
        self.v_max = f64::NEG_INFINITY;
        let mut bad_wire = 0i32;

        let surf = ds.face_surface(a_face);

        // If the face has several wires and one of them is bad, all are
        // processed so Umin/Umax/Vmin/Vmax are correct (OCCT comment L115-118).
        // The DS face's sub_shapes are flattened to edge+vertex indices by
        // prepare_faces (BOPDS_DS.cxx L1767-1773). The ordered boundary wires
        // live in the face TShape (outer_wire + inner_wires).
        let face_wire_shapes: Vec<Shape> = match loop_edges {
            Some(edges) => {
                // OCCT BOPAlgo_BuilderFace.cxx L437-439:
                //   aBB.MakeFace(aFace, aS, aLoc, aTol); aBB.Add(aFace, aWire);
                // The loop wire is the temporary face's only (outer) wire.
                let wire = Shape::new(
                    std::sync::Arc::new(TShape::Wire(TWireData {
                        my_shapes: vec![],
                        flags: tshape_flags::DEFAULT,
                        edges: edges.to_vec(),
                    })),
                    0,
                    Orientation::Forward,
                );
                vec![wire]
            }
            None => {
                let fshape = ds.shape_at(a_face);
                let face_data = match &*fshape.data {
                    TShape::Face(fd) => fd,
                    _ => return,
                };
                std::iter::once(face_data.outer_wire.clone())
                    .chain(face_data.inner_wires.iter().cloned())
                    .collect()
            }
        };

        for (wire_idx, wire_shape) in face_wire_shapes.iter().enumerate() {
            let wire_edge_shapes = match &*wire_shape.data {
                TShape::Wire(w) => w.edges.clone(),
                _ => continue,
            };
            // Edge DS indices (via the ptr_id → index map), in wire order.
            let mut wire_edge_idxs: Vec<usize> = Vec::with_capacity(wire_edge_shapes.len());
            for eshape in &wire_edge_shapes {
                if let Some(ei) = ds.map_shape_index(eshape.ptr_id(), eshape.location) {
                    wire_edge_idxs.push(ei);
                }
            }

            let mut firstpoint = 1usize; // 1 or 2
            let mut fleche_u = 0.0;
            let mut fleche_v = 0.0;
            let mut wire_is_not_empty = false;
            let mut ancien_pnt3d_initialise = false;
            let mut ancien_pnt3d = DVec3::ZERO;
            let mut seq_pnt2d: Vec<DVec2> = Vec::new();

            let mut nb_edges = wire_edge_idxs.len() as i32;
            // OCCT BRepTools_WireExplorer: reorder the wire edges into a
            // continuous loop before sampling. BRepPrim_GWedge stores the
            // min-face wires REVERSED with a permuted edge list; the
            // WireExplorer walks the vertex-adjacency chain to recover the
            // traversal order. Sampling in stored order yields a self-crossing
            // UV polygon. The face orientation is forced FORWARD by
            // IntTools_FClass2d::Init (L105), so only the wire orientation is
            // compounded into each edge (TopoDS_Iterator semantics).
            let wire_edges_ordered: Vec<(usize, Orientation)> = {
                let mut pairs: Vec<(usize, Orientation)> = Vec::new();
                for (k, &ei) in wire_edge_idxs.iter().enumerate() {
                    if ei >= ds.nb_shapes() {
                        continue;
                    }
                    let stored_ori = wire_edge_shapes
                        .get(k)
                        .map(|s| s.orientation)
                        .unwrap_or(Orientation::Forward);
                    let ori = if wire_shape.orientation == Orientation::Reversed {
                        match stored_ori {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            other => other,
                        }
                    } else {
                        stored_ori
                    };
                    if ori != Orientation::Forward && ori != Orientation::Reversed {
                        continue;
                    }
                    pairs.push((ei, ori));
                }
                order_wire_edges(ds, &pairs)
            };

            for (_k, &(ei, ori)) in wire_edges_ordered.iter().enumerate() {
                nb_edges -= 1;
                if ei >= ds.nb_shapes() {
                    continue;
                }

                // OCCT BRep_Tool::CurveOnSurface(edge, Face, pfbid, plbid).
                // rcad: pcurve may be absent (incremental DS). rcad falls back
                // to projecting the edge 3D curve (see EdgeEval::Projected).
                let eshape = ds.shape_at(ei);
                let edge_data = match &*eshape.data {
                    TShape::Edge(ed) => ed,
                    _ => continue,
                };
                // Seam (CurveOnClosedSurface) pcurves are selected by the edge
                // orientation (FORWARD → u=2*PI, REVERSED → u=0).
                let pcurve = edge_pcurve_on_face(ds, ei, a_face, ori);
                let edge_curve = edge_data.curve.clone();

                // OCCT L166-190: degenerated edge checks.
                let mut degenerated = ds.is_edge_degenerated(ei);
                if edge_data.first.is_null() || edge_data.last.is_null() {
                    degenerated = true;
                }

                let (pfbid, plbid): (f64, f64);
                let eval: EdgeEval;
                match pcurve {
                    Some((c2d, f, l)) => {
                        pfbid = f;
                        plbid = l;
                        eval = EdgeEval::OnPcurve {
                            c2d,
                            c3d: edge_curve.clone(),
                        };
                        // OCCT L167: BRep_Tool::IsClosed(edge, Face) marks the
                        // edge degenerated. Approximated by a closed pcurve.
                        if eval.is_closed_pcurve() {
                            degenerated = true;
                        }
                    }
                    None => {
                        // rcad: no pcurve — use the 3D curve projected onto the
                        // face surface. OCCT would `return;` here (always-IN),
                        // which is incorrect for the incremental rcad DS.
                        self.has_pcurves = false;
                        let Some(c3d) = edge_curve.clone() else {
                            continue;
                        };
                        let Some(sf) = surf.clone() else {
                            continue;
                        };
                        let a_r = edge_data.range;
                        pfbid = a_r[0];
                        plbid = a_r[1];
                        eval = EdgeEval::Projected { c3d, surf: sf, ori };
                    }
                }

                // OCCT L193-196: vertex tolerances (TolVertex, dead in
                // IntTools_FClass2d) — skipped.

                // OCCT L198-220: sample-based 3D degeneracy check.
                if !degenerated {
                    if let Some(c3d) = &edge_curve {
                        let p3da = c3d.point_at(0.5 * (pfbid + plbid));
                        let du = plbid - pfbid;
                        const NBSTEPS: usize = 10;
                        let a_prec2 = 0.25 * a_pr_cf * a_pr_cf;
                        degenerated = true;
                        for i in 0..=NBSTEPS {
                            let u = pfbid + i as f64 * du / NBSTEPS as f64;
                            let p3db = c3d.point_at(u);
                            let a_r2 = p3da.distance_squared(p3db);
                            if a_r2 > a_prec2 {
                                degenerated = false;
                                break;
                            }
                        }
                    }
                }

                // OCCT L228-233: nbs = NbSamples(C); if (nbs > 2) nbs *= 4.
                let mut nbs = eval.nb_samples(pfbid, plbid);
                if nbs > 2 {
                    nbs *= 4;
                }
                let mut du = (plbid - pfbid) / (nbs - 1) as f64;
                let (mut u, u_first, u_last): (f64, f64, f64);
                if ori == Orientation::Forward {
                    u = pfbid;
                    u_first = pfbid;
                    u_last = plbid;
                } else {
                    u = plbid;
                    u_first = plbid;
                    u_last = pfbid;
                    du = -du;
                }

                // OCCT L251-270: aPrms parameter array.
                let mut a_nbs1 = nbs + 1;
                let mut a_prms: Vec<f64> = Vec::with_capacity(a_nbs1 as usize + 1);
                a_prms.push(0.0); // 1-based filler
                if nbs == 2 {
                    let a_coef = 0.0025;
                    a_prms.push(u_first);
                    a_prms.push(u_first + a_coef * (u_last - u_first));
                    a_prms.push(u_last);
                } else if nbs > 2 {
                    a_prms.push(u_first);
                    for i_x in 2..a_nbs1 {
                        a_prms.push(u + (i_x - 1) as f64 * du);
                    }
                    a_prms.push(u_last);
                } else {
                    // nbs < 2: degenerate parameter list.
                    a_nbs1 = 1;
                    a_prms.push(u_first);
                }

                // OCCT L277-365: sample loop.
                let avant = seq_pnt2d.len();
                let mut prev_u = a_prms[1];
                for i_x in firstpoint..=a_nbs1 as usize {
                    u = a_prms[i_x];
                    let (p2d, p3d_opt) = eval.point(u);
                    if p2d.x < self.u_min {
                        self.u_min = p2d.x;
                    }
                    if p2d.x > self.u_max {
                        self.u_max = p2d.x;
                    }
                    if p2d.y < self.v_min {
                        self.v_min = p2d.y;
                    }
                    if p2d.y > self.v_max {
                        self.v_max = p2d.y;
                    }

                    let mut a_dst_x = f64::MAX; // RealLast
                    let p3d = if !degenerated {
                        p3d_opt
                    } else {
                        None
                    };
                    if let Some(p3) = &p3d {
                        if !seq_pnt2d.is_empty() && ancien_pnt3d_initialise {
                            a_dst_x = p3.distance_squared(ancien_pnt3d);
                        }
                    }

                    // OCCT L318-333: IsRealCurve3d.
                    let mut is_real_curve3d = true;
                    if a_dst_x < a_pr_cf2 {
                        if i_x > 1 {
                            let mid_p3d = match &eval {
                                EdgeEval::OnPcurve { c3d, .. } => c3d
                                    .as_ref()
                                    .map(|c| c.point_at(0.5 * (u + prev_u))),
                                EdgeEval::Projected { c3d, .. } => {
                                    Some(c3d.point_at(0.5 * (u + prev_u)))
                                }
                            };
                            if let Some(mp3) = mid_p3d {
                                if let Some(p3) = &p3d {
                                    let a_dst_x1 = p3.distance_squared(mp3);
                                    if a_dst_x1 < a_pr_cf2 {
                                        is_real_curve3d = false;
                                    }
                                }
                            }
                        }
                    }

                    if is_real_curve3d {
                        if !degenerated {
                            if let Some(p3) = &p3d {
                                ancien_pnt3d = *p3;
                                ancien_pnt3d_initialise = true;
                            }
                        }
                        seq_pnt2d.push(p2d);
                    }

                    let ii = seq_pnt2d.len();
                    if ii > avant + 4 {
                        let a = seq_pnt2d[ii - 3];
                        let b = seq_pnt2d[ii - 1];
                        let mid = seq_pnt2d[ii - 2];
                        let chord = b - a;
                        let len = chord.length();
                        if len > 1e-15 {
                            let lin_dir = chord / len;
                            let ul = elclib_parameter(a, lin_dir, mid);
                            let pp = elclib_value(ul, a, lin_dir);
                            let d_u = (pp.x - mid.x).abs();
                            let d_v = (pp.y - mid.y).abs();
                            if d_u > fleche_u {
                                fleche_u = d_u;
                            }
                            if d_v > fleche_v {
                                fleche_v = d_v;
                            }
                        }
                    }
                    prev_u = u;
                } // for(iX=firstpoint; iX<=aNbs1; iX++)

                if bad_wire != 0 {
                    continue; // OCCT L367-372
                }
                if firstpoint == 1 {
                    firstpoint = 2;
                }
                wire_is_not_empty = true;
                // OCCT L379-419 (dead aD1Prev/aD1Next/anIndexMap derivative
                // code) — skipped: computed but never read.
            } // for each edge

            // OCCT L423-575: wire post-processing.
            if nb_edges != 0 {
                // Count with normal explorer and with the wire explorer differs.
                let p_class = vec![DVec2::ZERO, DVec2::ZERO];
                self.tab_class.push(Class2d::new(
                    &p_class,
                    fleche_u,
                    fleche_v,
                    self.u_min,
                    self.v_min,
                    self.u_max,
                    self.v_max,
                ));
                self.uv_polygons.push(Vec::new());
                bad_wire = 1;
                self.tab_orien.push(-1);
            } else if wire_is_not_empty {
                if seq_pnt2d.len() > 3 {
                    let (a_s, a_per) = polygon_properties(&seq_pnt2d);
                    let mut an_exp_thick = (2.0 * a_s.abs() / a_per).max(1e-7);
                    let mut a_defl = fleche_u.max(fleche_v);
                    let mut a_discr_defl = (a_defl * 0.1).min(an_exp_thick * 10.0);
                    let mut is_changed = false;
                    while a_defl > an_exp_thick && a_discr_defl > 1e-7 {
                        // Deflection of the polygon is too much for this ratio
                        // of area and perimeter — discretize the wire more
                        // tightly. OCCT L467-529.
                        firstpoint = 1;
                        is_changed = true;
                        seq_pnt2d.clear();
                        fleche_u = 0.0;
                        fleche_v = 0.0;
                        // Same WireExplorer loop order as the main sampling.
                        for &(ei2, ori2) in wire_edges_ordered.iter() {
                            if ei2 >= ds.nb_shapes() {
                                continue;
                            }
                            let eshape2 = ds.shape_at(ei2);
                            let ed2 = match &*eshape2.data {
                                TShape::Edge(ed) => ed,
                                _ => continue,
                            };
                            // OCCT BRep_Tool::Range(edge, Face, pfbid, plbid).
                            let Some((c2d, f, l)) = edge_pcurve_on_face(ds, ei2, a_face, ori2) else {
                                // rcad: no pcurve to re-discretize — stop the
                                // refinement loop (documented adaptation).
                                is_changed = false;
                                break;
                            };
                            if (l - f).abs() < 1e-9 {
                                continue;
                            }
                            let a_discr = QuasiUniformDeflection::new(&c2d, a_discr_defl, f, l);
                            if !a_discr.is_done() {
                                break;
                            }
                            let nbp = a_discr.nb_points() as i32;
                            let (mut i_step, mut i, mut i_end): (i32, i32, i32) = (1, 1, nbp + 1);
                            if ori2 == Orientation::Reversed {
                                i_step = -1;
                                i = nbp;
                                // OCCT L494-499: iEnd = 0 for reversed edges —
                                // the loop stops before Parameter(0).
                                i_end = 0;
                            }
                            if firstpoint == 2 {
                                i += i_step;
                            }
                            while i != i_end {
                                let a_p2d = c2d.point_at(a_discr.parameter(i as usize));
                                seq_pnt2d.push(a_p2d);
                                i += i_step;
                            }
                            if nbp > 2 {
                                let ii = seq_pnt2d.len();
                                if ii >= 3 {
                                    let a = seq_pnt2d[ii - 3];
                                    let b = seq_pnt2d[ii - 1];
                                    let mid = seq_pnt2d[ii - 2];
                                    let chord = b - a;
                                    let len = chord.length();
                                    if len > 1e-15 {
                                        let lin_dir = chord / len;
                                        let ul = elclib_parameter(a, lin_dir, mid);
                                        let pp = elclib_value(ul, a, lin_dir);
                                        let d_u = (pp.x - mid.x).abs();
                                        let d_v = (pp.y - mid.y).abs();
                                        if d_u > fleche_u {
                                            fleche_u = d_u;
                                        }
                                        if d_v > fleche_v {
                                            fleche_v = d_v;
                                        }
                                    }
                                }
                            }
                            firstpoint = 2;
                        }
                        if !is_changed {
                            break;
                        }
                        an_exp_thick = (2.0 * a_s.abs() / a_per).max(1e-7);
                        a_defl = fleche_u.max(fleche_v);
                        a_discr_defl = (a_discr_defl * 0.1).min(an_exp_thick * 10.0);
                    }

                    let (mut a_s, _) = if is_changed {
                        polygon_properties(&seq_pnt2d)
                    } else {
                        (a_s, a_per)
                    };

                    // OCCT derives the wire role from the polygon area sign
                    // (outer wire CCW → area>0 → TabOrien=1). rcad primitive
                    // faces are wound CCW from the outward side, which the
                    // surface (u,v) frame projects to CW (negative area);
                    // OCCT's faces are wound so outer wires give positive
                    // area. Normalize by wire role so TabOrien matches OCCT
                    // semantics: outer → CCW, inner (hole) → CW.
                    // The first wire in the face TShape list is the outer
                    // boundary (OCCT: outer_wire then inner_wires).
                    let is_outer = wire_idx == 0;
                    // OCCT IntTools_FClass2d.cxx L556-563: myIsHole follows the
                    // RAW signed area (aS>0 not-hole, aS<0 hole). rcad keeps
                    // the role normalization below for TabOrien (the DS
                    // classifier relies on the normalized value), but IsHole
                    // must use the raw sign, so record it before normalizing.
                    let a_s_raw_positive = a_s > 0.0;
                    if (is_outer && a_s < 0.0) || (!is_outer && a_s > 0.0) {
                        seq_pnt2d.reverse();
                        let (a_s2, _) = polygon_properties(&seq_pnt2d);
                        a_s = a_s2;
                    }

                    if fleche_u < self.toluv {
                        fleche_u = self.toluv;
                    }
                    if fleche_v < self.toluv {
                        fleche_v = self.toluv;
                    }

                    self.uv_polygons.push(seq_pnt2d.clone());
                    self.tab_class.push(Class2d::new(
                        &seq_pnt2d,
                        fleche_u,
                        fleche_v,
                        self.u_min,
                        self.v_min,
                        self.u_max,
                        self.v_max,
                    ));
                    if a_s.abs() < SQUARE_CONFUSION {
                        bad_wire = 1;
                        self.tab_orien.push(-1);
                    } else {
                        // OCCT: myIsHole from the raw sign (before the TabOrien
                        // role normalization).
                        if a_s_raw_positive {
                            self.my_is_hole = false;
                        } else {
                            self.my_is_hole = true;
                        }
                        if a_s > 0.0 {
                            self.tab_orien.push(1);
                        } else {
                            self.tab_orien.push(0);
                        }
                    }
                } else {
                    bad_wire = 1;
                    self.tab_orien.push(-1);
                    seq_pnt2d.clear();
                    self.uv_polygons.push(Vec::new());
                    self.tab_class.push(Class2d::new(
                        &seq_pnt2d,
                        fleche_u,
                        fleche_v,
                        self.u_min,
                        self.v_min,
                        self.u_max,
                        self.v_max,
                    ));
                }
            } // else if(wire_is_not_empty)
        } // for each wire

        // OCCT L578-620.
        let nbtabclass = self.tab_class.len();
        if nbtabclass > 0 {
            // If an error was detected on a wire: set all TabOrien to -1.
            if bad_wire != 0 {
                self.tab_orien[0] = -1;
            }
            let (is_u_per, is_v_per) = surf.as_ref().map(surface_periodic).unwrap_or((false, false));
            if is_u_per {
                let mut uuu = std::f64::consts::PI + std::f64::consts::PI - (self.u_max - self.u_min);
                if uuu < 0.0 {
                    uuu = 0.0;
                }
                self.u1 = self.u_min - uuu * 0.5;
                self.u2 = self.u1 + std::f64::consts::PI + std::f64::consts::PI;
            } else {
                self.u1 = 0.0;
                self.u2 = 0.0;
            }
            if is_v_per {
                let mut uuu = std::f64::consts::PI + std::f64::consts::PI - (self.v_max - self.v_min);
                if uuu < 0.0 {
                    uuu = 0.0;
                }
                self.v1 = self.v_min - uuu * 0.5;
                self.v2 = self.v1 + std::f64::consts::PI + std::f64::consts::PI;
            } else {
                self.v1 = 0.0;
                self.v2 = 0.0;
            }
        }
    }

    /// OCCT IntTools_FClass2d::Perform (IntTools_FClass2d.cxx L637-804).
    pub fn perform(&self, ds: &dyn ShapeSource, _puv: DVec2, recadre_on_periodic: bool) -> State {
        let nbtabclass = self.tab_class.len();
        if nbtabclass == 0 {
            return State::In;
        }

        // U1 is the First Param and U2 in this case is U1+Period.
        let mut u = _puv.x;
        let mut v = _puv.y;
        let mut a_status = State::Unknown;

        let surf = ds.face_surface(self.face);
        let (is_u_per, is_v_per) = surf.as_ref().map(surface_periodic).unwrap_or((false, false));
        let uperiod = if is_u_per {
            std::f64::consts::PI + std::f64::consts::PI
        } else {
            0.0
        };
        let vperiod = if is_v_per {
            std::f64::consts::PI + std::f64::consts::PI
        } else {
            0.0
        };

        let mut urecadre = false;
        let mut vrecadre = false;

        // OCCT L666-678: periodic adjustment of the query point into
        // [Umin,Umax] / [Vmin,Vmax]. The adjusted uu/vv are used on retry.
        let (uu, _du) = if recadre_on_periodic && is_u_per {
            adjust_periodic(u, self.u_min, self.u_max, uperiod)
        } else {
            (u, 0.0)
        };
        let (vv, _dv) = if recadre_on_periodic && is_v_per {
            adjust_periodic(v, self.v_min, self.v_max, vperiod)
        } else {
            (v, 0.0)
        };

        loop {
            let mut dedans = 1;
            let puv = DVec2::new(u, v);
            let mut b_use_classifier = self.tab_orien[0] == -1;
            if !b_use_classifier {
                for n in 0..nbtabclass {
                    let cur = self.tab_class[n].si_dans(puv);
                    let tab_orien_n = self.tab_orien[n];
                    if cur == Class2dResult::Inside {
                        if tab_orien_n == 0 {
                            dedans = -1;
                            break;
                        }
                    } else if cur == Class2dResult::Outside {
                        if tab_orien_n == 1 {
                            dedans = -1;
                            break;
                        }
                    } else {
                        dedans = 0;
                        break;
                    }
                }
                if dedans == 0 {
                    b_use_classifier = true;
                } else if dedans == 1 {
                    a_status = State::In;
                } else {
                    a_status = State::Out;
                }
            }
            // Compute state of the point using the face classifier.
            if b_use_classifier {
                if self.has_pcurves {
                    // OCCT IntTools_FClass2d.cxx L726-756: BRepClass_FClassifier
                    // on a BRepClass_FaceExplorer.
                    let a_fc_tol = self.classifier_tol(ds, u, v);
                    let mut a_classifier = FClassifier::new();
                    a_classifier.perform(ds, self.face, puv, a_fc_tol);
                    a_status = a_classifier.state();
                } else if dedans == 0 {
                    // rcad: no pcurves — CSLib_Class2d "uncertain" means ON.
                    a_status = State::On;
                } else {
                    // rcad: no pcurves — bad wire → surface domain.
                    a_status = self.classify_fallback(ds, puv);
                }
            }

            if !recadre_on_periodic || (!is_u_per && !is_v_per) {
                return a_status;
            }
            if a_status == State::In || a_status == State::On {
                return a_status;
            }

            if !urecadre {
                u = uu;
                urecadre = true;
            } else if is_u_per {
                u += uperiod;
            }
            if u > self.u_max || !is_u_per {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if is_v_per {
                    v += vperiod;
                }
                u = uu;
                if v > self.v_max || !is_v_per {
                    return a_status;
                }
            }
        }
    }

    /// OCCT IntTools_FClass2d::TestOnRestriction (IntTools_FClass2d.cxx L808-943).
    pub fn test_on_restriction(&self, ds: &dyn ShapeSource, _puv: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        let nbtabclass = self.tab_class.len();
        if nbtabclass == 0 {
            return State::In;
        }

        let mut u = _puv.x;
        let mut v = _puv.y;

        let surf = ds.face_surface(self.face);
        let (is_u_per, is_v_per) = surf.as_ref().map(surface_periodic).unwrap_or((false, false));
        let uperiod = if is_u_per {
            std::f64::consts::PI + std::f64::consts::PI
        } else {
            0.0
        };
        let vperiod = if is_v_per {
            std::f64::consts::PI + std::f64::consts::PI
        } else {
            0.0
        };
        let mut a_status = State::Unknown;
        let mut urecadre = false;
        let mut vrecadre = false;

        // OCCT L833-845: periodic adjustment of the query point.
        let (uu, _du) = if recadre_on_periodic && is_u_per {
            adjust_periodic(u, self.u_min, self.u_max, uperiod)
        } else {
            (u, 0.0)
        };
        let (vv, _dv) = if recadre_on_periodic && is_v_per {
            adjust_periodic(v, self.v_min, self.v_max, vperiod)
        } else {
            (v, 0.0)
        };

        loop {
            let mut dedans = 1;
            let puv = DVec2::new(u, v);

            if self.tab_orien[0] != -1 {
                for n in 0..nbtabclass {
                    let cur = self.tab_class[n].si_dans_on_mode(puv, tol);
                    if cur == Class2dResult::Inside {
                        if self.tab_orien[n] == 0 {
                            dedans = -1;
                            break;
                        }
                    } else if cur == Class2dResult::Outside {
                        if self.tab_orien[n] == 1 {
                            dedans = -1;
                            break;
                        }
                    } else {
                        dedans = 0;
                        break;
                    }
                }
                if dedans == 0 {
                    a_status = State::On;
                } else if dedans == 1 {
                    a_status = State::In;
                } else {
                    a_status = State::Out;
                }
            } else if self.has_pcurves {
                // OCCT L892-903: wrong wire → face classifier.
                let mut a_classifier = FClassifier::new();
                a_classifier.perform(ds, self.face, puv, tol);
                a_status = a_classifier.state();
            } else {
                // rcad: no pcurves — wrong wire → surface-domain fallback.
                a_status = self.classify_fallback(ds, puv);
            }

            if !recadre_on_periodic || (!is_u_per && !is_v_per) {
                return a_status;
            }
            if a_status == State::In || a_status == State::On {
                return a_status;
            }

            if !urecadre {
                u = uu;
                urecadre = true;
            } else if is_u_per {
                u += uperiod;
            }
            if u > self.u_max || !is_u_per {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if is_v_per {
                    v += vperiod;
                }
                u = uu;
                if v > self.v_max || !is_v_per {
                    return a_status;
                }
            }
        }
    }

    /// rcad fallback for faces without pcurves: a bad/degenerate wire is
    /// resolved against the face surface's natural domain (a full sphere /
    /// cylinder lateral collapses to a seam in UV, where the polygon
    /// classification is unavailable). The periodic u is wrapped into
    /// [u0, u1]; seam points are mapped to On.
    fn classify_fallback(&self, ds: &dyn ShapeSource, puv: DVec2) -> State {
        let surf = ds.face_surface(self.face);
        let Some([u0, u1, v0, v1]) = surf.as_ref().map(|s| s.default_domain()) else {
            return State::Unknown;
        };
        let (mut u, v) = (puv.x, puv.y);
        let periodic = matches!(
            surf,
            Some(rcad_kernel::geom::Surface3::Sphere(_))
                | Some(rcad_kernel::geom::Surface3::Cylinder(_))
                | Some(rcad_kernel::geom::Surface3::Cone(_))
                | Some(rcad_kernel::geom::Surface3::Torus(_))
        );
        // On the seam (u at the domain boundary): ON → outside for
        // IsPointInFace, inside for IsPointInOnFace. This u-domain cut is the
        // seam of a surface of revolution (cylinder / cone / torus: a generator
        // line at u=0/2*PI). The sphere is different: rcad make_sphere builds
        // the seam as the meridian circle in the XZ plane at v = PI/2 (equator)
        // in the sphere's (u, v) with axis = Y, so the sphere's u-domain cut is
        // a parametrization seam, NOT the face boundary — the u check is
        // skipped for spheres (the v = mid check below handles the real seam).
        let is_sphere = matches!(surf, Some(rcad_kernel::geom::Surface3::Sphere(_)));
        if periodic && (u1 - u0) > 1e-9 && !is_sphere {
            let period = u1 - u0;
            u = u0 + (u - u0) - period * ((u - u0) / period).floor();
            let seam_tol = 1e-6;
            if (u - u0).abs() < seam_tol || (u - u1).abs() < seam_tol {
                return State::On;
            }
        }
        // Sphere seam: the meridian circle in the XZ plane at v = PI/2.
        if is_sphere {
            let mid_v = (v0 + v1) * 0.5;
            if (v - mid_v).abs() < 1e-6 {
                return State::On;
            }
        }
        if u >= u0 - 1e-6 && u <= u1 + 1e-6 && v >= v0 - 1e-6 && v <= v1 + 1e-6 {
            State::In
        } else {
            State::Out
        }
    }

    /// OCCT IntTools_FClass2d.cxx L728-745 — the tolerance used by the
    /// BRepClass_FaceClassifier fallback, derived from the surface resolution
    /// and whether the point lies inside the face's UV bounds.
    fn classifier_tol(&self, ds: &dyn ShapeSource, u: f64, v: f64) -> f64 {
        let (a_u_res, a_v_res) = match ds.face_surface(self.face) {
            Some(surf) => (
                rcad_kernel::topo::topods::u_resolution_for_surface(&surf, self.toluv),
                rcad_kernel::topo::topods::v_resolution_for_surface(&surf, self.toluv),
            ),
            None => (self.toluv, self.toluv),
        };
        let b_u_in = u >= self.u_min && u <= self.u_max;
        let b_v_in = v >= self.v_min && v <= self.v_max;
        if b_u_in == b_v_in {
            a_u_res.min(a_v_res)
        } else if !b_u_in {
            a_u_res
        } else {
            a_v_res
        }
    }
}

