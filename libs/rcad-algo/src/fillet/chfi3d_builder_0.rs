//! OCCT ChFi3d_Builder_0.cxx free functions — 1:1 translation of the
//! utilities used by the builder (PerformElement / PerformExtremity /
//! ExtentAnalyse / corners).
//!
//! rcad architecture note: OCCT TopoDS shapes carry their geometry through
//! global handle graphs, while rcad stores TShapes in a BRep pool, so the
//! functions here take the owning `&BRep` as an extra first argument where
//! OCCT reads geometry straight from the shape handles.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{CurveEval as _, Curve2dEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape, TShape};
use rcad_kernel::topods;

use super::chfi3d::is_tangent_faces;
use super::chfi_ds::{ChFiDSMap, ChFiDS_State, ChFiDSStripe, ChFiDSSurfData};
use super::topopebrepds::{
    TopOpeBRepDSCurvePointInterference, TopOpeBRepDSCurve, TopOpeBRepDSHDataStructure,
    TopOpeBRepDSInterference, TopOpeBRepDSKind, TopOpeBRepDSSurfaceCurveInterference,
};
use rcad_kernel::core::precision::CONFUSION;

/// OCCT TopTools_ShapeMapHasher identity: TShape + Location + Orientation.
#[inline]
pub fn shape_key(s: &Shape) -> (u64, u32, u8) {
    (s.ptr_id(), s.location, s.orientation as u8)
}

/// OCCT TopExp::Vertices(E, V1, V2) — the edge TShape stores its ends.
pub fn topexp_vertices(e: &Shape) -> (Shape, Shape) {
    let ed = e.as_edge().expect("TopExp::Vertices: not an edge");
    (ed.first.clone(), ed.last.clone())
}

/// OCCT gp_Vec::Angle — angle in [0, PI] between two vectors.
pub fn vec_angle(a: DVec3, b: DVec3) -> f64 {
    let dot = a.dot(b).clamp(-1.0, 1.0);
    // OCCT divides by magnitudes; rcad callers pass non-unit derivatives.
    let (la, lb) = (a.length(), b.length());
    if la <= 0.0 || lb <= 0.0 {
        return 0.0;
    }
    (a.dot(b) / (la * lb)).clamp(-1.0, 1.0).acos()
}

/// OCCT gp_Vec::IsParallel(Other, AngularTolerance) — same or opposite
/// direction within the angular tolerance.
pub fn vec_is_parallel(a: DVec3, b: DVec3, angular_tolerance: f64) -> bool {
    let ang = vec_angle(a, b);
    ang <= angular_tolerance || (std::f64::consts::PI - ang) <= angular_tolerance
}

/// OCCT BRep_Tool::Parameter(V, E) — vertex parameter on the edge curve.
/// rcad stores the map on the edge TShape (vertex_params keyed by vertex
/// TShape pointer); falls back to the range end matching the vertex.
pub fn brep_tool_parameter(brep: &topods::BRep, v: &Shape, e: &Shape) -> f64 {
    let ed = e.as_edge().expect("BRep_Tool::Parameter: not an edge");
    if let Some(p) = ed.vertex_params.get(&v.ptr_id()) {
        return *p;
    }
    // Architecture fallback: rcad primitive builders do not always fill
    // vertex_params; the vertex is one of the edge ends.
    if ed.first.ptr_id() == v.ptr_id() {
        return ed.range[0];
    }
    if ed.last.ptr_id() == v.ptr_id() {
        return ed.range[1];
    }
    let _ = brep;
    0.5 * (ed.range[0] + ed.range[1])
}

/// OCCT BRepTools::OriEdgeInFace(E, F) — the orientation of the edge as it
/// appears in the face's wires (Forward if absent).
pub fn brep_tools_ori_edge_in_face(brep: &topods::BRep, e: &Shape, f: &Shape) -> Orientation {
    let fd = f.as_face().expect("OriEdgeInFace: not a face");
    if let Some(wt) = brep.tshapes.get(fd.outer_wire.index) {
        if let TShape::Wire(wd) = wt.as_ref() {
            for we in &wd.edges {
                if we.is_same(e) {
                    return we.orientation;
                }
            }
        }
    }
    Orientation::Forward
}

/// OCCT TopOpeBRepTool_TOOL::Nt(uv, F, N) — face normal at UV honoring the
/// face orientation.
pub fn topopebreptool_nt(brep: &topods::BRep, uv: glam::DVec2, f: &Shape) -> Option<DVec3> {
    let fd = f.as_face()?;
    let surf = fd.surface.clone()?;
    let (_p, du, dv) = surf.derivatives(uv.x, uv.y);
    let mut n = du.cross(dv);
    if f.orientation == Orientation::Reversed {
        n = -n;
    }
    let len = n.length();
    if len <= 0.0 {
        None
    } else {
        Some(n / len)
    }
}

/// OCCT TopExp::CommonVertex(E1, E2, V) — the shared vertex of two edges
/// (ChFi3d_cherche_vertex equivalent).
pub fn topexp_common_vertex(e1: &Shape, e2: &Shape) -> Option<Shape> {
    let ed1 = e1.as_edge()?;
    let ed2 = e2.as_edge()?;
    for v1 in [&ed1.first, &ed1.last] {
        for v2 in [&ed2.first, &ed2.last] {
            if v1.is_same(v2) {
                return Some(v1.clone());
            }
        }
    }
    None
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L473-504 — ChFi3d_conexfaces.
// =========================================================================
pub fn chfi3d_conexfaces(e: &Shape, efmap: &ChFiDSMap) -> (Shape, Shape) {
    let mut f1 = Shape::null();
    let mut f2 = Shape::null();
    let list = if efmap.contains(e) { efmap.find(e).clone() } else { Vec::new() };
    for s in list {
        if f1.is_null() {
            f1 = s;
        } else {
            f2 = s;
            if !f2.is_same(&f1) {
                break;
            } else {
                f2 = Shape::null();
            }
        }
    }
    (f1, f2)
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L506-570 — ChFi3d_EdgeState.
// =========================================================================
pub fn chfi3d_edge_state(e: &[Shape; 3], efmap: &ChFiDSMap, brep: &topods::BRep) -> ChFiDS_State {
    let (f1, f2) = chfi3d_conexfaces(&e[0], efmap);
    let (f3, f4) = chfi3d_conexfaces(&e[1], efmap);
    let (f5, f6) = chfi3d_conexfaces(&e[2], efmap);

    if f2.is_null() || f4.is_null() || f6.is_null() {
        ChFiDS_State::FreeBoundary
    } else {
        let (mut o01, mut o02) = (Orientation::Forward, Orientation::Forward);
        let (mut o11, mut o12) = (Orientation::Forward, Orientation::Forward);
        let (mut o21, mut o22) = (Orientation::Forward, Orientation::Forward);
        let i = super::chfi3d::concave_side(brep, &f1, &f2, &e[0], &mut o01, &mut o02);
        let _ = i;
        let _i2 = super::chfi3d::concave_side(brep, &f3, &f4, &e[1], &mut o11, &mut o12);
        let j = super::chfi3d::concave_side(brep, &f5, &f6, &e[2], &mut o21, &mut o22);

        if o01 == o11 && o02 == o21 && o12 == o22 {
            ChFiDS_State::AllSame
        } else if o12 == o22 || i == 10 || j == 10 {
            ChFiDS_State::OnDiff
        } else {
            ChFiDS_State::OnSame
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5557-5598 — ChFi3d_ChercheBordsLibres.
// =========================================================================
pub fn chfi3d_cherche_bords_libres(
    ve_map: &ChFiDSMap,
    v1: &Shape,
) -> (bool, Shape, Shape) {
    let mut bordlibre = false;
    let mut edgelibre1 = Shape::null();
    let mut edgelibre2 = Shape::null();
    let edges = if ve_map.contains(v1) { ve_map.find(v1).clone() } else { Vec::new() };

    for cur in &edges {
        if bordlibre {
            break;
        }
        let ed = cur.as_edge().expect("not an edge");
        if ed.degenerated {
            continue;
        }
        let mut nboccur = 0;
        for cur1 in &edges {
            if cur1.is_same(cur) {
                nboccur += 1;
            }
        }
        if nboccur == 1 {
            edgelibre1 = cur.clone();
            bordlibre = true;
        }
    }
    if bordlibre {
        bordlibre = false;
        for cur in &edges {
            if bordlibre {
                break;
            }
            let ed = cur.as_edge().expect("not an edge");
            if ed.degenerated || cur.is_same(&edgelibre1) {
                continue;
            }
            let mut nboccur = 0;
            for cur1 in &edges {
                if cur1.is_same(cur) {
                    nboccur += 1;
                }
            }
            if nboccur == 1 {
                edgelibre2 = cur.clone();
                bordlibre = true;
            }
        }
    }
    (bordlibre, edgelibre1, edgelibre2)
}

/// OCCT ChFi3d_Builder_0.cxx L5604-5620 — ChFi3d_NbNotDegeneratedEdges.
pub fn chfi3d_nb_not_degenerated_edges(vtx: &Shape, ve_map: &ChFiDSMap) -> usize {
    let mut nba = if ve_map.contains(vtx) { ve_map.find(vtx).len() } else { 0 };
    let edges = if ve_map.contains(vtx) { ve_map.find(vtx).clone() } else { Vec::new() };
    for cur in &edges {
        let ed = cur.as_edge().expect("not an edge");
        if ed.degenerated {
            nba -= 1;
        }
    }
    nba
}

/// OCCT ChFi3d_Builder_0.cxx L5639-5658 — ChFi3d_NbSharpEdges.
pub fn chfi3d_nb_sharp_edges(
    vtx: &Shape,
    ve_map: &ChFiDSMap,
    efmap: &ChFiDSMap,
    brep: &topods::BRep,
) -> usize {
    let mut nba = if ve_map.contains(vtx) { ve_map.find(vtx).len() } else { 0 };
    let edges = if ve_map.contains(vtx) { ve_map.find(vtx).clone() } else { Vec::new() };
    for cur in &edges {
        let ed = cur.as_edge().expect("not an edge");
        if ed.degenerated {
            nba -= 1;
        } else {
            let (f1, f2) = chfi3d_conexfaces(cur, efmap);
            if !f2.is_null()
                && is_tangent_faces(brep, cur, &f1, &f2, crate::geomalgo::gtests_stubs::GeomAbsShape::G2)
            {
                nba -= 1;
            }
        }
    }
    nba
}

/// OCCT ChFi3d_Builder_0.cxx L5663-5678 — ChFi3d_NumberOfEdges.
pub fn chfi3d_number_of_edges(vtx: &Shape, ve_map: &ChFiDSMap) -> usize {
    let mut nba = chfi3d_nb_not_degenerated_edges(vtx, ve_map);
    let (bordlibre, _, _) = chfi3d_cherche_bords_libres(ve_map, vtx);
    if bordlibre {
        nba = (nba - 2) / 2 + 2;
    } else {
        nba /= 2;
    }
    nba
}

/// OCCT ChFi3d_Builder_0.cxx L5691-5707 — ChFi3d_NumberOfSharpEdges.
pub fn chfi3d_number_of_sharp_edges(
    vtx: &Shape,
    ve_map: &ChFiDSMap,
    efmap: &ChFiDSMap,
    brep: &topods::BRep,
) -> usize {
    let mut nba = chfi3d_nb_sharp_edges(vtx, ve_map, efmap, brep);
    let (bordlibre, _, _) = chfi3d_cherche_bords_libres(ve_map, vtx);
    if bordlibre {
        nba = (nba - 2) / 2 + 2;
    } else {
        nba /= 2;
    }
    nba
}

/// OCCT IntTools_Tools::IntermediatePoint — the mid parameter.
pub fn intermediate_point(i1: f64, i2: f64) -> f64 {
    0.5 * (i1 + i2)
}

/// OCCT Precision::PConfusion stand-in (rcad CONFUSION).
pub const P_CONFUSION: f64 = CONFUSION;

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1690-1745 — ChFi3d_ComputePCurv (the 2-point
// pcurve: axis-aligned lines or a 2-pole BSpline guaranteeing the
// parameterization Pardeb..Parfin).
// =========================================================================
pub fn chfi3d_compute_pcurv_2pt(
    uv1: glam::DVec2,
    uv2: glam::DVec2,
    pardeb: f64,
    parfin: f64,
    reverse: bool,
) -> rcad_kernel::geom::Curve2d {
    use rcad_kernel::geom::{Curve2d, Line2d};
    use glam::DVec2 as P2;
    let tol = P_CONFUSION;
    let (p1, p2) = if !reverse { (uv1, uv2) } else { (uv2, uv1) };

    if (p1.x - p2.x).abs() <= tol && ((p2.y - p1.y) - (parfin - pardeb)).abs() <= tol {
        // vertical line, growing v
        Curve2d::Line(Line2d {
            origin: P2::new(p1.x, p1.y - pardeb),
            direction: DVec2::new(0.0, 1.0),
        })
    } else if (p1.x - p2.x).abs() <= tol && ((p1.y - p2.y) - (parfin - pardeb)).abs() <= tol {
        // vertical line, decreasing v
        Curve2d::Line(Line2d {
            origin: P2::new(p1.x, p1.y + pardeb),
            direction: DVec2::new(0.0, -1.0),
        })
    } else if (p1.y - p2.y).abs() <= tol && ((p2.x - p1.x) - (parfin - pardeb)).abs() <= tol {
        // horizontal line, growing u
        Curve2d::Line(Line2d {
            origin: P2::new(p1.x - pardeb, p1.y),
            direction: DVec2::new(1.0, 0.0),
        })
    } else if (p1.y - p2.y).abs() <= tol && ((p1.x - p2.x) - (parfin - pardeb)).abs() <= tol {
        // horizontal line, decreasing u
        Curve2d::Line(Line2d {
            origin: P2::new(p1.x + pardeb, p1.y),
            direction: DVec2::new(-1.0, 0.0),
        })
    } else {
        // OCCT: 2-pole Bezier/BSpline with the imposed parameters.
        Curve2d::Bezier(rcad_kernel::geom::BezierCurve2 {
            control_points: vec![p1, p2],
            weights: vec![1.0, 1.0],
        })
    }
}

/// OCCT ChFi3d_Builder_0.cxx ChFi3d_SameParameter — the pcurve
/// same-parameter verification/correction against the surface.  Pending:
/// the identity keeps the 2-point parameterization (exact for the axis-
/// aligned pcurves produced above on analytic surfaces).
pub fn chfi3d_same_parameter(
    _c3d: &rcad_kernel::geom::Curve3,
    _pcurv: &mut rcad_kernel::geom::Curve2d,
    _s: &rcad_kernel::geom::Surface3,
    _tol3d: f64,
    tolreached: &mut f64,
) {
    *tolreached = _tol3d;
}

// =========================================================================
// OCCT ElSLib — cylinder iso curve construction used by ComputeArete.
// =========================================================================

/// OCCT Geom_CylindricalSurface::VIso(V) — the circle at height v.
pub fn cylinder_v_iso(
    origin: DVec3,
    xdir: DVec3,
    axis: DVec3,
    radius: f64,
    v: f64,
) -> rcad_kernel::geom::Circle3 {
    let ydir = axis.cross(xdir).normalize();
    let center = origin + axis * v;
    let mut c = rcad_kernel::geom::Circle3::new(center, axis, radius);
    let _ = (xdir, ydir);
    c
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1984-2180 — ChFi3d_ComputeArete.
// IFlag=0 pcurve et courbe 3d; IFlag>0 pcurve (parametrage impose si 2).
// Returns (C3d, Pcurv, Pardeb, Parfin, tolreached).
// =========================================================================
#[allow(clippy::too_many_arguments)]
pub fn chfi3d_compute_arete(
    brep: &topods::BRep,
    p1: &super::chfi_ds::ChFiDS_CommonPoint,
    uv1: glam::DVec2,
    p2: &super::chfi_ds::ChFiDS_CommonPoint,
    uv2: glam::DVec2,
    surf: &rcad_kernel::geom::Surface3,
    tol3d: f64,
    tol2d: f64,
    iflag: i32,
) -> (Option<rcad_kernel::geom::Curve3>, rcad_kernel::geom::Curve2d, f64, f64, f64) {
    use rcad_kernel::geom::{Curve2d, Curve3, Line2d};
    let mut c3d: Option<Curve3> = None;
    let pcurv;
    let mut pardeb = 0.0;
    let mut parfin = 0.0;
    let tolreached;

    if (uv1.x - uv2.x).abs() <= tol2d {
        // iso u
        if iflag == 0 {
            pardeb = uv1.y;
            parfin = uv2.y;
            // OCCT: C3d = Surf->UIso(UV1.X()) — the u-isocurve of the
            // surface at u = UV1.X(); rcad resolves the iso curve per
            // surface kind.
            let reversed = pardeb > parfin;
            if reversed {
                std::mem::swap(&mut pardeb, &mut parfin);
            }
            match surf {
                rcad_kernel::geom::Surface3::Cylinder(c) => {
                    let iso = cylinder_v_iso(c.origin, c.ref_dir, c.axis, c.radius, uv1.x);
                    c3d = Some(Curve3::Circle(iso));
                }
                _ => {
                    // pending: u-iso of non-cylindrical surfaces.
                }
            }
            if reversed {
                // OCCT reverses the curve; rcad records the range as
                // (Pardeb, Parfin) with Pardeb > Parfin dropped by the
                // swap above, so nothing more is needed here.
            }
        }
        if iflag != 1 {
            // OCCT: ChFi3d_ComputePCurv(hc, UV1, UV2, Pcurv, hs, Pardeb,
            // Parfin, tol3d, tolreached, false);
            let mut pc = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, false);
            let mut tr = tol3d;
            if let Some(c) = &c3d {
                chfi3d_same_parameter(c, &mut pc, surf, tol3d, &mut tr);
            }
            pcurv = pc;
            tolreached = tr;
        } else {
            pcurv = Curve2d::Line(Line2d {
                origin: uv1,
                direction: (uv2 - uv1).normalize(),
            });
            tolreached = tol3d;
        }
    } else if (uv1.y - uv2.y).abs() <= tol2d {
        // iso v
        if iflag == 0 {
            pardeb = uv1.x;
            parfin = uv2.x;
            let reversed = pardeb > parfin;
            if reversed {
                std::mem::swap(&mut pardeb, &mut parfin);
            }
            match surf {
                rcad_kernel::geom::Surface3::Cylinder(c) => {
                    let iso = cylinder_v_iso(c.origin, c.ref_dir, c.axis, c.radius, uv1.y);
                    c3d = Some(Curve3::Circle(iso));
                }
                _ => {
                    // pending: v-iso of non-cylindrical surfaces.
                }
            }
        }
        if iflag != 1 {
            let mut pc = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, false);
            let mut tr = tol3d;
            if let Some(c) = &c3d {
                chfi3d_same_parameter(c, &mut pc, surf, tol3d, &mut tr);
            }
            pcurv = pc;
            tolreached = tr;
        } else {
            pcurv = Curve2d::Line(Line2d {
                origin: uv1,
                direction: (uv2 - uv1).normalize(),
            });
            tolreached = tol3d;
        }
    } else if iflag == 0 {
        // OCCT L2036-2058: straight-line pcurve when a vertex is involved
        // or the points are not on arcs; otherwise the tangent-matched
        // BuildPCurve with the in-surface a-posteriori check — pending.
        pcurv = Curve2d::Bezier(rcad_kernel::geom::BezierCurve2 {
            control_points: vec![uv1, uv2],
            weights: vec![1.0, 1.0],
        });
        tolreached = tol3d;
        let _ = (p1, p2, brep);
    } else {
        // OCCT: hs->Load(Surf); hc->Load(C3d, Pardeb, Parfin);
        // ChFi3d_ProjectPCurv(...) — pending.
        pcurv = chfi3d_compute_pcurv_2pt(uv1, uv2, pardeb, parfin, false);
        tolreached = tol3d;
        let _ = (brep, p1, p2);
    }
    (c3d, pcurv, pardeb, parfin, tolreached)
}

// =========================================================================
// OCCT TopExp_Explorer stand-ins over a face (TopExp_Explorer(F, EDGE /
// VERTEX)).  rcad wire entries may be ptr-only Shapes (index =
// usize::MAX, see chfi3d architecture note 8); matching callers use
// IsSame, geometry callers resolve the BRep index by TShape pointer.
// =========================================================================

/// Resolve a possibly ptr-only wire-entry Shape to a BRep-indexed Shape,
/// keeping the entry's location and orientation.
pub fn resolve_brep_shape(brep: &topods::BRep, entry: &Shape) -> Shape {
    if entry.index != usize::MAX && brep.tshapes.get(entry.index).is_some() {
        return entry.clone();
    }
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if std::sync::Arc::as_ptr(ts) as u64 == entry.ptr_id() {
            return Shape::from_parts(ts.clone(), i, entry.location, entry.orientation);
        }
    }
    entry.clone()
}

/// OCCT TopExp_Explorer(S, TopAbs_EDGE) over the face's wires (outer +
/// inner), each edge carrying its in-face orientation.
pub fn topexp_face_edges(brep: &topods::BRep, f: &Shape) -> Vec<Shape> {
    let mut out = Vec::new();
    let Some(fd) = f.as_face() else {
        return out;
    };
    let mut wires: Vec<&Shape> = Vec::new();
    wires.push(&fd.outer_wire);
    for w in &fd.inner_wires {
        wires.push(w);
    }
    for w in wires {
        let wdata = w
            .as_wire();
        let Some(wdata) = wdata else { continue };
        for we in &wdata.edges {
            out.push(resolve_brep_shape(brep, we));
        }
    }
    out
}

/// OCCT TopExp_Explorer(S, TopAbs_VERTEX) over the face's edge vertices.
pub fn topexp_face_vertices(brep: &topods::BRep, f: &Shape) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    for e in topexp_face_edges(brep, f) {
        let ed = e.as_edge().expect("not an edge");
        for v in [&ed.first, &ed.last] {
            if !out.iter().any(|o| o.is_same(v)) {
                out.push(v.clone());
            }
        }
    }
    for v in &f.as_face().expect("not a face").internal_vertices {
        if !out.iter().any(|o| o.is_same(v)) {
            out.push(v.clone());
        }
    }
    out
}

/// OCCT Bnd_Box (Bnd_Box.hxx) — the min/max box the corner machinery
/// accumulates points into.
#[derive(Debug, Clone)]
pub struct BndBox {
    pub is_void: bool,
    pub min: DVec3,
    pub max: DVec3,
}

impl Default for BndBox {
    fn default() -> Self {
        BndBox {
            is_void: true,
            min: DVec3::ZERO,
            max: DVec3::ZERO,
        }
    }
}

impl BndBox {
    /// OCCT Bnd_Box::Add(P).
    pub fn add(&mut self, p: DVec3) {
        if self.is_void {
            self.is_void = false;
            self.min = p;
            self.max = p;
        } else {
            self.min = self.min.min(p);
            self.max = self.max.max(p);
        }
    }

    /// OCCT Bnd_Box::Get(...) — the corner coordinates.
    pub fn get(&self) -> (f64, f64, f64, f64, f64, f64) {
        (
            self.min.x,
            self.min.y,
            self.min.z,
            self.max.x,
            self.max.y,
            self.max.z,
        )
    }
}

// =========================================================================
// OCCT GeomAbs_SurfaceType (TKG3d/GeomAbs) — the surface kinds the builder
// switches on.
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BSplineSurface,
    BezierSurface,
    Other,
}

/// Natural (underlying) bounds and periodicity per OCCT Geom_* analytic
/// surfaces (Geom_CylindricalSurface::Bounds(0, 2PI, -inf, +inf) etc.).
fn surface_natural_bounds(s: &Surface3) -> ([f64; 4], bool, bool, f64, f64) {
    match s {
        Surface3::Plane(_) => (
            [f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY],
            false,
            false,
            0.0,
            0.0,
        ),
        Surface3::Cylinder(_) => ([0.0, 2.0 * std::f64::consts::PI, f64::NEG_INFINITY, f64::INFINITY], true, false, 2.0 * std::f64::consts::PI, 0.0),
        Surface3::Cone(_) => ([0.0, 2.0 * std::f64::consts::PI, f64::NEG_INFINITY, f64::INFINITY], true, false, 2.0 * std::f64::consts::PI, 0.0),
        Surface3::Sphere(_) => (
            [0.0, 2.0 * std::f64::consts::PI, -0.5 * std::f64::consts::PI, 0.5 * std::f64::consts::PI],
            true,
            false,
            2.0 * std::f64::consts::PI,
            0.0,
        ),
        Surface3::Torus(_) => (
            [0.0, 2.0 * std::f64::consts::PI, 0.0, 2.0 * std::f64::consts::PI],
            true,
            true,
            2.0 * std::f64::consts::PI,
            2.0 * std::f64::consts::PI,
        ),
        _ => ([0.0, 1.0, 0.0, 1.0], false, false, 0.0, 0.0),
    }
}

pub fn surface_type_of(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::BSpline(_) => GeomAbsSurfaceType::BSplineSurface,
        Surface3::Bezier(_) => GeomAbsSurfaceType::BezierSurface,
        _ => GeomAbsSurfaceType::Other,
    }
}

// =========================================================================
// OCCT GeomAdaptor_Surface (TKG3d) — the trimmed view of a DS surface used
// by ChFi3d_BoundSurf / ChFi3d_ComputeCurves.
// =========================================================================
#[derive(Debug, Clone)]
pub struct GeomAdaptorSurface {
    pub surface: Surface3,
    pub ufirst: f64,
    pub ulast: f64,
    pub vfirst: f64,
    pub vlast: f64,
    /// False while the adaptor only carries the surface (OCCT distinguishes
    /// the two Load() overloads).
    pub bounds_set: bool,
}

impl GeomAdaptorSurface {
    /// OCCT GeomAdaptor_Surface::Load(S).
    pub fn new(surface: Surface3) -> Self {
        let (b, _, _, _, _) = surface_natural_bounds(&surface);
        GeomAdaptorSurface {
            surface,
            ufirst: b[0],
            ulast: b[1],
            vfirst: b[2],
            vlast: b[3],
            bounds_set: false,
        }
    }

    /// OCCT GeomAdaptor_Surface::Load(S, U1, U2, V1, V2).
    pub fn load_bounded(&mut self, surface: Surface3, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.surface = surface;
        self.ufirst = u1;
        self.ulast = u2;
        self.vfirst = v1;
        self.vlast = v2;
        self.bounds_set = true;
    }

    pub fn surface(&self) -> &Surface3 {
        &self.surface
    }

    pub fn get_type(&self) -> GeomAbsSurfaceType {
        surface_type_of(&self.surface)
    }

    /// OCCT Adaptor3d_Surface::Value(U, V).
    pub fn value(&self, u: f64, v: f64) -> DVec3 {
        use rcad_kernel::geom::SurfaceEval as _;
        self.surface.point_at(u, v)
    }

    pub fn first_u_parameter(&self) -> f64 {
        self.ufirst
    }

    pub fn last_u_parameter(&self) -> f64 {
        self.ulast
    }

    pub fn first_v_parameter(&self) -> f64 {
        self.vfirst
    }

    pub fn last_v_parameter(&self) -> f64 {
        self.vlast
    }

    pub fn is_u_periodic(&self) -> bool {
        surface_natural_bounds(&self.surface).1
    }

    pub fn is_v_periodic(&self) -> bool {
        surface_natural_bounds(&self.surface).2
    }

    pub fn u_period(&self) -> f64 {
        surface_natural_bounds(&self.surface).3
    }

    pub fn v_period(&self) -> f64 {
        surface_natural_bounds(&self.surface).4
    }

    /// OCCT GeomAdaptor_Surface::UResolution (GeomAdaptor_Surface.cxx
    /// L1818-1875): plane T, cylinder/sphere R3d/(2R), cone R3d/Rmax,
    /// torus R3d/(2(Rmaj+Rmin)).
    pub fn u_resolution(&self, r3d: f64) -> f64 {
        match &self.surface {
            Surface3::Plane(_) => r3d,
            Surface3::Cylinder(c) => {
                if c.radius > CONFUSION {
                    r3d / (2.0 * c.radius)
                } else {
                    0.0
                }
            }
            Surface3::Sphere(s) => {
                if s.radius > CONFUSION {
                    r3d / (2.0 * s.radius)
                } else {
                    0.0
                }
            }
            Surface3::Cone(c) => {
                // OCCT uses the VIso radii at VLast/VFirst; the adaptor
                // carries the current (bounded) V range.
                let r1 = c.radius_at_slant(self.vlast.max(self.vfirst));
                let r2 = c.radius_at_slant(self.vfirst.min(self.vlast));
                let r = r1.max(r2);
                if r > CONFUSION {
                    r3d / r
                } else {
                    0.0
                }
            }
            Surface3::Torus(t) => {
                let r = t.major_radius + t.minor_radius;
                if r > CONFUSION {
                    r3d / (2.0 * r)
                } else {
                    0.0
                }
            }
            _ => r3d,
        }
    }

    /// OCCT VResolution — plane T, cylinder T, cone T, sphere R3d/(2R),
    /// torus R3d/(2*minor).
    pub fn v_resolution(&self, r3d: f64) -> f64 {
        match &self.surface {
            Surface3::Plane(_) => r3d,
            Surface3::Cylinder(_) => r3d,
            Surface3::Cone(_) => r3d,
            Surface3::Sphere(s) => {
                if s.radius > CONFUSION {
                    r3d / (2.0 * s.radius)
                } else {
                    0.0
                }
            }
            Surface3::Torus(t) => {
                if t.minor_radius > CONFUSION {
                    r3d / (2.0 * t.minor_radius)
                } else {
                    0.0
                }
            }
            _ => r3d,
        }
    }

    /// OCCT GeomAdaptor_Surface::Plane().
    pub fn plane(&self) -> rcad_kernel::geom::Plane {
        match &self.surface {
            Surface3::Plane(p) => p.clone(),
            _ => panic!("GeomAdaptor_Surface::Plane - not a plane"),
        }
    }

    /// OCCT GeomAdaptor_Surface::Cylinder().
    pub fn cylinder(&self) -> rcad_kernel::geom::CylindricalSurface {
        match &self.surface {
            Surface3::Cylinder(c) => c.clone(),
            _ => panic!("GeomAdaptor_Surface::Cylinder - not a cylinder"),
        }
    }
}

// =========================================================================
// OCCT BRepAdaptor_Surface (TKTopAlgo) — a face with its world surface and
// the UV bounds (Trim1/Trim2 equivalents used by BoundFac).
// =========================================================================
#[derive(Debug, Clone)]
pub struct BRepAdaptorSurface {
    pub brep: topods::BRep,
    pub face: Shape,
    pub surface: Surface3,
    pub ufirst: f64,
    pub ulast: f64,
    pub vfirst: f64,
    pub vlast: f64,
    pub bounds_set: bool,
}

impl BRepAdaptorSurface {
    /// OCCT BRepAdaptor_Surface::Initialize(F, Restriction=false).
    pub fn initialize(brep: &topods::BRep, face: &Shape) -> Self {
        let surface = brep
            .face_surface_world(face)
            .expect("BRepAdaptor_Surface: face without surface");
        let (b, _, _, _, _) = surface_natural_bounds(&surface);
        BRepAdaptorSurface {
            brep: brep.clone(),
            face: face.clone(),
            surface,
            ufirst: b[0],
            ulast: b[1],
            vfirst: b[2],
            vlast: b[3],
            bounds_set: false,
        }
    }

    pub fn surface(&self) -> &Surface3 {
        &self.surface
    }

    /// OCCT L834-851: BRE.MakeFace(FFv, Sface, tol) + Bs.Initialize(FFv,
    /// false) builds a bare face over the (possibly extended) surface and
    /// loads the adaptor from it.  rcad carries the surface directly on the
    /// adaptor; the bare-face construction collapses into this constructor
    /// (architecture note, chfi3d_builder_c1).
    pub fn initialize_surface(surface: Surface3) -> Self {
        let (b, _, _, _, _) = surface_natural_bounds(&surface);
        BRepAdaptorSurface {
            brep: topods::BRep::default(),
            face: Shape::null(),
            surface,
            ufirst: b[0],
            ulast: b[1],
            vfirst: b[2],
            vlast: b[3],
            bounds_set: false,
        }
    }

    /// OCCT BRepAdaptor_Surface::GeomSurfaceOriginal()/FirstUParameter...
    pub fn value(&self, u: f64, v: f64) -> DVec3 {
        use rcad_kernel::geom::SurfaceEval as _;
        self.surface.point_at(u, v)
    }

    pub fn get_type(&self) -> GeomAbsSurfaceType {
        surface_type_of(&self.surface)
    }

    pub fn is_u_periodic(&self) -> bool {
        surface_natural_bounds(&self.surface).1
    }

    pub fn is_v_periodic(&self) -> bool {
        surface_natural_bounds(&self.surface).2
    }

    pub fn u_period(&self) -> f64 {
        surface_natural_bounds(&self.surface).3
    }

    pub fn v_period(&self) -> f64 {
        surface_natural_bounds(&self.surface).4
    }

    /// OCCT Adaptor3d_Surface::UResolution (per-kind formulas above).
    pub fn u_resolution(&self, r3d: f64) -> f64 {
        let g = GeomAdaptorSurface::new(self.surface.clone());
        g.u_resolution(r3d)
    }

    pub fn v_resolution(&self, r3d: f64) -> f64 {
        let g = GeomAdaptorSurface::new(self.surface.clone());
        g.v_resolution(r3d)
    }

    /// OCCT BRepAdaptor_Surface::Load(S, U1, U2, V1, V2, Trsf, TolU, TolV)
    /// via ChFi3d_ApplyBounds — the (possibly trimmed-basis) surface and
    /// bounds are reloaded.
    pub fn apply_bounds(&mut self, surface: Surface3, u1: f64, u2: f64, v1: f64, v2: f64) {
        self.surface = surface;
        self.ufirst = u1;
        self.ulast = u2;
        self.vfirst = v1;
        self.vlast = v2;
        self.bounds_set = true;
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L242-253 — ChFi3d_Boite (2-point variant).
// =========================================================================
pub fn chfi3d_boite(p1: DVec2, p2: DVec2) -> (f64, f64, f64, f64) {
    let mu = p1.x.min(p2.x);
    let big_m = p1.x.max(p2.x);
    let mv = p1.y.min(p2.y);
    let big_mv = p1.y.max(p2.y);
    (mu, big_m, mv, big_mv)
}

/// OCCT ChFi3d_Builder_0.cxx L259-285 — ChFi3d_Boite (4-point variant).
#[allow(clippy::too_many_arguments)]
pub fn chfi3d_boite4(p1: DVec2, p2: DVec2, p3: DVec2, p4: DVec2) -> (f64, f64, f64, f64, f64, f64) {
    let a = p1.x.min(p2.x);
    let b = p3.x.min(p4.x);
    let mu = a.min(b);
    let a = p1.x.max(p2.x);
    let b = p3.x.max(p4.x);
    let big_m = a.max(b);
    let a = p1.y.min(p2.y);
    let b = p3.y.min(p4.y);
    let mv = a.min(b);
    let a = p1.y.max(p2.y);
    let b = p3.y.max(p4.y);
    let big_mv = a.max(b);
    let du = big_m - mu;
    let dv = big_mv - mv;
    (du, dv, mu, big_m, mv, big_mv)
}

/// OCCT ChFi3d_Builder_0.cxx L289-315 — Geometry(DStr, ind): the adaptor of
/// the shape (ind > 0) or surface (ind < 0) at index |ind|.  rcad passes
/// the owning BRep for the face branch (OCCT reads geometry from the
/// TopoDS handle graph).
pub fn chfi3d_geometry(
    brep: &topods::BRep,
    dstr: &TopOpeBRepDSHDataStructure,
    ind: i32,
) -> Option<GeomAdaptorSurface> {
    if ind == 0 {
        return None;
    }
    if ind > 0 {
        let f = dstr.shape(ind);
        if f.is_null() {
            return None;
        }
        let hs = BRepAdaptorSurface::initialize(brep, f);
        Some(GeomAdaptorSurface::new(hs.surface.clone()))
    } else {
        let s = dstr.surface(-ind);
        Some(GeomAdaptorSurface::new(s.surface.clone()))
    }
}

/// OCCT ChFi3d_Builder_0.cxx L319-331 — ChFi3d_SetPointTolerance.
pub fn chfi3d_set_point_tolerance(dstr: &mut TopOpeBRepDSHDataStructure, bx: &BndBox, ip: i32) {
    let (a, b, c, d, e, f) = bx.get();
    let mut d = d - a;
    let mut e = e - b;
    let mut f = f - c;
    d *= d;
    e *= e;
    f *= f;
    let vtol = (d + e + f).sqrt() * 1.5;
    dstr.change_point(ip).set_tolerance(vtol);
}

/// OCCT ChFi3d_Builder_0.cxx L335-343 — ChFi3d_EnlargeBox(C, wd, wf, ...).
pub fn chfi3d_enlarge_box_curve(
    c: &rcad_kernel::geom::Curve3,
    wd: f64,
    wf: f64,
    box1: &mut BndBox,
    box2: &mut BndBox,
) {
    use rcad_kernel::geom::CurveEval as _;
    box1.add(c.point_at(wd));
    box2.add(c.point_at(wf));
}

/// OCCT ChFi3d_Builder_0.cxx L347-359 — ChFi3d_EnlargeBox(S, PC, wd, wf, ...).
pub fn chfi3d_enlarge_box_surf_pc(
    s: &GeomAdaptorSurface,
    pc: &rcad_kernel::geom::Curve2d,
    wd: f64,
    wf: f64,
    box1: &mut BndBox,
    box2: &mut BndBox,
) {
    use rcad_kernel::geom::Curve2dEval as _;
    let uv = pc.point_at(wd);
    box1.add(s.value(uv.x, uv.y));
    let uv = pc.point_at(wf);
    box2.add(s.value(uv.x, uv.y));
}

/// OCCT ChFi3d_Builder_0.cxx L363-381 — ChFi3d_EnlargeBox(E, LF, w, box).
pub fn chfi3d_enlarge_box_edge_faces(
    brep: &topods::BRep,
    e: &Shape,
    lf: &[Shape],
    w: f64,
    bx: &mut BndBox,
) {
    use rcad_kernel::geom::CurveEval as _;
    let ed = e.as_edge().expect("not an edge");
    let Some(curve) = ed.curve.as_ref() else {
        return;
    };
    bx.add(curve.point_at(w));
    for f in lf {
        if f.is_null() {
            continue;
        }
        if let Some((pc, f_, l_)) = brep.curve_on_surface(e, f) {
            // OCCT: BC.Initialize(E, F); box.Add(BC.Value(w)) — the point
            // via the face pcurve parameterization.
            let _ = (f_, l_);
            let uv = pc.point_at(w);
            if let Some(s) = brep.face_surface_world(f) {
                use rcad_kernel::geom::SurfaceEval as _;
                bx.add(s.point_at(uv.x, uv.y));
            }
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L157-235 — ChFi3d_BoundSrfImpl + L691-699
// ChFi3d_BoundFac: resize the face adaptor bounds around [umin..umax x
// vmin..vmax].  rcad: one implementation per adaptor kind.
// =========================================================================

/// OCCT ChFi3d_BoundSrfImpl (Builder_0.cxx L157-235) over a face adaptor.
pub fn chfi3d_bound_fac(
    s: &mut BRepAdaptorSurface,
    uumin: f64,
    uumax: f64,
    vvmin: f64,
    vvmax: f64,
    checknaturalbounds: bool,
) {
    let mut a_umin = uumin;
    let mut a_umax = uumax;
    let mut a_vmin = vvmin;
    let mut a_vmax = vvmax;
    let mut a_surface = s.surface.clone();
    if let Surface3::Trimmed(tr) = &a_surface {
        a_surface = *tr.basis.clone();
    }

    let (nb, uper_u, vper_v, period_u, period_v) = surface_natural_bounds(&a_surface);
    let (a_u1, a_u2, a_v1, a_v2) = (nb[0], nb[1], nb[2], nb[3]);

    let mut a_step_u = a_umax - a_umin;
    let mut a_step_v = a_vmax - a_vmin;
    let a_scale_u = s.u_resolution(1.0);
    let a_scale_v = s.v_resolution(1.0);
    let a_step3d_u = a_step_u / a_scale_u;
    let a_step3d_v = a_step_v / a_scale_v;

    if a_step3d_u > a_step3d_v {
        a_step_v = a_step3d_u * a_scale_v;
    }
    if a_step3d_v > a_step3d_u {
        a_step_u = a_step3d_v * a_scale_u;
    }

    if period_u > 0.0 {
        a_step_u = 0.1 * (period_u - (a_umax - a_umin));
    }
    if period_v > 0.0 {
        a_step_v = 0.1 * (period_v - (a_vmax - a_vmin));
    }

    let mut a_uu1 = a_umin - a_step_u;
    let mut a_uu2 = a_umax + a_step_u;
    let mut a_vv1 = a_vmin - a_step_v;
    let mut a_vv2 = a_vmax + a_step_v;
    if checknaturalbounds {
        if !uper_u {
            a_uu1 = a_uu1.max(a_u1);
            a_uu2 = a_uu2.min(a_u2);
        }
        if !vper_v {
            a_vv1 = a_vv1.max(a_v1);
            a_vv2 = a_vv2.min(a_v2);
        }
    }

    s.apply_bounds(a_surface, a_uu1, a_uu2, a_vv1, a_vv2);
}

/// OCCT ChFi3d_BoundSrf (Builder_0.cxx L706-714) over a GeomAdaptor.
pub fn chfi3d_bound_srf(
    s: &mut GeomAdaptorSurface,
    uumin: f64,
    uumax: f64,
    vvmin: f64,
    vvmax: f64,
    checknaturalbounds: bool,
) {
    let mut a_umin = uumin;
    let mut a_umax = uumax;
    let mut a_vmin = vvmin;
    let mut a_vmax = vvmax;
    let mut a_surface = s.surface.clone();
    if let Surface3::Trimmed(tr) = &a_surface {
        a_surface = *tr.basis.clone();
    }

    let (nb, uper_u, vper_v, period_u, period_v) = surface_natural_bounds(&a_surface);
    let (a_u1, a_u2, a_v1, a_v2) = (nb[0], nb[1], nb[2], nb[3]);

    let mut a_step_u = a_umax - a_umin;
    let mut a_step_v = a_vmax - a_vmin;
    let a_scale_u = s.u_resolution(1.0);
    let a_scale_v = s.v_resolution(1.0);
    let a_step3d_u = a_step_u / a_scale_u;
    let a_step3d_v = a_step_v / a_scale_v;

    if a_step3d_u > a_step3d_v {
        a_step_v = a_step3d_u * a_scale_v;
    }
    if a_step3d_v > a_step3d_u {
        a_step_u = a_step3d_v * a_scale_u;
    }

    if period_u > 0.0 {
        a_step_u = 0.1 * (period_u - (a_umax - a_umin));
    }
    if period_v > 0.0 {
        a_step_v = 0.1 * (period_v - (a_vmax - a_vmin));
    }

    let mut a_uu1 = a_umin - a_step_u;
    let mut a_uu2 = a_umax + a_step_u;
    let mut a_vv1 = a_vmin - a_step_v;
    let mut a_vv2 = a_vmax + a_step_v;
    if checknaturalbounds {
        if !uper_u {
            a_uu1 = a_uu1.max(a_u1);
            a_uu2 = a_uu2.min(a_u2);
        }
        if !vper_v {
            a_vv1 = a_vv1.max(a_v1);
            a_vv2 = a_vv2.min(a_v2);
        }
    }

    s.load_bounded(a_surface, a_uu1, a_uu2, a_vv1, a_vv2);
}

/// OCCT ChFi3d_Builder_0.cxx L4453-4505 — ChFi3d_BoundSurf: the adaptor of
/// the SurfData surface trimmed to the interference pcurve box.
pub fn chfi3d_bound_surf(
    dstr: &TopOpeBRepDSHDataStructure,
    fd1: &ChFiDSSurfData,
    ifaco1: i32,
    ifaarc1: i32,
) -> GeomAdaptorSurface {
    // rmq : as in fact 2 interferences of Fd1 serve only to set limits,
    // indexes IFaCo1 and IFaArc1 are not useful (kept as an option).
    let mut hs1 = GeomAdaptorSurface::new(dstr.surface(fd1.surf()).surface.clone());

    if ifaco1 == 0 || ifaarc1 == 0 {
        return hs1;
    }

    let fi_co1 = fd1.interference(ifaco1);
    let fi_arc1 = fd1.interference(ifaarc1);

    use rcad_kernel::geom::Curve2dEval as _;
    let uvf1 = fi_co1.pcurve_on_surf().map(|pc| pc.point_at(fi_co1.parameter_first())).unwrap_or(DVec2::ZERO);
    let uvl1 = fi_co1.pcurve_on_surf().map(|pc| pc.point_at(fi_co1.parameter_last())).unwrap_or(DVec2::ZERO);
    let uvf2 = fi_arc1.pcurve_on_surf().map(|pc| pc.point_at(fi_arc1.parameter_first())).unwrap_or(DVec2::ZERO);
    let uvl2 = fi_arc1.pcurve_on_surf().map(|pc| pc.point_at(fi_arc1.parameter_last())).unwrap_or(DVec2::ZERO);
    let (du0, dv0, mu, big_mu, mv, big_mv) = chfi3d_boite4(uvf1, uvf2, uvl1, uvl2);
    let styp = hs1.get_type();
    match styp {
        GeomAbsSurfaceType::Cylinder => {
            let radius = match &hs1.surface {
                rcad_kernel::geom::Surface3::Cylinder(c) => c.radius,
                _ => 0.0,
            };
            let mut dv = dv0;
            dv = (0.5 * dv).max(4.0 * radius);
            let du = 0.0;
            let surf = dstr.surface(fd1.surf()).surface.clone();
            hs1.load_bounded(surf, mu, big_mu, mv - dv, big_mv + dv);
            let _ = du;
        }
        GeomAbsSurfaceType::Torus | GeomAbsSurfaceType::Cone => {
            let du = (std::f64::consts::PI - 0.5 * du0).min(0.1 * du0);
            let dv = 0.0;
            let surf = dstr.surface(fd1.surf()).surface.clone();
            hs1.load_bounded(surf, mu - du, big_mu + du, mv, big_mv);
            let _ = dv;
        }
        GeomAbsSurfaceType::Plane => {
            let du = (0.5 * du0).max(4.0 * dv0);
            let dv = 0.0;
            let surf = dstr.surface(fd1.surf()).surface.clone();
            hs1.load_bounded(surf, mu - du, big_mu + du, mv, big_mv);
            let _ = dv;
        }
        _ => {}
    }
    hs1
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L5750-5769 — ChFi3d_Couture: determine if F has
// a sewing (seam) edge.
// =========================================================================
pub fn chfi3d_couture(brep: &topods::BRep, f: &Shape) -> (bool, Shape) {
    let mut couture = false;
    let mut edgecouture = Shape::null();
    for ecur in topexp_face_edges(brep, f) {
        if couture {
            break;
        }
        if brep.is_edge_closed_on_face(&ecur, f) {
            couture = true;
            edgecouture = ecur;
        }
    }
    (couture, edgecouture)
}

/// OCCT ChFi3d_Builder_0.cxx L5773-5801 — ChFi3d_CoutureOnVertex.
pub fn chfi3d_couture_on_vertex(brep: &topods::BRep, f: &Shape, v: &Shape) -> (bool, Shape) {
    let mut couture = false;
    let mut edgecouture = Shape::null();
    for ecur in topexp_face_edges(brep, f) {
        if brep.is_edge_closed_on_face(&ecur, f) {
            let ed = ecur.as_edge().expect("not an edge");
            if ed.first.is_same(v) || ed.last.is_same(v) {
                couture = true;
                edgecouture = ecur;
                break;
            }
        }
    }
    (couture, edgecouture)
}

/// OCCT BRepTools::IsReallyClosed(E, F) — the edge occurs twice among the
/// face's wires.
fn brep_tools_is_really_closed(brep: &topods::BRep, e: &Shape, f: &Shape) -> bool {
    let edges = topexp_face_edges(brep, f);
    let mut n = 0usize;
    for we in &edges {
        if we.is_same(e) {
            n += 1;
        }
    }
    n > 1
}

/// OCCT ChFi3d_Builder_0.cxx L5805-5831 — ChFi3d_IsPseudoSeam: a closed
/// edge whose both vertices also touch another really-closed edge of F.
pub fn chfi3d_is_pseudo_seam(brep: &topods::BRep, e: &Shape, f: &Shape) -> bool {
    if !brep.is_edge_closed_on_face(e, f) {
        return false;
    }
    let mut neighbor_seam_found = false;
    let ed = e.as_edge().expect("not an edge");
    let (vf, vl) = (ed.first.clone(), ed.last.clone());
    for ecur in topexp_face_edges(brep, f) {
        if !ecur.is_same(e) {
            let edc = ecur.as_edge().expect("not an edge");
            let (v1, v2) = (edc.first.clone(), edc.last.clone());
            if (v1.is_same(&vf) || v1.is_same(&vl) || v2.is_same(&vf) || v2.is_same(&vl))
                && brep_tools_is_really_closed(brep, &ecur, f)
            {
                neighbor_seam_found = true;
                break;
            }
        }
    }
    neighbor_seam_found
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2156-2170 — ChFi3d_FilCurveInDS.
// =========================================================================
pub fn chfi3d_fil_curve_in_ds(
    icurv: i32,
    isurf: i32,
    pcurv: Option<rcad_kernel::geom::Curve2d>,
    et: Orientation,
) -> TopOpeBRepDSInterference {
    // OCCT: SurfaceCurveInterference(Transition(Et), SURFACE, Isurf,
    //       CURVE, Icurv, Pcurv).  The support keeps the TopOpeBRepDS_SURFACE
    // kind even when the index designates a shape (OCCT callers pass face
    // shape indices through this same constructor).
    TopOpeBRepDSInterference::SurfaceCurve(TopOpeBRepDSSurfaceCurveInterference::new(
        et,
        TopOpeBRepDSKind::Surface,
        isurf,
        TopOpeBRepDSKind::Curve,
        icurv,
        pcurv,
    ))
}

/// OCCT ChFi3d_Builder_0.cxx L2341-2367 — ChFi3d_FilPointInDS.
pub fn chfi3d_fil_point_in_ds(
    et: Orientation,
    ic: i32,
    ip: i32,
    par: f64,
    is_vertex: bool,
) -> TopOpeBRepDSInterference {
    let kind_g = if is_vertex {
        TopOpeBRepDSKind::Vertex
    } else {
        TopOpeBRepDSKind::Point
    };
    TopOpeBRepDSInterference::CurvePoint(TopOpeBRepDSCurvePointInterference::new(
        et,
        TopOpeBRepDSKind::Curve,
        ic,
        kind_g,
        ip,
        par,
    ))
}

/// OCCT ChFi3d_Builder_0.cxx L2371-2385 — ChFi3d_FilVertexInDS.
pub fn chfi3d_fil_vertex_in_ds(
    et: Orientation,
    ic: i32,
    ip: i32,
    par: f64,
) -> TopOpeBRepDSInterference {
    TopOpeBRepDSInterference::CurvePoint(TopOpeBRepDSCurvePointInterference::new(
        et,
        TopOpeBRepDSKind::Curve,
        ic,
        TopOpeBRepDSKind::Vertex,
        ip,
        par,
    ))
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L3540-3547 — ChFi3d_ConvTol2dToTol3d.
// =========================================================================
pub fn chfi3d_conv_tol2d_to_tol3d(s: &GeomAdaptorSurface, tol2d: f64) -> f64 {
    let ures = s.u_resolution(1.0e-7);
    let vres = s.v_resolution(1.0e-7);
    let uresto3d = 1.0e-7 * tol2d / ures;
    let vresto3d = 1.0e-7 * tol2d / vres;
    uresto3d.max(vresto3d)
}

/// OCCT ChFi3d_Builder_0.cxx L3555-3596 — ChFi3d_EvalTolReached.
pub fn chfi3d_eval_tol_reached(
    s1: &GeomAdaptorSurface,
    pc1: &rcad_kernel::geom::Curve2d,
    s2: &GeomAdaptorSurface,
    pc2: &rcad_kernel::geom::Curve2d,
    c: &rcad_kernel::geom::Curve3,
) -> f64 {
    use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _};
    let mut distmax = 0.0f64;
    let [f, l] = c.default_domain();
    let nbp = 45usize;
    let step = 1.0 / (nbp as f64 - 1.0);
    for i in 0..nbp {
        let mut t = step * i as f64;
        t = (1.0 - t) * f + t * l;
        let uv = pc1.point_at(t);
        let ps1 = s1.value(uv.x, uv.y);
        let uv = pc2.point_at(t);
        let ps2 = s2.value(uv.x, uv.y);
        let pc = c.point_at(t);
        let mut d = (ps1 - pc).length_squared();
        if d > distmax {
            distmax = d;
        }
        d = (ps2 - pc).length_squared();
        if d > distmax {
            distmax = d;
        }
        d = (ps1 - ps2).length_squared();
        if d > distmax {
            distmax = d;
        }
    }
    let distmax = 1.5 * distmax.sqrt();
    distmax.max(CONFUSION)
}

/// OCCT ChFi3d_Builder_0.cxx L385-469 — ChFi3d_EnlargeBox(DStr, st, sd,
/// b1, b2, isfirst).
pub fn chfi3d_enlarge_box_dstr(
    brep: &topods::BRep,
    dstr: &TopOpeBRepDSHDataStructure,
    st: Option<&ChFiDSStripe>,
    sd: &ChFiDSSurfData,
    b1: &mut BndBox,
    b2: &mut BndBox,
    isfirst: bool,
) {
    use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _};
    let cp1 = sd.vertex(isfirst, 1);
    let cp2 = sd.vertex(isfirst, 2);
    b1.add(cp1.point());
    b2.add(cp2.point());
    let fi1 = sd.interference_on_s1();
    let fi2 = sd.interference_on_s2();
    let s = &dstr.surface(sd.surf()).surface;
    let pcs1 = fi1.pcurve_on_surf();
    let pcs2 = fi2.pcurve_on_surf();
    let c3d1 = dstr
        .curve(fi1.line_index())
        .curve
        .clone();
    let c3d2 = dstr
        .curve(fi2.line_index())
        .curve
        .clone();
    let f1 = chfi3d_geometry(brep, dstr, sd.index_of_s1);
    let f2 = chfi3d_geometry(brep, dstr, sd.index_of_s2);
    let p1 = fi1.parameter(isfirst);
    if let Some(c3d1) = &c3d1 {
        b1.add(c3d1.point_at(p1));
    }
    if let Some(pcs1) = pcs1 {
        let uv = pcs1.point_at(p1);
        let sv = GeomAdaptorSurface::new(s.clone()).value(uv.x, uv.y);
        b1.add(sv);
    }
    if let Some(f1) = &f1 {
        let pcf1 = fi1.pcurve_on_face();
        if let Some(pcf1) = pcf1 {
            let uv = pcf1.point_at(p1);
            b1.add(f1.value(uv.x, uv.y));
        }
    }
    let p2 = fi2.parameter(isfirst);
    if let Some(c3d2) = &c3d2 {
        b2.add(c3d2.point_at(p2));
    }
    if let Some(pcs2) = pcs2 {
        let uv = pcs2.point_at(p2);
        let sv = GeomAdaptorSurface::new(s.clone()).value(uv.x, uv.y);
        b2.add(sv);
    }
    if let Some(f2) = &f2 {
        let pcf2 = fi2.pcurve_on_face();
        if let Some(pcf2) = pcf2 {
            let uv = pcf2.point_at(p2);
            b2.add(f2.value(uv.x, uv.y));
        }
    }
    if let Some(st) = st {
        let (icurv, ipcurve, orint, (pa, pb)) = if isfirst {
            (
                st.first_curve(),
                st.first_pcurve(),
                st.orientation(true),
                st.first_parameters(),
            )
        } else {
            (
                st.last_curve(),
                st.last_pcurve(),
                st.orientation(false),
                st.last_parameters(),
            )
        };
        let c3d = dstr.curve(icurv).curve.clone();
        let c2d = ipcurve.cloned();
        let (mut p1, mut p2) = (pa, pb);
        if orint != Orientation::Forward {
            p2 = pa;
            p1 = pb;
        }
        if let Some(c3d) = &c3d {
            b1.add(c3d.point_at(p1));
            b2.add(c3d.point_at(p2));
        }
        if let Some(c2d) = &c2d {
            let uv = c2d.point_at(p1);
            let sv = GeomAdaptorSurface::new(s.clone()).value(uv.x, uv.y);
            b1.add(sv);
            let uv = c2d.point_at(p2);
            let sv = GeomAdaptorSurface::new(s.clone()).value(uv.x, uv.y);
            b2.add(sv);
        }
    }
}

/// OCCT ChFi3d_Builder_0.cxx L1472-1507 — ChFi3d_ReparamPcurv.  The BSpline
/// branch reparameterizes the knots; non-BSpline pcurves pass unchanged
/// (OCCT returns them untouched as well).
pub fn chfi3d_reparam_pcurv(uf: f64, ul: f64, pcurv: rcad_kernel::geom::Curve2d) -> rcad_kernel::geom::Curve2d {
    match pcurv {
        rcad_kernel::geom::Curve2d::BSpline(mut pc) => {
            // BSpline parametric bounds: knots[degree] .. knots[nb-1-degree].
            let upcf = pc.knots[pc.degree];
            let upcl = pc.knots[pc.knots.len() - 1 - pc.degree];
            // OCCT L1491-1505: Segment on the trimmed range, then
            // BSplCLib::Reparametrize over the knots to [Uf, Ul].  rcad
            // reparameterizes by a linear knot map; when the BSpline already
            // spans [Uf, Ul] nothing is done.
            if (uf - upcf).abs() > P_CONFUSION || (ul - upcl).abs() > P_CONFUSION {
                let span = upcl - upcf;
                if span.abs() > P_CONFUSION {
                    let a = (ul - uf) / span;
                    let b = uf - a * upcf;
                    for k in pc.knots.iter_mut() {
                        *k = a * *k + b;
                    }
                }
            }
            rcad_kernel::geom::Curve2d::BSpline(pc)
        }
        rcad_kernel::geom::Curve2d::Trimmed(_) => {
            // OCCT L1481-1485: a trimmed pcurve is unwrapped to its basis
            // before the BSpline test; rcad TrimmedCurve2 carries the basis.
            pcurv
        }
        other => other,
    }
}

/// OCCT ChFi3d_Builder_0.cxx L1514-1559 — ChFi3d_ProjectPCurv: the pcurve
/// of the intersection line on an analytic surface via ProjLib.  Returns
/// (Pcurv, tolreached); None encodes the OCCT non-analytic no-op (the
/// handle stays null) / the raised NotImplenented for the default case.
pub fn chfi3d_project_pcurv(
    hcg: &rcad_kernel::geom::Curve3,
    hsg: &GeomAdaptorSurface,
    tol: f64,
) -> Option<(rcad_kernel::geom::Curve2d, f64)> {
    use rcad_kernel::base::proj_lib::project_on_surface;
    let _ = tol;
    if hsg.get_type() != GeomAbsSurfaceType::BezierSurface
        && hsg.get_type() != GeomAbsSurfaceType::BSplineSurface
    {
        let pcurv = project_on_surface(hcg, &hsg.surface)?;
        // OCCT reads Projc.GetTolerance(); the rcad projector does not
        // track it — the analytic projections are exact, CONFUSION stands.
        let tolreached = P_CONFUSION;
        Some((pcurv, tolreached))
    } else {
        None
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L3694-4294 — ChFi3d_ComputeCurves.
//
// Calculates the intersection between two surfaces with known extremities
// (Pardeb/Parfin 4-tuples: S1 uv then S2 uv).  The analytic branches
// (cylinder/plane and plane/plane via IntAna_QuadQuadGeo) are translated
// 1:1; the generic surface branch (GeomInt_IntSS after trsfsurf, Builder_0
// L3888-4070) and the walking fallback (IntWalk_PWalking + GeomInt_WLApprox,
// L4075-4293) are pending TKGeomAlgo/TKTopAlgo translations and report the
// OCCT failure path (return None).
// =========================================================================

/// The intersection result of ChFi3d_ComputeCurves (the C3d/Pc1/Pc2 handle
/// out-parameters plus tolreached).
#[derive(Debug, Clone)]
pub struct ComputedCurves {
    pub c3d: rcad_kernel::geom::Curve3,
    pub pc1: rcad_kernel::geom::Curve2d,
    pub pc2: rcad_kernel::geom::Curve2d,
    pub tolreached: f64,
}

/// OCCT ElCLib::Parameter(L, P).
pub fn elclib_parameter_line(l: &rcad_kernel::geom::Line3, p: DVec3) -> f64 {
    (p - l.origin).dot(l.direction)
}

/// OCCT ElCLib::Value(U, L).
pub fn elclib_line_value(u: f64, l: &rcad_kernel::geom::Line3) -> DVec3 {
    l.origin + l.direction * u
}

/// OCCT ElCLib::Parameter(C, P) — the angle of the projected point.
pub fn elclib_parameter_circle(c: &rcad_kernel::geom::Circle3, p: DVec3) -> f64 {
    let d = p - c.center;
    let u = d - d.dot(c.normal) * c.normal;
    u.dot(c.x_dir).atan2(u.dot(c.y_dir))
}

/// OCCT ElCLib::Parameter(E, P).
pub fn elclib_parameter_ellipse(e: &rcad_kernel::geom::Ellipse3, p: DVec3) -> f64 {
    let d = p - e.center;
    let u = d - d.dot(e.normal) * e.normal;
    u.dot(e.major_dir).atan2(e.normal.cross(e.major_dir).dot(u))
}

/// OCCT ElCLib::D1(U, C, P, V) — the derivative of a circle.
pub fn circle_d1(c: &rcad_kernel::geom::Circle3, u: f64) -> DVec3 {
    c.radius * (-u.sin() * c.x_dir + u.cos() * c.y_dir)
}

/// OCCT ElCLib::D1(U, E, P, V) — the derivative of an ellipse.
pub fn ellipse_d1(e: &rcad_kernel::geom::Ellipse3, u: f64) -> DVec3 {
    let ydir = e.normal.cross(e.major_dir).normalize();
    -e.major_radius * u.sin() * e.major_dir + e.minor_radius * u.cos() * ydir
}

/// OCCT Geom_Line::Reverse / Geom_Circle::Reverse / Geom_Ellipse::Reverse
/// (the parameterization is flipped in place).
pub fn reverse_curve(c: &rcad_kernel::geom::Curve3) -> rcad_kernel::geom::Curve3 {
    match c {
        rcad_kernel::geom::Curve3::Line(l) => {
            rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: l.origin,
                direction: -l.direction,
            })
        }
        rcad_kernel::geom::Curve3::Circle(ci) => {
            rcad_kernel::geom::Curve3::Circle(rcad_kernel::geom::Circle3 {
                center: ci.center,
                normal: -ci.normal,
                x_dir: ci.x_dir,
                y_dir: ci.y_dir,
                radius: ci.radius,
            })
        }
        rcad_kernel::geom::Curve3::Ellipse(e) => {
            rcad_kernel::geom::Curve3::Ellipse(rcad_kernel::geom::Ellipse3 {
                center: e.center,
                normal: -e.normal,
                major_dir: e.major_dir,
                major_radius: e.major_radius,
                minor_radius: e.minor_radius,
            })
        }
        other => other.clone(),
    }
}

#[allow(clippy::too_many_arguments)]
pub fn chfi3d_compute_curves(
    s1: &GeomAdaptorSurface,
    s2: &GeomAdaptorSurface,
    pardeb: [f64; 4],
    parfin: [f64; 4],
    tol3d: f64,
    tol2d: f64,
    tolreached: &mut f64,
) -> Option<ComputedCurves> {
    use rcad_kernel::core::precision::ANGULAR;
    use rcad_kernel::geom::{CurveEval as _, SurfaceEval as _};

    let pdeb1 = s1.value(pardeb[0], pardeb[1]);
    let pfin1 = s1.value(parfin[0], parfin[1]);
    let pdeb2 = s2.value(pardeb[2], pardeb[3]);
    let pfin2 = s2.value(parfin[2], parfin[3]);

    let mut distrefdeb = pdeb1.distance(pdeb2); // checks the worthiness
    let mut distreffin = pfin1.distance(pfin2); // of input data
    if distrefdeb < tol3d {
        distrefdeb = tol3d;
    }
    if distreffin < tol3d {
        distreffin = tol3d;
    }

    let pdeb = 0.5 * (pdeb1 + pdeb2);
    let pfin = 0.5 * (pfin1 + pfin2);

    let mut distref = 0.005 * pdeb.distance(pfin);
    if distref < distrefdeb {
        distref = distrefdeb;
    }
    if distref < distreffin {
        distref = distreffin;
    }

    // To reorientate the result of the analytic intersection, the beginning
    // of the tangent should be in the direction of the start/end line.
    let vref = pfin - pdeb;
    *tolreached = tol3d;

    let type1 = s1.get_type();
    let type2 = s2.get_type();
    if (type1 == GeomAbsSurfaceType::Cylinder && type2 == GeomAbsSurfaceType::Plane)
        || (type1 == GeomAbsSurfaceType::Plane && type2 == GeomAbsSurfaceType::Cylinder)
    {
        let (pl, cyl) = if type1 == GeomAbsSurfaceType::Plane {
            (s1.plane(), s2.cylinder())
        } else {
            (s2.plane(), s1.cylinder())
        };

        // OCCT L3760: IntAna_QuadQuadGeo ImpKK(pl, cyl, Angular, tol3d)
        // (the H height defaults to 0 — infinite cylinder).
        let quad_pl = rcad_kernel::geom::Surface3::Plane(pl.clone());
        let quad_cyl = rcad_kernel::geom::Surface3::Cylinder(cyl.clone());
        let q1 = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(&quad_pl)
            .expect("plane quadric");
        let q2 = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(&quad_cyl)
            .expect("cylinder quadric");
        let mut imp_kk = crate::geomalgo::int_patch::quad_quad_geo::QuadQuadGeo::new();
        imp_kk.perform_plane_cylinder(&q1, &q2, ANGULAR, tol3d, 0.0);
        let mut is_int_done = imp_kk.is_done();

        if imp_kk.type_inter() == crate::geomalgo::int_patch::quad_quad_geo::AnaResultType::Ellipse {
            let an_el = imp_kk.ellipse();
            let a_major_r = an_el.major_radius;
            let a_minor_r = an_el.minor_radius;
            is_int_done = a_major_r < 100000.0 * a_minor_r;
        }

        if !is_int_done {
            return None;
        }

        use crate::geomalgo::int_patch::quad_quad_geo::AnaResultType as T;
        let (mut c3d, udeb0, ufin0, vint, c1line): (
            rcad_kernel::geom::Curve3,
            f64,
            f64,
            DVec3,
            bool,
        ) = match imp_kk.type_inter() {
            T::Line => {
                let nbsol = imp_kk.nb_solutions();
                let mut c1 = imp_kk.line(1);
                let mut udeb = 0.0f64;
                for ilin in 1..=nbsol {
                    c1 = imp_kk.line(ilin);
                    udeb = elclib_parameter_line(&c1, pdeb);
                    let ptest = elclib_line_value(udeb, &c1);
                    if ptest.distance(pdeb) < tol3d {
                        break;
                    }
                }
                let ufin = elclib_parameter_line(&c1, pfin);
                let vint = c1.direction;
                (
                    rcad_kernel::geom::Curve3::Line(c1.clone()),
                    udeb,
                    ufin,
                    vint,
                    true,
                )
            }
            T::Circle => {
                let c1 = imp_kk.circle();
                let udeb = elclib_parameter_circle(&c1, pdeb);
                let ufin = elclib_parameter_circle(&c1, pfin);
                let vint = circle_d1(&c1, udeb);
                (
                    rcad_kernel::geom::Curve3::Circle(c1.clone()),
                    udeb,
                    ufin,
                    vint,
                    false,
                )
            }
            T::Ellipse => {
                let c1 = imp_kk.ellipse();
                let udeb = elclib_parameter_ellipse(&c1, pdeb);
                let ufin = elclib_parameter_ellipse(&c1, pfin);
                let vint = ellipse_d1(&c1, udeb);
                (
                    rcad_kernel::geom::Curve3::Ellipse(c1.clone()),
                    udeb,
                    ufin,
                    vint,
                    false,
                )
            }
            _ => {
                // OCCT default: C3d stays null — the caller raises.
                return None;
            }
        };

        let mut udeb = udeb0;
        let mut ufin = ufin0;
        if vint.dot(vref) < 0.0 {
            let c3dr = reverse_curve(&c3d);
            if c1line {
                udeb = -udeb;
                ufin = -ufin;
            } else {
                udeb = 2.0 * std::f64::consts::PI - udeb;
                ufin = 2.0 * std::f64::consts::PI - ufin;
            }
            c3d = c3dr;
        }
        if !c1line {
            // OCCT L3830: ElCLib::AdjustPeriodic(0, 2PI, Angular, Udeb, Ufin).
            let (au, bu) = elclib_adjust_periodic(0.0, 2.0 * std::f64::consts::PI, ANGULAR, udeb, ufin);
            udeb = au;
            ufin = bu;
        }

        // OCCT L3834-3857: ProjectPCurv on S1/S2 with the cylinder pcurve
        // translation (Translate(gp_Vec2d)) when the start point deviates.
        let (pc1, tolr1) = chfi3d_project_pcurv(&c3d, s1, tol3d)?;
        let pc1 = translate_pcurve_if_needed(pc1, s1, pardeb[0], pardeb[1], udeb, tol2d);
        let (pc2, tolr2) = chfi3d_project_pcurv(&c3d, s2, tol3d)?;
        let pc2 = translate_pcurve_if_needed(pc2, s2, pardeb[2], pardeb[3], udeb, tol2d);

        let c3dt = rcad_kernel::geom::Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            c3d, udeb, ufin,
        ));
        *tolreached = 1.5 * tolr1.max(tolr2);
        let ev = chfi3d_eval_tol_reached(s1, &pc1, s2, &pc2, &c3dt);
        *tolreached = (*tolreached).min(ev);
        return Some(ComputedCurves {
            c3d: c3dt,
            pc1,
            pc2,
            tolreached: *tolreached,
        });
    } else if type1 == GeomAbsSurfaceType::Plane && type2 == GeomAbsSurfaceType::Plane {
        // OCCT L3864-3887: IntAna_QuadQuadGeo LInt(S1->Plane(), S2->Plane(),
        // Angular, tol3d).
        let quad1 = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(&rcad_kernel::geom::Surface3::Plane(s1.plane()))
            .expect("plane quadric");
        let quad2 = crate::geomalgo::int_surf::quadric::Quadric::from_surface3(&rcad_kernel::geom::Surface3::Plane(s2.plane()))
            .expect("plane quadric");
        let mut lint = crate::geomalgo::int_patch::quad_quad_geo::QuadQuadGeo::new();
        lint.perform_plane_plane(&quad1, &quad2, ANGULAR, tol3d);
        if lint.is_done() {
            let l = lint.line(1);
            let mut c3d = rcad_kernel::geom::Curve3::Line(l.clone());
            let mut udeb = elclib_parameter_line(&l, pdeb);
            let mut ufin = elclib_parameter_line(&l, pfin);
            let vint = l.direction;
            if vint.dot(vref) < 0.0 {
                c3d = reverse_curve(&c3d);
                udeb = -udeb;
                ufin = -ufin;
            }
            let (pc1, _) = chfi3d_project_pcurv(&c3d, s1, tol3d)?;
            let (pc2, _) = chfi3d_project_pcurv(&c3d, s2, tol3d)?;
            let c3dt = rcad_kernel::geom::Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
                c3d, udeb, ufin,
            ));
            return Some(ComputedCurves {
                c3d: c3dt,
                pc1,
                pc2,
                tolreached: *tolreached,
            });
        }
        None
    } else {
        // OCCT L3888-4070: generic branch — trsfsurf + GeomInt_IntSS with
        // tolap = 2.e-7, plus the L4075-4293 IntWalk_PWalking fallback.
        // Pending TKTopAlgo/TKGeomAlgo translations: report the OCCT
        // failure path.
        None
    }
}

/// OCCT ElCLib::AdjustPeriodic(UFirst, ULast, Eps, U1, U2) — brings (U1, U2)
/// into the period with U2 > U1.
pub fn elclib_adjust_periodic(
    ufirst: f64,
    ulast: f64,
    _eps: f64,
    u1: f64,
    u2: f64,
) -> (f64, f64) {
    let period = ulast - ufirst;
    let mut a1 = u1;
    let mut a2 = u2;
    while a1 < ufirst {
        a1 += period;
    }
    while a1 >= ulast {
        a1 -= period;
    }
    while a2 <= a1 {
        a2 += period;
    }
    (a1, a2)
}

/// OCCT L3835-3845 / L3846-3857: when the projected pcurve start point
/// deviates from Pardeb by more than tol2d, translate the pcurve by the
/// difference (cylinder-parameter periodicity shift).
fn translate_pcurve_if_needed(
    pc: rcad_kernel::geom::Curve2d,
    s: &GeomAdaptorSurface,
    refu: f64,
    refv: f64,
    at: f64,
    tol2d: f64,
) -> rcad_kernel::geom::Curve2d {
    use rcad_kernel::geom::Curve2dEval as _;
    if s.get_type() == GeomAbsSurfaceType::Cylinder {
        let uv = pc.point_at(at);
        let x = refu - uv.x;
        let y = refv - uv.y;
        if x.abs() >= tol2d || y.abs() >= tol2d {
            return translate_curve2d(&pc, DVec2::new(x, y));
        }
    }
    pc
}

/// OCCT Geom2d_Curve::Translate(gp_Vec2d).
pub fn translate_curve2d(pc: &rcad_kernel::geom::Curve2d, t: DVec2) -> rcad_kernel::geom::Curve2d {
    match pc {
        rcad_kernel::geom::Curve2d::Line(l) => {
            rcad_kernel::geom::Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: l.origin + t,
                direction: l.direction,
            })
        }
        rcad_kernel::geom::Curve2d::Circle(c) => {
            rcad_kernel::geom::Curve2d::Circle(rcad_kernel::geom::Circle2d {
                center: c.center + t,
                ..c.clone()
            })
        }
        rcad_kernel::geom::Curve2d::BSpline(b) => {
            let mut b = b.clone();
            for p in b.control_points.iter_mut() {
                *p += t;
            }
            rcad_kernel::geom::Curve2d::BSpline(b)
        }
        rcad_kernel::geom::Curve2d::Bezier(b) => {
            let mut b = b.clone();
            for p in b.control_points.iter_mut() {
                *p += t;
            }
            rcad_kernel::geom::Curve2d::Bezier(b)
        }
        other => {
            // Wrap into a translated BSpline is not applicable; the
            // remaining analytic kinds carry their own frame and are not
            // produced by ProjectPCurv here.
            let _ = t;
            other.clone()
        }
    }
}
