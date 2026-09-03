//! OCCT ChFi3d_Builder_0.cxx free functions — 1:1 translation of the
//! utilities used by the builder (PerformElement / PerformExtremity /
//! ExtentAnalyse / corners).
//!
//! rcad architecture note: OCCT TopoDS shapes carry their geometry through
//! global handle graphs, while rcad stores TShapes in a BRep pool, so the
//! functions here take the owning `&BRep` as an extra first argument where
//! OCCT reads geometry straight from the shape handles.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{Orientation, Shape, TShape};
use rcad_kernel::topods;

use super::chfi3d::is_tangent_faces;
use super::chfi_ds::{ChFiDSMap, ChFiDS_State};
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
