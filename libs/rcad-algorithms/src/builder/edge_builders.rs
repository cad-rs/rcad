use glam::DVec2;
use rcad_kernel::geom::*;
use crate::bopds::ds::*;
use crate::tolerance::*;
use crate::builder::types::{WireSegment, WireEdgeSource, WireOrientation};
use super::wire_splitter::{world_to_uv, edge_uv_tangent, compute_seam_tangent_angles};

/// ✅ OCCT-aligned: degenerate edge segments (BOPAlgo_Builder_2.cxx L408-412).
pub fn build_degenerate_edge_segments(ds: &DS, ei: usize, sv: usize, ev: usize, face: &DSFace, face_idx: usize) -> Vec<WireSegment> {
    let deg_pcurve = match &face.surface {
        Surface3::Sphere(_) => {
            let pole_v = world_to_uv(&face.surface, ds.vertices[sv].point)
                .map(|uv| uv.y).unwrap_or(0.0);
            let mut ic_uvs: Vec<f64> = Vec::new();
            for &ci in &face.face_info.curves_sc {
                let ic = &ds.intersection_curves[ci];
                if let Curve3::Circle(c) = &ic.curve { if c.radius < 1e-3 { continue; } }
                let pole_pt = ds.vertices[sv].point;
                let tol_sq = TOLERANCE_ABS_SQ * 1_000_000.0;
                let at_s = ds.vertices[ic.start_vertex].point.distance_squared(pole_pt) <= tol_sq;
                let at_e = ds.vertices[ic.end_vertex].point.distance_squared(pole_pt) <= tol_sq;
                if !at_s && !at_e { continue; }
                let t = if at_s { ic.t_range[0] } else { ic.t_range[1] };
                if let Some(pc) = ic.pcurve_on_b.as_ref().or(ic.pcurve_on_a.as_ref()) {
                    let uv: DVec2 = pc.point_at(t);
                    let u = uv.x;
                    if u.abs() > 0.01 && (u - std::f64::consts::PI).abs() > 0.01
                        && (u - std::f64::consts::TAU).abs() > 0.01 { ic_uvs.push(uv.x); }
                }
            }
            let ic_u = if ic_uvs.is_empty() { std::f64::consts::TAU }
                       else { ic_uvs.iter().sum::<f64>() / ic_uvs.len() as f64 };
            Some(Curve2d::Line(Line2d { origin: DVec2::new(ic_u, pole_v), direction: DVec2::new(-ic_u, 0.0) }))
        }
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Torus(_) => {
            world_to_uv(&face.surface, ds.vertices[sv].point).map(|uv| {
                Curve2d::Line(Line2d { origin: DVec2::new(0.0, uv.y), direction: DVec2::new(std::f64::consts::TAU, 0.0) })
            })
        }
        _ => world_to_uv(&face.surface, ds.vertices[sv].point).map(|uv| {
            Curve2d::Line(Line2d { origin: DVec2::new(0.0, uv.y), direction: DVec2::new(std::f64::consts::TAU, 0.0) })
        }),
    };
    let tangent = compute_seam_tangent_angles(ds, ei, sv, ev, &face.surface);
    let fwd = WireSegment {
        start_vertex: sv, end_vertex: sv,
        source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
        is_seam: true, second_pcurve: deg_pcurve.clone(), first_pcurve: None, t_range: [0.0, 1.0],
        tangent_start: tangent.0, tangent_end: tangent.1,
    };
    let deg_pcurve_rev = match &deg_pcurve {
        Some(Curve2d::Line(l)) => {
            let ic_u = l.origin.x; let pole_v = l.origin.y;
            if (ic_u - std::f64::consts::TAU).abs() < 1e-10 {
                Some(Curve2d::Line(Line2d { origin: DVec2::new(0.0, pole_v), direction: DVec2::new(std::f64::consts::TAU, 0.0) }))
            } else {
                Some(Curve2d::Line(Line2d { origin: DVec2::new(std::f64::consts::TAU, pole_v), direction: DVec2::new(-(std::f64::consts::TAU - ic_u), 0.0) }))
            }
        }
        _ => None,
    };
    let rev = WireSegment {
        start_vertex: sv, end_vertex: sv,
        source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
        is_seam: true, second_pcurve: deg_pcurve_rev, first_pcurve: None, t_range: [0.0, 1.0],
        tangent_start: Some(std::f64::consts::PI), tangent_end: Some(std::f64::consts::PI),
    };
    vec![fwd, rev]
}

/// ✅ OCCT: DoSplitSEAMOnFace — build second pcurve shifted by surface period.
pub fn build_seam_second_pcurve(ds: &DS, surface: &Surface3, sv: usize, ev: usize, edge_tol: f64) -> Option<Curve2d> {
    let (is_periodic, period, u_min, u_max) = match surface {
        Surface3::Sphere(_) | Surface3::Cylinder(_) => (true, std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
        _ => (false, 0.0, 0.0, 0.0),
    };
    if !is_periodic { return None; }
    let mid_3d = (ds.vertices[sv].point + ds.vertices[ev].point) * 0.5;
    let uv_mid = world_to_uv(surface, mid_3d)?;
    let dU = match surface {
        Surface3::Sphere(sph) => edge_tol / sph.radius.max(1e-15),
        Surface3::Cylinder(cyl) => edge_tol / cyl.radius.max(1e-15),
        _ => TOLERANCE_ABS,
    };
    let shift_u = if (uv_mid.x - u_min).abs() < dU { period } else if (uv_mid.x - u_max).abs() < dU { -period } else { return None; };
    let uv_s = world_to_uv(surface, ds.vertices[sv].point)?;
    let uv_e = world_to_uv(surface, ds.vertices[ev].point)?;
    Some(Curve2d::Line(Line2d { origin: DVec2::new(uv_s.x + shift_u, uv_s.y), direction: DVec2::new(0.0, uv_e.y - uv_s.y) }))
}

/// ✅ OCCT-aligned: split seam sub-edges for sphere (DoSplitSEAMOnFace).
pub fn build_sphere_seam_segments(ds: &DS, ei: usize, sv: usize, ev: usize, face: &DSFace, _face_idx: usize) -> Vec<WireSegment> {
    let ds_edge = &ds.edges[ei];
    let mut segs: Vec<WireSegment> = Vec::new();
    if ds_edge.pave_blocks.len() > 1 {
        for pb in &ds_edge.pave_blocks {
            let sv_seg = pb.pave1.vertex_idx; let ev_seg = pb.pave2.vertex_idx;
            if sv_seg == ev_seg { continue; }
            let (t_start, t_end) = compute_seam_tangent_angles(ds, ei, sv_seg, ev_seg, &face.surface);
            let second_pcurve = build_seam_second_pcurve(ds, &face.surface, sv_seg, ev_seg, ds_edge.geom_tol);
            let first_pcurve = world_to_uv(&face.surface, ds.vertices[sv_seg].point).and_then(|uv_s| {
                world_to_uv(&face.surface, ds.vertices[ev_seg].point).map(|uv_e| {
                    Curve2d::Line(Line2d { origin: DVec2::new(uv_s.x, uv_s.y), direction: DVec2::new(0.0, uv_e.y - uv_s.y) })
                })
            });
            segs.push(WireSegment {
                start_vertex: sv_seg, end_vertex: ev_seg,
                source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
                is_seam: true, second_pcurve: second_pcurve.clone(), first_pcurve, t_range: [0.0, 1.0],
                tangent_start: t_start, tangent_end: t_end,
            });
            let (t_start_rev, t_end_rev) = compute_seam_tangent_angles(ds, ei, ev_seg, sv_seg, &face.surface);
            let second_pcurve_rev = second_pcurve.map(|pc| match pc {
                Curve2d::Line(l) => Curve2d::Line(Line2d { origin: l.origin + l.direction, direction: -l.direction }),
                _ => pc,
            });
            segs.push(WireSegment {
                start_vertex: ev_seg, end_vertex: sv_seg,
                source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
                is_seam: true, second_pcurve: second_pcurve_rev, first_pcurve: None, t_range: [0.0, 1.0],
                tangent_start: t_start_rev, tangent_end: t_end_rev,
            });
        }
    } else {
        let second_pcurve = build_seam_second_pcurve(ds, &face.surface, sv, ev, ds_edge.geom_tol);
        let first_pcurve = world_to_uv(&face.surface, ds.vertices[sv].point).and_then(|uv_sv| {
            world_to_uv(&face.surface, ds.vertices[ev].point).map(|uv_ev| {
                Curve2d::Line(Line2d { origin: DVec2::new(uv_sv.x, uv_sv.y), direction: DVec2::new(0.0, uv_ev.y - uv_sv.y) })
            })
        });
        let (ts, te) = compute_seam_tangent_angles(ds, ei, sv, ev, &face.surface);
        segs.push(WireSegment {
            start_vertex: sv, end_vertex: ev,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
            is_seam: true, second_pcurve: second_pcurve, first_pcurve,
            t_range: [0.0, 1.0], tangent_start: ts, tangent_end: te,
        });
    }
    segs
}

/// ✅ OCCT-aligned: cylinder/cone seam edge — FWD+REV with shifted pcurves.
pub fn build_cylinder_seam_segments(ds: &DS, ei: usize, sv: usize, ev: usize, face: &DSFace) -> Vec<WireSegment> {
    let (t_start, t_end) = compute_seam_tangent_angles(ds, ei, sv, ev, &face.surface);
    let uv_a = world_to_uv(&face.surface, ds.vertices[sv].point);
    let uv_b = world_to_uv(&face.surface, ds.vertices[ev].point);
    let (pcurve_opt, second_pcurve_opt) = match (uv_a, uv_b) {
        (Some(ua), Some(ub)) => {
            let p0 = DVec2::new(ua.x, ua.y); let p1 = DVec2::new(ub.x, ub.y);
            let dir = p1 - p0;
            let first = Curve2d::Line(Line2d { origin: p0, direction: dir });
            let is_periodic = matches!(face.surface, Surface3::Cylinder(_) | Surface3::Sphere(_));
            let second = if is_periodic { Curve2d::Line(Line2d { origin: p0 + DVec2::new(std::f64::consts::TAU, 0.0), direction: dir }) } else { first.clone() };
            (Some(first), Some(second))
        }
        _ => (None, None),
    };
    vec![
        WireSegment {
            start_vertex: sv, end_vertex: ev,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
            is_seam: true, second_pcurve: None, first_pcurve: pcurve_opt, t_range: [0.0, 1.0],
            tangent_start: t_start, tangent_end: t_end,
        },
        WireSegment {
            start_vertex: ev, end_vertex: sv,
            source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
            is_seam: true, second_pcurve: second_pcurve_opt, first_pcurve: None, t_range: [0.0, 1.0],
            tangent_start: t_end.map(|a| (a + std::f64::consts::PI) % std::f64::consts::TAU),
            tangent_end: t_start.map(|a| (a + std::f64::consts::PI) % std::f64::consts::TAU),
        },
    ]
}
