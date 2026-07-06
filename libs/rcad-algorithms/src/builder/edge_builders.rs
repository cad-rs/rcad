use glam::DVec2;
use rcad_kernel::geom::*;
use crate::bopds::ds::*;
use crate::tolerance::*;
use crate::builder::types::{WireSegment, WireEdgeSource, WireOrientation};
use super::wire_splitter::{world_to_uv, edge_uv_tangent};
use crate::builder::curve_eq;

///  ?OCCT: DoSplitSEAMOnFace  ?build second pcurve shifted by surface period.
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

///  ?OCCT-aligned: split seam sub-edges for sphere (DoSplitSEAMOnFace).
pub fn build_sphere_seam_segments(ds: &DS, ei: usize, sv: usize, ev: usize, face: &DSFace, _face_idx: usize) -> Vec<WireSegment> {
 let ds_edge = &ds.edges[ei];
 let mut segs: Vec<WireSegment> = Vec::new();
 if ds_edge.pave_blocks.len() > 1 {
 for pb in &ds_edge.pave_blocks {
 let sv_seg = pb.pave1.vertex_idx; let ev_seg = pb.pave2.vertex_idx;
 if sv_seg == ev_seg { continue; }
 let second_pcurve = build_seam_second_pcurve(ds, &face.surface, sv_seg, ev_seg, ds_edge.geom_tol);
 let first_pcurve = world_to_uv(&face.surface, ds.vertices[sv_seg].point).and_then(|uv_s| {
 world_to_uv(&face.surface, ds.vertices[ev_seg].point).map(|uv_e| {
 Curve2d::Line(Line2d { origin: DVec2::new(uv_s.x, uv_s.y), direction: DVec2::new(0.0, uv_e.y - uv_s.y) })
 })
 });
 segs.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
 is_closed_on_face: true, second_pcurve: second_pcurve.clone(), first_pcurve, t_range: [0.0, 1.0],
 });
 let second_pcurve_rev = second_pcurve.map(|pc| match pc {
 Curve2d::Line(l) => Curve2d::Line(Line2d { origin: l.origin + l.direction, direction: -l.direction }),
 _ => pc,
 });
 segs.push(WireSegment {
 start_vertex: ev_seg, end_vertex: sv_seg,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
 is_closed_on_face: true, second_pcurve: second_pcurve_rev, first_pcurve: None, t_range: [0.0, 1.0],
 });
 }
 } else {
 let second_pcurve = build_seam_second_pcurve(ds, &face.surface, sv, ev, ds_edge.geom_tol);
 let first_pcurve = world_to_uv(&face.surface, ds.vertices[sv].point).and_then(|uv_sv| {
 world_to_uv(&face.surface, ds.vertices[ev].point).map(|uv_ev| {
 Curve2d::Line(Line2d { origin: DVec2::new(uv_sv.x, uv_sv.y), direction: DVec2::new(0.0, uv_ev.y - uv_sv.y) })
 })
 });
 let (ts, te): (Option<f64>, Option<f64>) = (None, None);
 segs.push(WireSegment {
 start_vertex: sv, end_vertex: ev,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
 is_closed_on_face: true, second_pcurve: None, first_pcurve,
 t_range: [0.0, 1.0],
 });
 let (ts_rev, te_rev): (Option<f64>, Option<f64>) = (None, None);
 segs.push(WireSegment {
 start_vertex: ev, end_vertex: sv,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
 is_closed_on_face: true, second_pcurve, first_pcurve: None, t_range: [0.0, 1.0],
 });
 }
 segs
}

///  ?OCCT-aligned: cylinder/cone seam edge  ?FWD+REV with shifted pcurves.
pub fn build_cylinder_seam_segments(ds: &DS, ei: usize, sv: usize, ev: usize, face: &DSFace) -> Vec<WireSegment> {
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
 is_closed_on_face: true, second_pcurve: None, first_pcurve: pcurve_opt, t_range: [0.0, 1.0],
 },
 WireSegment {
 start_vertex: ev, end_vertex: sv,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
 is_closed_on_face: true, second_pcurve: second_pcurve_opt, first_pcurve: None, t_range: [0.0, 1.0],
 },
 ]
}

///  ?OCCT-aligned: IsSplitToReverseWithWarn (BOPTools_AlgoTools.cxx L1432-1523).
/// Compares the direction of a split sub-edge against its original edge.
pub fn is_split_to_reverse(
 ds: &DS,
 sub_ei: usize,
 orig_ei: usize,
) -> bool {
 let sub_edge = &ds.edges[sub_ei];
 let orig_edge = &ds.edges[orig_ei];
 // OCCT L1441-1448: skip degenerated edges
 if ds.is_edge_degenerated(sub_ei) || ds.is_edge_degenerated(orig_ei) {
 return false;
 }
 // OCCT L1462-1465: same 3D curve -> compare parameter direction
 if curve_eq(&sub_edge.curve, &orig_edge.curve) {
 let sub_dir = sub_edge.t_range[1] - sub_edge.t_range[0];
 let orig_dir = orig_edge.t_range[1] - orig_edge.t_range[0];
 return (sub_dir > 0.0) != (orig_dir > 0.0);
 }
 // OCCT L1479-1514: sample midpoint, project onto original, compare tangents
 let sub_mid = (sub_edge.t_range[0] + sub_edge.t_range[1]) * 0.5;
 let sub_mid_pt = sub_edge.curve.point_at(sub_mid);
 let n = 30;
 let mut best_t = orig_edge.t_range[0];
 let mut best_d = f64::MAX;
 for i in 0..=n {
 let t = orig_edge.t_range[0]
 + (orig_edge.t_range[1] - orig_edge.t_range[0]) * (i as f64 / n as f64);
 let p = orig_edge.curve.point_at(t);
 let d = sub_mid_pt.distance_squared(p);
 if d < best_d {
 best_d = d;
 best_t = t;
 }
 }
 if best_d > 1e-6 {
 return false;
 }
 let eps = 1e-8;
 let sub_t1 = sub_edge.curve.point_at((sub_mid + eps).min(sub_edge.t_range[1]));
 let sub_t2 = sub_edge.curve.point_at((sub_mid - eps).max(sub_edge.t_range[0]));
 let tangent_sub = sub_t1 - sub_t2;
 let orig_t1 = orig_edge.curve.point_at((best_t + eps).min(orig_edge.t_range[1]));
 let orig_t2 = orig_edge.curve.point_at((best_t - eps).max(orig_edge.t_range[0]));
 let tangent_orig = orig_t1 - orig_t2;
 tangent_sub.dot(tangent_orig) < 0.0
}
