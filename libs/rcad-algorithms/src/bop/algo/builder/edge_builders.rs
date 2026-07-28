use super::wire_splitter::{edge_uv_tangent, world_to_uv};
use crate::bop::algo::builder::curve_eq;
use crate::bop::algo::builder::types::{WireEdgeSource, WireOrientation, WireSegment};
use crate::bop::ds::ds::*;
use crate::tolerance::*;
use glam::DVec2;
use rcad_kernel::geom::*;

///  ?OCCT: DoSplitSEAMOnFace  ?build second pcurve shifted by surface period.
pub fn build_seam_second_pcurve(
    ds: &DS,
    surface: &Surface3,
    sv: usize,
    ev: usize,
    edge_tol: f64,
) -> Option<Curve2d> {
    let (is_periodic, period, u_min, u_max) = match surface {
        Surface3::Sphere(_) | Surface3::Cylinder(_) => {
            (true, std::f64::consts::TAU, 0.0, std::f64::consts::TAU)
        }
        _ => (false, 0.0, 0.0, 0.0),
    };
    if !is_periodic {
        return None;
    }
    let mid_3d = (ds.vertex_point(sv) + ds.vertex_point(ev)) * 0.5;
    let uv_mid = world_to_uv(surface, mid_3d)?;
    let dU = match surface {
        Surface3::Sphere(sph) => edge_tol / sph.radius.max(TOLERANCE_CLAMP_MIN),
        Surface3::Cylinder(cyl) => edge_tol / cyl.radius.max(TOLERANCE_CLAMP_MIN),
        _ => TOLERANCE_ABS,
    };
    let shift_u = if (uv_mid.x - u_min).abs() < dU {
        period
    } else if (uv_mid.x - u_max).abs() < dU {
        -period
    } else {
        return None;
    };
    let uv_s = world_to_uv(surface, ds.vertex_point(sv))?;
    let uv_e = world_to_uv(surface, ds.vertex_point(ev))?;
    Some(Curve2d::Line(Line2d {
        origin: DVec2::new(uv_s.x + shift_u, uv_s.y),
        direction: DVec2::new(0.0, uv_e.y - uv_s.y),
    }))
}


///  IsSplitToReverseWithWarn (BOPTools_AlgoTools.cxx L1432-1523).
/// Compares the direction of a split sub-edge against its original edge.
pub fn is_split_to_reverse(ds: &DS, sub_ei: usize, orig_ei: usize) -> bool {
    let sub_curve = ds.edge_curve(sub_ei);
    let orig_curve = ds.edge_curve(orig_ei);
    let (Some(sub_curve), Some(orig_curve)) = (sub_curve, orig_curve) else {
        return false;
    };
    // Skip degenerated edges
    if ds.is_edge_degenerated(sub_ei) || ds.is_edge_degenerated(orig_ei) {
        return false;
    }
    let sub_range = ds.edge_range(sub_ei);
    let orig_range = ds.edge_range(orig_ei);
    // Same 3D curve  ?compare parameter direction
    if curve_eq(sub_curve, orig_curve) {
        let sub_dir = sub_range[1] - sub_range[0];
        let orig_dir = orig_range[1] - orig_range[0];
        return (sub_dir > 0.0) != (orig_dir > 0.0);
    }
    // Sample midpoint, project onto original, compare tangents
    let sub_mid = (sub_range[0] + sub_range[1]) * 0.5;
    let sub_mid_pt = sub_curve.point_at(sub_mid);
    let n = 30;
    let mut best_t = orig_range[0];
    let mut best_d = f64::MAX;
    for i in 0..=n {
        let t = orig_range[0]
            + (orig_range[1] - orig_range[0]) * (i as f64 / n as f64);
        let p = orig_curve.point_at(t);
        let d = sub_mid_pt.distance_squared(p);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    if best_d > TOLERANCE_MESH_LEGACY {
        return false;
    }
    let eps = 1e-8;
    let sub_t1 = sub_curve
        .point_at((sub_mid + eps).min(sub_range[1]));
    let sub_t2 = sub_curve
        .point_at((sub_mid - eps).max(sub_range[0]));
    let tangent_sub = sub_t1 - sub_t2;
    let orig_t1 = orig_curve
        .point_at((best_t + eps).min(orig_range[1]));
    let orig_t2 = orig_curve
        .point_at((best_t - eps).max(orig_range[0]));
    let tangent_orig = orig_t1 - orig_t2;
    tangent_sub.dot(tangent_orig) < 0.0
}
