use crate::tolerance::*;
use std::sync::Arc;

use glam::DVec3;
use rcad_kernel::PCurve;
use rcad_kernel::geom::*;
use rcad_kernel::topods::{self, TShape};

use crate::inttools::pcurve_derive::fallback_pcurve_by_projection;

/// Populates `brep.geom` with analytic geometry for a box BRep.
///
/// After this call, every edge has a `Curve3::Line` and every face has a `Surface3::Plane`.
/// Precondition: brep was created by `BRep::from_primitive(Box{..})`.
pub fn populate_box_geom(brep: &mut rcad_kernel::BRep) {
    // Collect edge indices (ts_index, first_idx, last_idx) for edges missing a curve.
    let edge_data: Vec<(usize, usize, usize)> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(ti, ts)| {
            if let TShape::Edge(ed) = ts.as_ref() {
                if ed.curve.is_none() {
                    Some((ti, ed.first.index, ed.last.index))
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect();

    // Set line curves in a separate pass (avoids borrow conflict with vertex_point).
    for (ei, first_idx, last_idx) in &edge_data {
        let p0 = brep.vertex_point(*first_idx).unwrap_or(DVec3::ZERO);
        let p1 = brep.vertex_point(*last_idx).unwrap_or(DVec3::ZERO);
        let delta = p1 - p0;
        let len = delta.length();
        let dir = if len > TOLERANCE_LEN_MIN {
            delta / len
        } else {
            DVec3::X
        };
        let ts = &mut brep.tshapes[*ei];
        let TShape::Edge(ed) = Arc::make_mut(ts) else {
            continue;
        };
        ed.curve = Some(Curve3::Line(Line3 {
            origin: p0,
            direction: dir,
        }));
        ed.range = [0.0, len.max(TOLERANCE_LEN_MIN)];
        ed.degenerated = len <= TOLERANCE_LEN_MIN;
    }

    // Collect face indices whose surfaces need to be set.
    let face_indices: Vec<usize> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(fi, ts)| {
            if let TShape::Face(fd) = ts.as_ref() {
                if fd.surface.is_none() { Some(fi) } else { None }
            } else {
                None
            }
        })
        .collect();

    // Compute plane from wire for each face in a separate pass.
    for &fi in &face_indices {
        let (origin, normal) = {
            let ts = &brep.tshapes[fi];
            let TShape::Face(fd) = ts.as_ref() else {
                continue;
            };
            compute_face_plane_from_wire(brep, fd)
        };
        let ts = &mut brep.tshapes[fi];
        let TShape::Face(fd) = Arc::make_mut(ts) else {
            continue;
        };
        fd.surface = Some(Surface3::Plane(Plane::new(origin, normal)));
    }
}

/// Compute a plane (origin, normal) from the first three distinct vertices
/// of a face's outer wire.
fn compute_face_plane_from_wire(
    brep: &rcad_kernel::BRep,
    fd: &topods::TFaceData,
) -> (DVec3, DVec3) {
    let get_wire_pts = |sr: &topods::Shape| -> Vec<DVec3> {
        let mut pts = Vec::new();
        let Some(wire_ts) = brep.tshapes.get(sr.index) else {
            return pts;
        };
        let TShape::Wire(wd) = wire_ts.as_ref() else {
            return pts;
        };
        for er in &wd.edges {
            let Some(edge_ts) = brep.tshapes.get(er.index) else {
                continue;
            };
            let TShape::Edge(ed) = edge_ts.as_ref() else {
                continue;
            };
            let p = brep.vertex_point(ed.first.index);
            if let Some(pt) = p {
                if !pts
                    .iter()
                    .any(|q: &DVec3| (*q - pt).length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE)
                {
                    pts.push(pt);
                }
            }
        }
        pts
    };

    let pts = get_wire_pts(&fd.outer_wire);
    if pts.len() >= 3 {
        let e1 = (pts[1] - pts[0]).normalize_or_zero();
        let e2 = (pts[2] - pts[0]).normalize_or_zero();
        let n = e1.cross(e2).normalize_or_zero();
        if n.length_squared() > 0.5 {
            return (pts[0], n);
        }
    }
    // Fallback
    (pts.first().copied().unwrap_or(DVec3::ZERO), DVec3::Z)
}

/// Fixes stale plane surface origins in `geom` after boolean operations.
///
/// The boolean builder already stores the correct `Surface3` for each face
/// (including cylinders, spheres, etc.) but plane origins may point to the
/// original operand's vertex positions rather than the result's vertices.
/// This function walks every face whose surface is a `Surface3::Plane` and
/// recomputes the origin from the first vertex of the face's outer wire.
///
/// Non-plane surfaces (cylinders, spheres, cones, tori) are left untouched.
/// Existing edge curves are preserved; only missing straight-line curves are added.
pub fn recompute_plane_surfaces(brep: &mut rcad_kernel::BRep) {
    // Collect face indices and their recomputed origins.
    let mut face_origins: Vec<(usize, DVec3)> = Vec::new();
    for (fi, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Face(fd) = ts.as_ref() else {
            continue;
        };
        // Compute origin from the first vertex of the outer wire.
        let origin = brep
            .tshapes
            .get(fd.outer_wire.index)
            .and_then(|wire_ts| {
                let TShape::Wire(wd) = wire_ts.as_ref() else {
                    return None;
                };
                let first_edge_ref = wd.edges.first()?;
                let edge_ts = brep.tshapes.get(first_edge_ref.index)?;
                let TShape::Edge(ed) = edge_ts.as_ref() else {
                    return None;
                };
                brep.vertex_point(ed.first.index)
            })
            .unwrap_or(DVec3::ZERO);
        face_origins.push((fi, origin));
    }

    // Update only Plane surface origins; leave all other surface types untouched.
    for (fi, origin) in &face_origins {
        let ts = &mut brep.tshapes[*fi];
        let TShape::Face(fd) = Arc::make_mut(ts) else {
            continue;
        };
        let Some(surf) = &mut fd.surface else {
            continue;
        };
        if let Surface3::Plane(p) = surf {
            p.origin = *origin;
        }
    }

    // Add straight-line edge curves only where missing.
    // Collect (ts_index, first_idx, last_idx) for edges missing a curve.
    let edge_data: Vec<(usize, usize, usize)> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(ti, ts)| {
            if let TShape::Edge(ed) = ts.as_ref() {
                if ed.curve.is_none() {
                    Some((ti, ed.first.index, ed.last.index))
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect();

    // Set curves in a separate pass (avoids borrow conflict).
    for (ei, first_idx, last_idx) in &edge_data {
        let p0 = brep.vertex_point(*first_idx).unwrap_or(DVec3::ZERO);
        let p1 = brep.vertex_point(*last_idx).unwrap_or(DVec3::ZERO);
        let delta = p1 - p0;
        let len = delta.length();
        let dir = if len > TOLERANCE_LEN_MIN {
            delta / len
        } else {
            DVec3::X
        };
        let ts = &mut brep.tshapes[*ei];
        let TShape::Edge(ed) = Arc::make_mut(ts) else {
            continue;
        };
        ed.curve = Some(Curve3::Line(Line3 {
            origin: p0,
            direction: dir,
        }));
        ed.range = [0.0, (p1 - p0).dot(dir)];
        ed.degenerated = len <= TOLERANCE_LEN_MIN;
    }
}

/// Populate `edge_pcurves` for edges adjacent to curved faces that currently
/// lack a PCurve entry.
///
/// After a boolean operation, intersection edges on curved surfaces (cylinder,
/// sphere, cone, torus) often have no PCurve.  This function uses
/// [`fallback_pcurve_by_projection`] to derive a 2D parameter-space curve on
/// each adjacent curved surface and stores it in `brep.geom.edge_pcurves`.
///
/// Call this after [`boolean_op`] when downstream code needs PCurves
/// (e.g. parametric queries, STEP export of trimmed surfaces).
pub fn populate_boolean_result_pcurves(brep: &mut rcad_kernel::BRep) {
    // Collect all (edge_idx, face_idx) pairs where the face has a curved surface.
    // Face index = tshapes index of the TShape::Face.
    let face_surfaces: Vec<(usize, Surface3)> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(fi, ts)| {
            let TShape::Face(fd) = ts.as_ref() else {
                return None;
            };
            fd.surface.clone().map(|s| (fi, s))
        })
        .collect();

    let face_edge_pairs: Vec<(usize, usize, Surface3)> = face_surfaces
        .iter()
        .flat_map(|(fi, surf)| {
            let ts = &brep.tshapes[*fi];
            let TShape::Face(fd) = ts.as_ref() else {
                return Vec::new();
            };
            // Get edge refs from outer wire and inner wires
            let mut edge_refs = Vec::new();
            if let Some(wts) = brep.tshapes.get(fd.outer_wire.index) {
                if let TShape::Wire(wd) = wts.as_ref() {
                    for er in &wd.edges {
                        edge_refs.push(er.index);
                    }
                }
            }
            for iw in &fd.inner_wires {
                if let Some(wts) = brep.tshapes.get(iw.index) {
                    if let TShape::Wire(wd) = wts.as_ref() {
                        for er in &wd.edges {
                            edge_refs.push(er.index);
                        }
                    }
                }
            }
            edge_refs
                .into_iter()
                .map(move |ei| (ei, *fi, surf.clone()))
                .collect::<Vec<_>>()
        })
        .collect();

    for (edge_idx, face_idx, surface) in face_edge_pairs {
        // Only fill for curved surfaces.
        if matches!(surface, Surface3::Plane(_)) {
            continue;
        }

        // Check if a PCurve for this face already exists on this edge.
        let already_has = {
            let Some(ts) = brep.tshapes.get(edge_idx) else {
                continue;
            };
            let TShape::Edge(ed) = ts.as_ref() else {
                continue;
            };
            ed.pcurves.contains_key(&face_idx)
        };
        if already_has {
            continue;
        }

        // Get the curve and range from the edge.
        let (curve_opt, t_range_opt): (Option<Curve3>, Option<[f64; 2]>) = {
            let Some(ts) = brep.tshapes.get(edge_idx) else {
                continue;
            };
            let TShape::Edge(ed) = ts.as_ref() else {
                continue;
            };
            (ed.curve.clone(), Some(ed.range))
        };

        // Derive the 2D PCurve.
        let pcurve2d = if let (Some(curve), Some(t_range)) = (curve_opt, t_range_opt) {
            fallback_pcurve_by_projection(&curve, &t_range, &surface)
        } else {
            // No analytic curve: sample from vertex endpoints as a straight line.
            let Some(ts) = brep.tshapes.get(edge_idx) else {
                continue;
            };
            let TShape::Edge(ed) = ts.as_ref() else {
                continue;
            };
            let Some(p0) = brep.vertex_point(ed.first.index) else {
                continue;
            };
            let Some(p1) = brep.vertex_point(ed.last.index) else {
                continue;
            };
            if (p1 - p0).length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
                continue; // degenerate
            }
            // Project a polyline of 17 equally-spaced points between the endpoints.
            let polyline: Vec<_> = (0..17)
                .map(|i| p0 + (p1 - p0) * (i as f64 / 16.0))
                .collect();
            match crate::inttools::pcurve_derive::polyline_pcurve_by_projection(&polyline, &surface)
            {
                Some(c2d) => c2d,
                None => continue,
            }
        };

        // Store in the edge's pcurves map, keyed by face index.
        let ts = &mut brep.tshapes[edge_idx];
        let TShape::Edge(ed) = Arc::make_mut(ts) else {
            continue;
        };
        // Use a reasonable default parameter range for the pcurve.
        let prange = t_range_opt.unwrap_or([0.0, 1.0]);
        ed.pcurves
            .insert(face_idx, (pcurve2d, prange[0], prange[1]));
    }
}
