//! OCCT `BRepPrimAPI_MakeRevol` over an ARC-BEARING planar meridian.
//!
//! The analytic counterpart of `algo_ext::features::revolve_polygon_solid`,
//! which only accepts a polygon and therefore has to facet DRAW's
//! `profile … C radius angle` arcs.  DRAW's `profile` builds real
//! `Geom_Circle` arcs, and `revol` sweeps them into `Geom_SurfaceOfRevolution`
//! faces (BRepSweep_Rotation::MakeEmptyFace); a faceted approximation gives a
//! different surface mix, so OCCT's boolean grids cannot match it.
//!
//! This entry reproduces the DRAW flow: the meridian is built as a planar FACE
//! (`profile` closes the wire and builds a face by default), and the face is
//! handed to the translated `BRepSweep_Revol` (BRepSweep_Revol.cxx L29-50,
//! BRepSweep_Rotation.cxx) exactly as `BRepPrimAPI_MakeRevol` does.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation};

use super::revol::BRepSweepRevol;
use crate::algo_ext::extrude_profile::ProfileSegment;

/// OCCT `BRepPrimAPI_MakeRevol(S, A, D)` for a closed planar meridian whose
/// segments are lines and circle arcs.
///
/// `segments` is the closed meridian in order (`p1` of one segment is `p0` of
/// the next; the last closes back to the first).  A full turn yields a SOLID,
/// a partial turn a SHELL plus its two caps.
pub fn revolve_profile_solid(
    segments: &[ProfileSegment],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
) -> Result<topods::BRep, String> {
    if segments.len() < 2 {
        return Err("meridian needs >= 2 segments".into());
    }
    let dir = axis_dir.normalize_or_zero();
    if dir.length_squared() < 1e-24 {
        return Err("axis_dir must be non-zero".into());
    }
    if !(angle_rad > 0.0 && angle_rad <= std::f64::consts::TAU + 1e-9) {
        return Err("angle_rad must be in (0, 2*PI]".into());
    }
    let full_turn = (angle_rad - std::f64::consts::TAU).abs() < 1e-9;

    let mut brep = topods::BRep::new();
    // One vertex per segment start; the meridian closes back to the first.
    let verts: Vec<topods::Shape> = segments
        .iter()
        .map(|s| brep.add_tvertex(s.p0()))
        .collect();
    let rev = |s: &topods::Shape| topods::Shape {
        orientation: Orientation::Reversed,
        ..s.clone()
    };
    let mut edges: Vec<topods::Shape> = Vec::with_capacity(segments.len());
    for (i, seg) in segments.iter().enumerate() {
        let j = (i + 1) % segments.len();
        let (curve, range) = match seg {
            ProfileSegment::Line { p0, p1 } => {
                let d = *p1 - *p0;
                let len = d.length();
                if len <= 1e-12 {
                    return Err("zero-length meridian segment".into());
                }
                (Curve3::Line(Line3::new(*p0, d / len)), [0.0, len])
            }
            ProfileSegment::Arc {
                center,
                normal,
                x_dir,
                y_dir,
                radius,
                t0,
                t1,
                ..
            } => (
                Curve3::Circle(Circle3 {
                    center: *center,
                    normal: *normal,
                    x_dir: *x_dir,
                    y_dir: *y_dir,
                    radius: *radius,
                }),
                [*t0, *t1],
            ),
        };
        edges.push(brep.add_tedge(
            Some(curve),
            verts[i].clone(),
            rev(&verts[j]),
            range,
        ));
    }

    // OCCT `BRepPrimAPI_MakeRevol` sweeps the meridian FACE (DRAW `profile`
    // closes its wire and builds a face by default), so the generator handed to
    // BRepSweep_Revol is a Face and the lateral faces come from its edges.
    //
    // For a Face generator `BRepSweep_Rotation::SplitShell`
    // (BRepSweep_Rotation.cxx L797-802) runs `BRepTools_Quilt Q; Q.Add(S);
    // Q.Shells();` — BRepTools_Quilt is translated in
    // crate::topalgo::brep_tools_quilt and wired into SplitShell, so the
    // genuine OCCT flow runs here.
    let plane = meridian_plane(segments);
    let wire = brep.add_twire(edges);
    let face = brep.add_tface(Some(plane), wire, vec![], None, None, vec![], false);

    let axis = (axis_origin, dir);
    let mut the_revol = if full_turn {
        BRepSweepRevol::with_full_circle(&face, axis, false)
    } else {
        BRepSweepRevol::with_angle(&face, axis, angle_rad, false)
    };
    let out = the_revol.shape();
    if out.shape_type() == rcad_kernel::topods::ShapeType::Vertex {
        // OCCT `BRepSweep_Revol::Shape()` returns the null dummy when the
        // generator is unused by the sweep.
        return Err("BRepSweep_Revol produced no shape".into());
    }
    Ok(crate::algo_ext::topods_ext::extract_result_brep(
        &out,
        brep.locations.clone(),
    ))
}

/// The meridian's plane — an arc's own frame when present (they all lie in one
/// plane for a DRAW `profile`), else from the first non-degenerate triple.
fn meridian_plane(segments: &[ProfileSegment]) -> Surface3 {
    for seg in segments {
        if let ProfileSegment::Arc { center, normal, .. } = seg {
            return Surface3::Plane(Plane::new(*center, *normal));
        }
    }
    let p0 = segments[0].p0();
    let n = segments.len();
    for i in 1..n {
        for j in (i + 1)..n {
            let nn = (segments[i].p0() - p0).cross(segments[j].p0() - p0);
            if nn.length_squared() > 1e-18 {
                return Surface3::Plane(Plane::new(p0, nn.normalize()));
            }
        }
    }
    Surface3::Plane(Plane::new(p0, DVec3::Z))
}
