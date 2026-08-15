//! Extrude a closed planar profile made of line and circular-arc segments
//! into a solid (topods).
//!
//! Mirrors OCCT `BRepPrimAPI_MakePrism(face, vec)` for profiles with circular
//! arcs (DRAW `profile 鈥?C <r> <angle> 鈥 + `prism`): line segments sweep into
//! planar lateral faces, arc segments sweep into cylindrical lateral faces,
//! and the caps are planar faces whose wires mix line and circular edges. The
//! faceted `extrude_polygon_solid` cannot reproduce the analytic topology of
//! the OCCT reference (e.g. `bfuse_simple` I/J series), so boolean tests with
//! arc profiles use this construction.

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, CylindricalSurface, Line3, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation};

use super::tolerance::TOLERANCE_LEN_MIN;

/// One segment of a closed planar profile.
#[derive(Debug, Clone)]
pub enum ProfileSegment {
    /// Straight segment from `p0` to `p1`.
    Line { p0: DVec3, p1: DVec3 },
    /// Circular arc from `p0` to `p1`: `point(t) = center + r*(cos t*x_dir +
    /// sin t*y_dir)` for `t` in `[t0, t1]`; `p0 = point(t0)`, `p1 = point(t1)`.
    Arc {
        p0: DVec3,
        p1: DVec3,
        center: DVec3,
        normal: DVec3,
        x_dir: DVec3,
        y_dir: DVec3,
        radius: f64,
        t0: f64,
        t1: f64,
    },
}

impl ProfileSegment {
    pub fn p0(&self) -> DVec3 {
        match self {
            ProfileSegment::Line { p0, .. } => *p0,
            ProfileSegment::Arc { p0, .. } => *p0,
        }
    }
    pub fn p1(&self) -> DVec3 {
        match self {
            ProfileSegment::Line { p1, .. } => *p1,
            ProfileSegment::Arc { p1, .. } => *p1,
        }
    }
}

/// Extrude a closed planar profile along `direction` by `depth`.
///
/// The profile segments must form a closed loop (`seg[i].p1 == seg[i+1].p0`).
/// The profile plane normal is derived from the segment directions; every arc
/// must lie in that plane (the DRAW `profile` guarantee).
pub fn extrude_profile_solid(
    profile: &[ProfileSegment],
    direction: DVec3,
    depth: f64,
) -> Result<topods::BRep, super::features::FeatureError> {
    let n = profile.len();
    if n < 3 {
        return Err(super::features::FeatureError::InvalidInput("profile needs >= 3 segments"));
    }
    let dir = direction.normalize_or_zero();
    if dir.length_squared() < 1e-24 {
        return Err(super::features::FeatureError::ZeroVector("direction"));
    }
    if !depth.is_finite() || depth <= 0.0 {
        return Err(super::features::FeatureError::NonPositiveInput("depth"));
    }
    // Profile plane normal from the first two non-collinear segment chords.
    let mut cap_n = DVec3::ZERO;
    for i in 0..n {
        let d0 = profile[i].p1() - profile[i].p0();
        for j in 0..n {
            if j == i {
                continue;
            }
            let d1 = profile[j].p1() - profile[j].p0();
            let c = d0.cross(d1);
            if c.length_squared() > TOLERANCE_LEN_MIN * TOLERANCE_LEN_MIN {
                cap_n = c.normalize();
                break;
            }
        }
        if cap_n.length_squared() > 0.5 {
            break;
        }
    }
    if cap_n.length_squared() < 0.5 {
        return Err(super::features::FeatureError::InvalidInput("profile segments are collinear"));
    }
    // The cap normal must point against the extrusion (the source cap is the
    // "bottom"): OCCT MakePrism keeps the source face orientation, and rcad's
    // make_prism_from_face_brep uses cap_n = -extr_dir with the source cap
    // FORWARD. If the profile winding gives cap_n = +extr_dir, reverse it so
    // the source cap's outward normal is -extr_dir.
    if cap_n.dot(dir) > 0.0 {
        cap_n = -cap_n;
    }

    let mut brep = topods::BRep::new();
    // Top cap: same TShapes + Location (OCCT BRepSweep_Trsf::Process).
    let loc = brep.add_location(glam::DAffine3::from_translation(dir * depth));
    let rev = |sr: topods::Shape| topods::Shape {
        orientation: Orientation::Reversed,
        ..sr
    };
    let rev_v = |v: &topods::Shape| topods::Shape {
        orientation: Orientation::Reversed,
        ..v.clone()
    };

    // Profile vertices: the segment start points.
    let v: Vec<topods::Shape> = profile
        .iter()
        .map(|s| brep.add_tvertex(s.p0()))
        .collect();
    let ve: Vec<topods::Shape> = v
        .iter()
        .map(|s| topods::Shape {
            location: loc,
            ..s.clone()
        })
        .collect();

    // Bottom edges (profile segments) and top edges (located copies).
    let mut b_ed: Vec<topods::Shape> = Vec::with_capacity(n);
    for (i, seg) in profile.iter().enumerate() {
        let j = (i + 1) % n;
        let curve = match seg {
            ProfileSegment::Line { p0, p1 } => {
                let d = *p1 - *p0;
                let len = d.length();
                Some(Curve3::Line(Line3 {
                    origin: *p0,
                    direction: if len > TOLERANCE_LEN_MIN { d / len } else { DVec3::X },
                }))
            }
            ProfileSegment::Arc {
                center,
                normal,
                x_dir,
                y_dir,
                radius,
                ..
            } => Some(Curve3::Circle(Circle3 {
                center: *center,
                normal: normal.normalize_or_zero(),
                x_dir: x_dir.normalize_or_zero(),
                y_dir: y_dir.normalize_or_zero(),
                radius: *radius,
            })),
        };
        let range = match seg {
            ProfileSegment::Line { p0, p1 } => {
                let len = (*p1 - *p0).length();
                [0.0, len]
            }
            ProfileSegment::Arc { t0, t1, .. } => [*t0, *t1],
        };
        b_ed.push(brep.add_tedge(
            curve,
            v[i].clone(),
            rev_v(&v[j]),
            range,
        ));
    }
    let t_ed: Vec<topods::Shape> = b_ed
        .iter()
        .map(|e| topods::Shape {
            location: loc,
            ..e.clone()
        })
        .collect();

    // Vertical sweep edges.
    let mut e_ver: Vec<topods::Shape> = Vec::with_capacity(n);
    for i in 0..n {
        let d = dir * depth;
        let len = d.length();
        e_ver.push(brep.add_tedge(
            Some(Curve3::Line(Line3 {
                origin: profile[i].p0(),
                direction: if len > TOLERANCE_LEN_MIN { dir } else { DVec3::X },
            })),
            v[i].clone(),
            rev_v(&ve[i]),
            [0.0, len],
        ));
    }

    // Side faces.
    let mut faces: Vec<topods::Shape> = Vec::with_capacity(n + 2);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = profile[i].p0();
        let b = profile[i].p1();
        // Outward normal from the profile centroid to the segment midpoint.
        let centroid = {
            let mut sum = DVec3::ZERO;
            for s in profile.iter() {
                sum += s.p0();
            }
            sum / n as f64
        };
        let mid = (a + b) * 0.5;
        let outward = (mid - centroid).normalize_or_zero();
        let wire = brep.add_twire(vec![
            b_ed[i].clone(),
            e_ver[j].clone(),
            rev(t_ed[i].clone()),
            rev(e_ver[i].clone()),
        ]);
        let (surf, side_reversed) = match &profile[i] {
            ProfileSegment::Line { .. } => {
                // Swept line: plane with u = edge direction, v = sweep
                // direction (the base edge at v=0, the extruded copy at
                // v=depth), matching OCCT BRepSweep lateral faces. The
                // natural normal u x v may point inward; the face is then
                // stored REVERSED (OCCT BRepPrim_GWedge stores the MIN faces
                // reversed) so the plane frame stays consistent: u x v ==
                // normal.
                let u_dir = (b - a).normalize_or_zero();
                let v_dir = dir;
                let n_nat = u_dir.cross(v_dir).normalize_or_zero();
                let inward = n_nat.dot(outward) < 0.0;
                let n_out = if inward { -n_nat } else { n_nat };
                (
                    Surface3::Plane(Plane {
                        origin: a,
                        normal: n_out,
                        u_dir,
                        v_dir,
                    }),
                    inward,
                )
            }
            ProfileSegment::Arc {
                center,
                x_dir,
                radius,
                ..
            } => {
                // Swept arc: cylinder through the arc center. ref_dir = x_dir
                // aligns u=0 with the profile frame; flip it when the cylinder
                // radial at the segment midpoint opposes the outward direction.
                let ref_dir = {
                    let rm = (mid - *center).normalize_or_zero();
                    let xd = x_dir.normalize_or_zero();
                    if rm.dot(outward) < 0.0 { -xd } else { xd }
                };
                (
                    Surface3::Cylinder(CylindricalSurface {
                        origin: *center,
                        axis: dir,
                        radius: *radius,
                        ref_dir,
                    }),
                    false,
                )
            }
        };
        let face = brep.add_tface(Some(surf), wire, vec![], None, None, vec![], false);
        faces.push(if side_reversed { rev(face) } else { face });
    }

    // Caps.
    // Profile-plane frame: u/v directions in the plane (u from the first
    // segment direction).
    let u_dir = {
        let d = profile[0].p1() - profile[0].p0();
        d.normalize_or_zero()
    };
    let v_dir = cap_n.cross(u_dir).normalize_or_zero();
    // Source (bottom) cap: outward -extr_dir (= cap_n), FORWARD.
    let bot_wire = brep.add_twire(b_ed.clone());
    faces.push(brep.add_tface(
        Some(Surface3::Plane(Plane {
            origin: profile[0].p0(),
            normal: cap_n,
            u_dir,
            v_dir,
        })),
        bot_wire,
        vec![],
        None,
        None,
        vec![],
        false,
    ));
    // Top cap: REVERSED + Location (OCCT reuses the source face TShape).
    let top_wire = brep.add_twire(t_ed.iter().map(|e| rev(e.clone())).collect());
    faces.push(rev(brep.add_tface(
        Some(Surface3::Plane(Plane {
            origin: profile[0].p0() + dir * depth,
            normal: cap_n,
            u_dir,
            v_dir,
        })),
        top_wire,
        vec![],
        None,
        None,
        vec![],
        false,
    )));

    let shell = brep.add_tshell(faces);
    brep.add_tsolid(vec![shell]);
    Ok(brep)
}

