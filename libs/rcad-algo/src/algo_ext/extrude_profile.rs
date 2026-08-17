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

/// OCCT `gp_Ax3(P, V)` X/Y axes of the right-handed frame with Z = V
/// (gp_Ax3.cxx L29-80).
fn gp_ax3_axes(v: DVec3) -> (DVec3, DVec3) {
    let a = v.x.abs();
    let b = v.y.abs();
    let c = v.z.abs();
    let x = if b <= a && b <= c {
        if a > c {
            DVec3::new(-v.z, 0.0, v.x)
        } else {
            DVec3::new(v.z, 0.0, -v.x)
        }
    } else if a <= b && a <= c {
        if b > c {
            DVec3::new(0.0, -v.z, v.y)
        } else {
            DVec3::new(0.0, v.z, -v.y)
        }
    } else if a > b {
        DVec3::new(-v.y, v.x, 0.0)
    } else {
        DVec3::new(v.y, -v.x, 0.0)
    };
    let x = x.normalize_or_zero();
    (x, v.cross(x).normalize_or_zero())
}

/// Sample points with weights of one profile segment, mirroring
/// BRepLib_FindSurface fillParams/fillPoints (BRepLib_FindSurface.cxx
/// L178-242): a line yields its two ends with weight = length each; an arc
/// yields four points at t0, t0+d/3, t0+2d/3, t1 with weights = sum of the
/// adjacent chord distances.
fn find_surface_samples(seg: &ProfileSegment) -> Vec<(DVec3, f64)> {
    match seg {
        ProfileSegment::Line { p0, p1 } => {
            let len = (*p1 - *p0).length();
            vec![(*p0, len), (*p1, len)]
        }
        ProfileSegment::Arc {
            center,
            x_dir,
            y_dir,
            radius,
            t0,
            t1,
            ..
        } => {
            let pt = |t: f64| *center + *radius * (t.cos() * *x_dir + t.sin() * *y_dir);
            let d = *t1 - *t0;
            let p = [pt(*t0), pt(*t0 + d / 3.0), pt(*t0 + 2.0 * d / 3.0), pt(*t1)];
            let c01 = (p[1] - p[0]).length();
            let c12 = (p[2] - p[1]).length();
            let c23 = (p[3] - p[2]).length();
            vec![
                (p[0], c01),
                (p[1], c01 + c12),
                (p[2], c12 + c23),
                (p[3], c23),
            ]
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
    // Cap frame from OCCT BRepLib_FindSurface (BRepLib_FindSurface.cxx
    // L248-578), as used by BRepLib_MakeFace/BRepPrimAPI_MakePrism: origin =
    // weighted barycenter of the sampled profile points, normal = profile
    // plane normal oriented so the wire is CCW seen from +normal, u/v = the
    // gp_Ax3(P, V) axes.
    let mut sum = DVec3::ZERO;
    let mut wsum = 0.0;
    let mut pts: Vec<DVec3> = Vec::new();
    for seg in profile.iter() {
        for (p, w) in find_surface_samples(seg) {
            sum += w * p;
            wsum += w;
            pts.push(p);
        }
    }
    let bary = if wsum > 0.0 { sum / wsum } else { DVec3::ZERO };
    // n_raw: profile plane normal from the first two non-collinear segment
    // chords (the covariance eigenvector direction for planar points).
    let mut n_raw = DVec3::ZERO;
    for i in 0..n {
        let d0 = profile[i].p1() - profile[i].p0();
        for j in 0..n {
            if j == i {
                continue;
            }
            let d1 = profile[j].p1() - profile[j].p0();
            let c = d0.cross(d1);
            if c.length_squared() > TOLERANCE_LEN_MIN * TOLERANCE_LEN_MIN {
                n_raw = c.normalize();
                break;
            }
        }
        if n_raw.length_squared() > 0.5 {
            break;
        }
    }
    if n_raw.length_squared() < 0.5 {
        return Err(super::features::FeatureError::InvalidInput("profile segments are collinear"));
    }
    // Orient n so the profile is CCW from +n (BRepLib_FindSurface L563-577:
    // the infinite point is IN for a CW wire, which reverses the normal).
    let (u0, v0) = gp_ax3_axes(n_raw);
    let mut area = 0.0;
    for i in 0..pts.len() {
        let a = pts[i] - bary;
        let b = pts[(i + 1) % pts.len()] - bary;
        area += (a.dot(u0)) * (b.dot(v0)) - (b.dot(u0)) * (a.dot(v0));
    }
    let cap_n = if area < 0.0 { -n_raw } else { n_raw };
    let (u_dir, v_dir) = gp_ax3_axes(cap_n);

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
        // Lateral face wire order matches OCCT BRepSweep_NumLinearRegularSweep
        // (the intermediate newWire): [directing edge at the profile edge's
        // first vertex, directing edge at the second vertex, the generating
        // edge, the translated generating edge]. OCCT (bcommon G1, BFS-SRC
        // dump of the DS y=x lateral): [sweep@V1 F, sweep@V2 R, profile R,
        // top-profile F] — the V1 directing edge keeps FORWARD (the first
        // vertex of a FORWARD edge is FORWARD), the V2 directing edge is
        // REVERSED (the last vertex), matching BuildShell L264-277
        // (Or = It.Orientation() of the profile-edge vertex).
        let wire = brep.add_twire(vec![
            e_ver[i].clone(),
            rev(e_ver[j].clone()),
            rev(b_ed[i].clone()),
            t_ed[i].clone(),
        ]);
        let surf = match &profile[i] {
            ProfileSegment::Line { .. } => {
                // Swept line: plane with u = edge direction, v = -sweep
                // direction (the base edge at v=0, the extruded copy at
                // v=-depth), matching OCCT BRepSweep_Translation::MakeEmptyFace
                // (D.Reverse(), BRepSweep_Translation.cxx L238-239) via
                // GeomAdaptor_SurfaceOfLinearExtrusion::Plane(). The surface
                // keeps the natural frame u x v == normal (OCCT stores the
                // surface without the face orientation); the shell REVERSED
                // provides the composed outward direction.
                let u_dir = (b - a).normalize_or_zero();
                let v_dir = -dir;
                let n_nat = u_dir.cross(v_dir).normalize_or_zero();
                Surface3::Plane(Plane {
                    origin: a,
                    normal: n_nat,
                    u_dir,
                    v_dir,
                })
            }
            ProfileSegment::Arc {
                center,
                normal,
                x_dir,
                radius,
                ..
            } => {
                // Swept arc: cylinder through the arc center. OCCT
                // GeomAdaptor_SurfaceOfLinearExtrusion::Cylinder() uses
                // D = -sweep as the axis (the circle axis is parallel to the
                // sweep, so ZReverse keeps the axis along D) and ZReverse
                // (which flips X) applies iff D . Z < 0, i.e. dir . normal > 0.
                let xd = x_dir.normalize_or_zero();
                let zrev = dir.dot(normal.normalize_or_zero()) > 0.0;
                let ref_dir = if zrev { -xd } else { xd };
                Surface3::Cylinder(CylindricalSurface {
                    origin: *center,
                    axis: -dir,
                    radius: *radius,
                    ref_dir,
                })
            }
        };
        let face = brep.add_tface(Some(surf), wire, vec![], None, None, vec![], false);
        // Lateral faces are stored natural (FORWARD); the shell carries the
        // REVERSED orientation (OCCT BRepSweep_NumLinearRegularSweep::BuildShell
        // adds the whole shell with ShellOri = DirectSolid). OCCT stores the
        // lateral faces FORWARD (the swept edge keeps its wire orientation) and
        // the composed (shell REVERSED x face FORWARD) points outward.
        faces.push(face);
    }

    // Caps. Bottom cap outward = -extr_dir, top cap outward = +extr_dir.
    // OCCT BRepSweep_Translation::DirectSolid (BRepSweep_Translation.cxx
    // L420-436) orients the shell REVERSED iff sweep . source-normal > 0, so
    // the composed cap orientations are REVERSED (bottom) / FORWARD (top) when
    // cap_n . dir > 0. The caps themselves are stored natural: the bottom cap
    // with the direction FIRST vertex orientation (FORWARD), the top cap with
    // the direction SECOND vertex orientation (REVERSED) — the shell REVERSED
    // then composes them to REVERSED (bottom) / FORWARD (top).
    let bot_rev = cap_n.dot(dir) > 0.0;
    let bot_wire = brep.add_twire(b_ed.clone());
    let bot_face = brep.add_tface(
        Some(Surface3::Plane(Plane {
            origin: bary,
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
    );
    faces.push(bot_face);
    // Top cap: same TShape + Location (OCCT reuses the source face TShape).
    // The translated profile wire keeps the source edge orientations (OCCT
    // BRepSweep_Trsf::Process, BRepSweep_Trsf.cxx L83-90: newShape.Move(loc) —
    // the wire is NOT reversed). OCCT's z1 cap TShape therefore stores the top
    // edges FORWARD; the FClass2d loop sampling (IntTools_FClass2d.cxx L105:
    // face forced FORWARD) then walks the top boundary CCW in UV. Storing the
    // top wire reversed would sample it CW and misclassify the split loops as
    // holes in BuilderFace::PerformAreas (bcommon G1: z1 cap -> 1 area).
    let top_wire = brep.add_twire(t_ed.clone());
    let top_face = brep.add_tface(
        Some(Surface3::Plane(Plane {
            origin: bary + dir * depth,
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
    );
    faces.push(rev(top_face));

    let mut shell = brep.add_tshell(faces);
    if bot_rev {
        shell.orientation = Orientation::Reversed;
    }
    if std::env::var("RCAD_SHELL_REV").is_ok() {
        eprintln!("[SHELL-REV] bot_rev={} shell_or={:?}", bot_rev, shell.orientation);
    }
    brep.add_tsolid(vec![shell]);
    Ok(brep)
}

