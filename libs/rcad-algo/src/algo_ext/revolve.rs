//! Revolve a closed planar polygon around an axis into a solid (topods).
//!
//! Replaces the legacy `rcad_modeling::revolve` (removed with the old builder
//! API). For a full turn the lateral faces are exact analytic surfaces
//! (Plane/Cylinder/Cone) matching OCCT `BRepPrimAPI_MakeRevol` output
//! (BRepSweep_Revol structure: one face per non-axis profile edge, shared
//! closed circle edges at each swept vertex, seam = the profile edge itself).
//! A partial sweep falls back to generic `RevolutionSurface` faces.
//! The polygon must lie in a plane containing the axis (the OCCT `revol`
//! usage).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Circle3, ConicalSurface, Curve2d, Curve3, CylindricalSurface, Line2d, Line3, Plane,
    RevolutionSurface, Surface3,
};
use rcad_kernel::topods::{self, CurveRepresentation, Orientation};

/// Rotate `p` around the axis `(origin, dir)` by `angle` radians.
fn rotate_point(p: DVec3, origin: DVec3, dir: DVec3, angle: f64) -> DVec3 {
    let v = p - origin;
    let v_para = dir * v.dot(dir);
    let v_perp = v - v_para;
    let (s, c) = angle.sin_cos();
    v_para + v_perp * c + dir.cross(v_perp) * s + origin
}

/// Revolve a closed planar polygon (coplanar with the axis) around an axis.
///
/// `angle_rad` in `(0, 2π]`. For `2π` the start/end profile coincide and the
/// result is a closed solid of revolution without caps.
pub fn revolve_polygon(
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    axis_dir: DVec3,
    angle_rad: f64,
) -> Result<topods::BRep, String> {
    let n = profile_verts.len();
    if n < 3 {
        return Err("profile_verts needs >= 3 vertices".into());
    }
    if !axis_origin.is_finite() || !axis_dir.is_finite() {
        return Err("axis must be finite".into());
    }
    let dir = axis_dir.normalize_or_zero();
    if dir.length_squared() < 1e-24 {
        return Err("axis_dir must be non-zero".into());
    }
    if !(angle_rad > 0.0 && angle_rad <= std::f64::consts::TAU + 1e-9) {
        return Err("angle_rad must be in (0, 2*PI]".into());
    }
    let angle = angle_rad.min(std::f64::consts::TAU);

    let full_turn = (angle - std::f64::consts::TAU).abs() < 1e-9;
    if full_turn {
        revolve_polygon_full_turn(profile_verts, axis_origin, dir)
    } else {
        revolve_polygon_partial(profile_verts, axis_origin, dir, angle)
    }
}

/// Full-turn revolution with OCCT-form analytic lateral faces.
///
/// OCCT `BRepPrimAPI_MakeRevol` (BRepSweep_Revol) result structure:
/// - one vertex per non-axis profile vertex (the axis vertices carry no edge);
/// - one closed circle edge per non-axis vertex (range [0, 2π]);
/// - one profile edge (seam) per non-axis, non-radial profile edge;
/// - one face per non-axis profile edge:
///   radial edge -> Plane (disk, or annulus with an inner wire),
///   edge parallel to the axis -> Cylinder,
///   diagonal edge -> Cone;
/// - the axis profile edge generates no face.
pub fn revolve_polygon_full_turn(
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    dir: DVec3,
) -> Result<topods::BRep, String> {
    const EPS: f64 = 1e-12;
    let n = profile_verts.len();
    let mut brep = topods::BRep::new();

    // Per-vertex axis projections, radii and radial x-directions.
    let centers: Vec<DVec3> = profile_verts
        .iter()
        .map(|&p| axis_origin + dir * (p - axis_origin).dot(dir))
        .collect();
    let radii: Vec<f64> = profile_verts
        .iter()
        .zip(&centers)
        .map(|(&p, &o)| (p - o).length())
        .collect();

    // Vertex shapes and closed circle edges, only for non-axis vertices.
    let mut vtx: Vec<Option<topods::Shape>> = vec![None; n];
    let mut circ: Vec<Option<topods::Shape>> = vec![None; n];
    for i in 0..n {
        if radii[i] <= EPS {
            continue;
        }
        let v = brep.add_tvertex(profile_verts[i]);
        let x_dir = (profile_verts[i] - centers[i]) / radii[i];
        let y_dir = dir.cross(x_dir).normalize();
        let curve = Curve3::Circle(Circle3 {
            center: centers[i],
            normal: dir,
            x_dir,
            y_dir,
            radius: radii[i],
        });
        // Closed circle: both endpoints are the same vertex (OCCT seam-vertex
        // convention, matching make_cylinder_brep).  The last vertex reference
        // is REVERSED (BRep_Builder::Add stores [V1(FWD), V2(REV)]) — the
        // WireSplitter's EdgeInfo in-flag reads the first vertex orientation,
        // so a REVERSED copy of this edge must iterate as [FWD, REV] to be
        // classifiable as an OUT edge (looped) rather than an IN edge
        // (relegated to an internal wire).
        let rev_v = |s: &topods::Shape| topods::Shape {
            data: s.data.clone(),
            index: s.index,
            orientation: Orientation::Reversed,
            location: s.location,
        };
        let e = brep.add_tedge(Some(curve), v.clone(), rev_v(&v), [0.0, std::f64::consts::TAU]);
        vtx[i] = Some(v);
        circ[i] = Some(e);
    }

    // Profile (seam) edges for the non-radial edges between two non-axis
    // vertices.  OCCT BRepSweep_Revol keeps the profile edge itself as the
    // seam of its swept face (u = 0 and u = 2π).  A radial edge sweeps a
    // planar face bounded by the vertex circles only — the profile edge is
    // NOT an edge of the result (OCCT cyla dump: no radial profile edges).
    let mut seam: Vec<Option<topods::Shape>> = vec![None; n];
    for i in 0..n {
        let j = (i + 1) % n;
        if radii[i] <= EPS || radii[j] <= EPS {
            continue;
        }
        let d = profile_verts[j] - profile_verts[i];
        let len = d.length();
        if len <= EPS {
            continue;
        }
        let d_n = d / len;
        // True radial edge (perpendicular to the axis through the axis):
        // no seam edge.
        let d_perp = d_n.cross(dir).length();
        let line_axis_dist = if d_perp > EPS {
            ((profile_verts[i] - axis_origin).dot(d_n.cross(dir).normalize())).abs()
        } else {
            (profile_verts[i] - axis_origin - dir * (profile_verts[i] - axis_origin).dot(dir)).length()
        };
        if line_axis_dist <= EPS && d_perp > 1.0 - 1e-9 {
            continue;
        }
        let curve = Curve3::Line(Line3 {
            origin: profile_verts[i],
            direction: d_n,
        });
        let e = brep.add_tedge(
            Some(curve),
            vtx[i].clone().unwrap(),
            vtx[j].clone().unwrap(),
            [0.0, len],
        );
        seam[i] = Some(e);
    }

    let rev = |sr: &topods::Shape| topods::Shape {
        data: sr.data.clone(),
        index: sr.index,
        orientation: Orientation::Reversed,
        location: sr.location,
    };

    // One face per non-axis profile edge.
    let mut faces: Vec<topods::Shape> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        if radii[i] <= EPS && radii[j] <= EPS {
            // Axis edge: no swept face (OCCT drops it).
            continue;
        }
        let (ri, rj) = (radii[i], radii[j]);
        // Classify the swept surface from the profile LINE geometry (OCCT
        // BRepSweep_Revol): a line through the axis sweeps a Plane (when
        // perpendicular) or a Cone (at an angle); a line parallel to the axis
        // sweeps a Cylinder; a line skew to the axis sweeps a general
        // revolution surface.
        let d = profile_verts[j] - profile_verts[i];
        let d_n = d.normalize_or_zero();
        let d_axis = d_n.dot(dir).abs();
        let d_perp = d_n.cross(dir).length();
        let line_axis_dist = if d_perp > EPS {
            // Skew distance between the line and the axis.
            ((profile_verts[i] - axis_origin).dot(d_n.cross(dir).normalize())).abs()
        } else {
            // Line parallel to the axis: distance of a point to the axis.
            (profile_verts[i] - axis_origin - dir * (profile_verts[i] - axis_origin).dot(dir)).length()
        };
        let hits_axis = line_axis_dist <= EPS;

        if hits_axis && d_perp > 1.0 - 1e-9 {
            // Radial edge (perpendicular to the axis through the axis): a
            // planar face bounded by the two vertex circles only (the profile
            // edge itself is not an edge of the face).
            let outer = if ri >= rj { i } else { j };
            let inner = if ri >= rj { j } else { i };
            let p0 = profile_verts[i];
            let p1 = profile_verts[j];
            // Plane frame: u along the reversed edge direction, v follows the
            // OCCT cyla dump (bottom/step/top planes).
            let u_dir = (p0 - p1).normalize_or_zero();
            // OCCT frame: v = u × dir when the edge runs outward (away from
            // the axis), v = dir × u when it runs inward.
            let outward = rj > ri;
            let v_dir = if outward {
                u_dir.cross(dir).normalize()
            } else {
                dir.cross(u_dir).normalize()
            };
            // OCCT gp_Pln is always right-handed (v = normal × u). For the
            // outward edge the frame above gives u × v = -dir, so the stored
            // normal must be u × v, not dir — otherwise the plane is
            // left-handed and the face's effective normal is flipped relative
            // to OCCT (K1: the annulus SD face then fails the
            // IsSplitToReverse check against the box top and the box's
            // SplitSolid shell cannot close).
            let plane = Surface3::Plane(Plane {
                origin: centers[i],
                normal: u_dir.cross(v_dir).normalize(),
                u_dir,
                v_dir,
            });
            let outer_edge = circ[outer].clone().unwrap();
            if radii[inner] <= EPS {
                // Disk: single outer circle.  OCCT orientation: reversed when
                // the edge runs from the axis/inner side outward, forward when
                // it runs inward (bottom disk [-circle], top disk [+circle]).
                let outer_fwd = ri > rj; // edge direction outer -> inner
                let wire = if outer_fwd {
                    brep.add_twire(vec![outer_edge])
                } else {
                    brep.add_twire(vec![rev(&outer_edge)])
                };
                faces.push(brep.add_tface(Some(plane), wire, vec![], None, None, vec![], false));
            } else {
                // Annulus: the swept face of a radial profile edge. OCCT
                // BRepSweep_Revol (K1 runner dump: outer=[r20 F], inner=[r25 R]
                // for the outward edge (50,80,100)->(50,75,100)) stores the
                // START-vertex circle as the outer wire FORWARD and the
                // END-vertex circle as the inner wire REVERSED — the wire roles
                // follow the profile-edge vertex order, not the radius. This
                // opposes the shared circles to the adjacent lateral walls
                // (r20cyl top rim R vs annulus outer F; annulus inner R vs
                // r25cyl bottom rim F) so the shell is traversable.
                let outer_wire = brep.add_twire(vec![circ[i].clone().unwrap()]);
                let inner_wire = brep.add_twire(vec![rev(&circ[j].clone().unwrap())]);
                faces.push(brep.add_tface(
                    Some(plane),
                    outer_wire,
                    vec![inner_wire],
                    None,
                    None,
                    vec![],
                    false,
                ));
            }
        } else if d_axis > 1.0 - 1e-9 && ri > EPS && rj > EPS {
            // Edge parallel to the axis: swept surface is a cylinder.
            let (lo, hi) = if (centers[i] - axis_origin).dot(dir) <= (centers[j] - axis_origin).dot(dir)
            {
                (i, j)
            } else {
                (j, i)
            };
            let x_dir = (profile_verts[lo] - centers[lo]) / radii[lo];
            let surf = Surface3::Cylinder(CylindricalSurface {
                origin: centers[lo],
                axis: dir,
                radius: radii[lo],
                ref_dir: x_dir,
            });
            let seam_lo = seam[i].clone().unwrap();
            let wire = brep.add_twire(vec![
                circ[lo].clone().unwrap(),
                rev(&circ[hi].clone().unwrap()),
                seam_lo.clone(),
                rev(&seam_lo),
            ]);
            let face = brep.add_tface(Some(surf), wire, vec![], None, None, vec![], false);
            // OCCT BRepSweep_Revol: the seam (u = 0 and u = 2π) is a closed
            // edge on the lateral face — it carries two pcurves stored as a
            // CurveOnClosedSurface representation (matching make_cylinder_brep
            // BRepPrim_OneAxis::LateralFace L434-438). Without it the pipeline
            // treats the seam as a free boundary and avoids the face.
            let lat_key = (face.ptr_id(), face.location);
            let seam_len = (profile_verts[j] - profile_verts[i]).length();
            let (pc1, pc2) = if lo == i {
                (Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, 0.0), DVec2::Y)),
                 Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)))
            } else {
                (Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, seam_len), DVec2::new(0.0, -1.0))),
                 Curve2d::Line(Line2d::new(DVec2::new(0.0, seam_len), DVec2::new(0.0, -1.0))))
            };
            brep.edge_mut_inplace(seam_lo.clone()).pcurves.insert(lat_key, (pc1.clone(), 0.0, seam_len));
            brep.edge_mut_inplace(seam_lo.clone()).representations
                .push(rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    face: lat_key,
                    pcurve1: pc1,
                    pcurve2: pc2,
                    range: [0.0, seam_len],
                });
            // Circle edges are V-isolines on the lateral face.
            brep.edge_mut_inplace(circ[lo].clone().unwrap()).pcurves.insert(
                lat_key,
                (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X)), 0.0, std::f64::consts::TAU),
            );
            brep.edge_mut_inplace(circ[hi].clone().unwrap()).pcurves.insert(
                lat_key,
                (Curve2d::Line(Line2d::new(DVec2::new(0.0, seam_len), DVec2::X)), 0.0, std::f64::consts::TAU),
            );
            faces.push(face);
        } else if hits_axis && ri > EPS && rj > EPS {
            // Diagonal edge through the axis: swept surface is a cone.
            let (lo, hi) = if (centers[i] - axis_origin).dot(dir) <= (centers[j] - axis_origin).dot(dir)
            {
                (i, j)
            } else {
                (j, i)
            };
            let dz = (centers[hi] - centers[lo]).length();
            let dr = radii[hi] - radii[lo];
            let half_angle = if dz.abs() > EPS && dr.abs() > EPS {
                dr.atan2(dz)
            } else {
                0.0
            };
            let x_dir = (profile_verts[lo] - centers[lo]) / radii[lo];
            let surf = Surface3::Cone(ConicalSurface {
                apex: centers[lo],
                axis: dir,
                radius: radii[lo],
                half_angle_rad: half_angle,
                ref_dir: x_dir,
            });
            let seam_lo = seam[i].clone().unwrap();
            let wire = brep.add_twire(vec![
                circ[lo].clone().unwrap(),
                rev(&circ[hi].clone().unwrap()),
                seam_lo.clone(),
                rev(&seam_lo),
            ]);
            let face = brep.add_tface(Some(surf), wire, vec![], None, None, vec![], false);
            // Seam + V-isoline pcurves on the cone lateral (same convention as
            // the cylinder branch above).
            let lat_key = (face.ptr_id(), face.location);
            let seam_len = (profile_verts[j] - profile_verts[i]).length();
            let (pc1, pc2) = if lo == i {
                (Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, 0.0), DVec2::Y)),
                 Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::Y)))
            } else {
                (Curve2d::Line(Line2d::new(DVec2::new(std::f64::consts::TAU, seam_len), DVec2::new(0.0, -1.0))),
                 Curve2d::Line(Line2d::new(DVec2::new(0.0, seam_len), DVec2::new(0.0, -1.0))))
            };
            brep.edge_mut_inplace(seam_lo.clone()).pcurves.insert(lat_key, (pc1.clone(), 0.0, seam_len));
            brep.edge_mut_inplace(seam_lo.clone()).representations
                .push(rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                    face: lat_key,
                    pcurve1: pc1,
                    pcurve2: pc2,
                    range: [0.0, seam_len],
                });
            brep.edge_mut_inplace(circ[lo].clone().unwrap()).pcurves.insert(
                lat_key,
                (Curve2d::Line(Line2d::new(DVec2::new(0.0, 0.0), DVec2::X)), 0.0, std::f64::consts::TAU),
            );
            brep.edge_mut_inplace(circ[hi].clone().unwrap()).pcurves.insert(
                lat_key,
                (Curve2d::Line(Line2d::new(DVec2::new(0.0, seam_len), DVec2::X)), 0.0, std::f64::consts::TAU),
            );
            faces.push(face);
        } else {
            // Skew line (or a diagonal through the axis with an on-axis
            // endpoint): general revolution surface.  An on-axis endpoint has
            // no circle edge; use a degenerate point edge so wires close.
            let mut e_i = circ[i].clone();
            if e_i.is_none() {
                let v = brep.add_tvertex(profile_verts[i]);
                e_i = Some(brep.add_tedge(None, v.clone(), v.clone(), [0.0, 0.0]));
            }
            let mut e_j = circ[j].clone();
            if e_j.is_none() {
                let v = brep.add_tvertex(profile_verts[j]);
                e_j = Some(brep.add_tedge(None, v.clone(), v.clone(), [0.0, 0.0]));
            }
            let seam_lo = match seam[i].clone() {
                Some(s) => s,
                None => {
                    // On-axis endpoint: build the seam from the two vertices.
                    let v0 = vtx[i].clone().or_else(|| Some(brep.add_tvertex(profile_verts[i]))).unwrap();
                    let v1 = vtx[j].clone().or_else(|| Some(brep.add_tvertex(profile_verts[j]))).unwrap();
                    let len = d.length();
                    brep.add_tedge(Some(Curve3::Line(Line3 {
                        origin: profile_verts[i],
                        direction: d_n,
                    })), v0, v1, [0.0, len])
                }
            };
            let profile_line = Curve3::Line(Line3 {
                origin: profile_verts[i],
                direction: d_n,
            });
            let surface = Surface3::Revolution(RevolutionSurface {
                profile: Box::new(profile_line),
                axis_origin,
                axis_dir: dir,
            });
            let wire = brep.add_twire(vec![
                e_i.unwrap(),
                rev(&e_j.unwrap()),
                seam_lo.clone(),
                rev(&seam_lo),
            ]);
            faces.push(brep.add_tface(Some(surface), wire, vec![], None, None, vec![], false));
        }
    }

    let shell = brep.add_tshell(faces);
    brep.add_tsolid(vec![shell]);
    Ok(brep)
}

/// Partial revolution (angle < 2π) with OCCT BRepSweep_Revol-form faces:
/// analytic Plane/Cylinder/Cone laterals (RevolutionSurface fallback for skew
/// lines) and two planar caps.  The end profile copies are the SAME TShapes
/// carried with a rotation Location (index `l1`), matching OCCT's TopLoc
/// structure (the end cap is the start cap shape at L1, reversed); this keeps
/// the vertex count at the profile-vertex count (nbshapes V=4 for the ring).
pub fn revolve_polygon_partial(
    profile_verts: &[DVec3],
    axis_origin: DVec3,
    dir: DVec3,
    angle: f64,
) -> Result<topods::BRep, String> {
    const EPS: f64 = 1e-12;
    let n = profile_verts.len();
    let mut brep = topods::BRep::new();

    // Per-vertex axis projections and radii.
    let centers: Vec<DVec3> = profile_verts
        .iter()
        .map(|&p| axis_origin + dir * (p - axis_origin).dot(dir))
        .collect();
    let radii: Vec<f64> = profile_verts
        .iter()
        .zip(&centers)
        .map(|(&p, &o)| (p - o).length())
        .collect();

    // The rotation Location L1 (OCCT TopLoc_Location of the end profile).
    let rot = glam::DAffine3::from_translation(axis_origin)
        * glam::DAffine3::from_axis_angle(dir, angle)
        * glam::DAffine3::from_translation(-axis_origin);
    let l1 = brep.add_location(rot);

    // Profile vertices (start positions); the end copies reference the same
    // TShapes with Location L1.
    let vtx: Vec<topods::Shape> = profile_verts
        .iter()
        .map(|&p| brep.add_tvertex(p))
        .collect();
    let with_l1 = |s: &topods::Shape| topods::Shape::from_parts(
        s.data.clone(), s.index, l1, s.orientation,
    );
    let rev = |sr: &topods::Shape| topods::Shape {
        data: sr.data.clone(),
        index: sr.index,
        orientation: Orientation::Reversed,
        location: sr.location,
    };

    // Arc edges: vertex i sweeps from u=0 to u=angle (range [0, angle]).
    // OCCT BRep_Builder::Add stores the edge endpoints as [V1(FWD), V2(REV)]
    // (BRep_Builder.cxx Add(TopoDS_Edge, TopoDS_Vertex, ...) — the second
    // vertex is normalized to REVERSED), so the swept (L1) endpoint is passed
    // reversed, matching make_cylinder_brep.
    let mut arcs: Vec<Option<topods::Shape>> = vec![None; n];
    for i in 0..n {
        let r = radii[i];
        if r < EPS {
            arcs[i] = Some(brep.add_tedge(None, vtx[i].clone(), rev(&with_l1(&vtx[i])), [0.0, 0.0]));
            continue;
        }
        let center = centers[i];
        let x_dir = (profile_verts[i] - center) / r;
        let y_dir = dir.cross(x_dir).normalize();
        let curve = Curve3::Circle(Circle3 {
            center,
            normal: dir,
            x_dir,
            y_dir,
            radius: r,
        });
        arcs[i] = Some(brep.add_tedge(Some(curve), vtx[i].clone(), rev(&with_l1(&vtx[i])), [0.0, angle]));
    }

    // Profile edges (lines at the start position); the u=angle copies are the
    // same edges with Location L1.
    let mut prof: Vec<Option<topods::Shape>> = vec![None; n];
    for i in 0..n {
        let j = (i + 1) % n;
        let d = profile_verts[j] - profile_verts[i];
        let len = d.length();
        if len <= EPS {
            continue;
        }
        let curve = Curve3::Line(Line3 {
            origin: profile_verts[i],
            direction: d / len,
        });
        prof[i] = Some(brep.add_tedge(Some(curve), vtx[i].clone(), rev(&vtx[j]), [0.0, len]));
    }

    // Lateral faces: one per profile edge.  Wire: [arc_i FWD, arc_j REV,
    // prof_i (u=0), prof_i (u=angle, Location L1)] — the OCCT ring dump
    // bottom-annulus wire [+23 -21 -20 +20(L1)].
    let mut faces: Vec<topods::Shape> = Vec::new();
    for i in 0..n {
        let j = (i + 1) % n;
        if radii[i] <= EPS && radii[j] <= EPS {
            continue;
        }
        let (ri, rj) = (radii[i], radii[j]);
        let d = profile_verts[j] - profile_verts[i];
        let d_n = d.normalize_or_zero();
        let d_axis = d_n.dot(dir).abs();
        let d_perp = d_n.cross(dir).length();
        let line_axis_dist = if d_perp > EPS {
            ((profile_verts[i] - axis_origin).dot(d_n.cross(dir).normalize())).abs()
        } else {
            (profile_verts[i] - axis_origin - dir * (profile_verts[i] - axis_origin).dot(dir)).length()
        };
        let hits_axis = line_axis_dist <= EPS;

        let prof_i = prof[i].clone().unwrap();
        // OCCT BRepSweep_Revol lateral wire: [arc_START FWD, arc_END REV,
        // prof(u=0) REV, prof(u=angle, L1) FWD] (the OCCT ring dump
        // "+23 -21 -20 +20(L1)"). The u=0 profile edge is REVERSED so it
        // opposes the cap face's FORWARD profile edge along the shared
        // boundary (the ShellSplitter's GetEdgeOff needs the reversed
        // orientation to connect the faces into a closed shell).
        let wire = brep.add_twire(vec![
            arcs[i].clone().unwrap(),
            rev(&arcs[j].clone().unwrap()),
            rev(&prof_i),
            with_l1(&prof_i),
        ]);

        let surface = if hits_axis && d_perp > 1.0 - 1e-9 {
            // Radial edge: partial planar sector. Right-handed normal (v =
            // normal × u); the outward edge's frame gives u × v = -dir, so the
            // stored normal must be u × v (see revolve_polygon_full_turn).
            let u_dir = (profile_verts[i] - profile_verts[j]).normalize_or_zero();
            let outward = rj > ri;
            let v_dir = if outward {
                u_dir.cross(dir).normalize()
            } else {
                dir.cross(u_dir).normalize()
            };
            Surface3::Plane(Plane {
                origin: centers[i],
                normal: u_dir.cross(v_dir).normalize(),
                u_dir,
                v_dir,
            })
        } else if d_axis > 1.0 - 1e-9 && ri > EPS && rj > EPS {
            // Edge parallel to the axis: partial cylinder.
            Surface3::Cylinder(CylindricalSurface {
                origin: centers[i],
                axis: dir,
                radius: ri,
                ref_dir: (profile_verts[i] - centers[i]) / ri,
            })
        } else if hits_axis && ri > EPS && rj > EPS {
            // Diagonal through the axis: partial cone.
            let dz = (centers[j] - centers[i]).length();
            let dr = rj - ri;
            let half_angle = if dz.abs() > EPS && dr.abs() > EPS {
                dr.atan2(dz)
            } else {
                0.0
            };
            Surface3::Cone(ConicalSurface {
                apex: centers[i],
                axis: dir,
                radius: ri,
                half_angle_rad: half_angle,
                ref_dir: (profile_verts[i] - centers[i]) / ri,
            })
        } else {
            // Skew line: general revolution surface.
            Surface3::Revolution(RevolutionSurface {
                profile: Box::new(Curve3::Line(Line3 {
                    origin: profile_verts[i],
                    direction: d_n,
                })),
                axis_origin,
                axis_dir: dir,
            })
        };
        faces.push(brep.add_tface(Some(surface), wire, vec![], None, None, vec![], false));
    }

    // Caps: the profile polygon at u=0 (forward) and at u=angle (reversed,
    // Location L1 — OCCT -3(L1)).  The end cap is a separate face whose wire
    // uses the profile edges with Location L1.
    let normal = {
        let a = profile_verts[1] - profile_verts[0];
        let b = profile_verts[2] - profile_verts[0];
        let n = a.cross(b);
        if n.length_squared() > 1e-24 {
            n.normalize()
        } else {
            dir
        }
    };
    let start_wire = brep.add_twire(prof.iter().map(|p| p.clone().unwrap()).collect());
    let start_cap = brep.add_tface(
        Some(Surface3::Plane(Plane::new(profile_verts[0], -normal))),
        start_wire,
        vec![],
        None,
        None,
        vec![],
        false,
    );
    faces.push(start_cap.clone());

    // The end cap is the START cap shape carried at the rotation Location L1,
    // REVERSED (OCCT BRepSweep_Revol ring dump: the end profile is the start
    // profile at L1 with REVERSED orientation, "face or=1" in the runner).
    // Using the same TFace at a different Location keeps the face count at the
    // profile-face count (nbshapes F=5 for the ring: 2 walls + 2 annuli +
    // 1 cap), matching OCCT.
    faces.push(rev(&with_l1(&start_cap)));

    let shell = brep.add_tshell(faces);
    brep.add_tsolid(vec![shell]);
    Ok(brep)
}
