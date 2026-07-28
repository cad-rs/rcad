// extra.rs — included into mod.rs via include!("extra.rs")
// All items from above the include line in mod.rs are accessible without imports.

fn compute_planar_section_properties(polylines: &[Vec<DVec3>], plane: &Plane) -> Option<SectionProperties> {
    if polylines.is_empty() {
        return None;
    }

    let (area, centroid, ixx, iyy, ixy) = compute_polygon_properties(polylines, plane);

    let perimeter: f64 = polylines
        .iter()
        .map(|pts| {
            pts.windows(2)
                .map(|w| (w[1] - w[0]).length())
                .sum::<f64>()
        })
        .sum();

    Some(SectionProperties {
        area,
        centroid,
        ixx,
        iyy,
        ixy,
        perimeter,
    })
}

fn compute_polygon_properties(polylines: &[Vec<DVec3>], plane: &Plane) -> (f64, DVec3, f64, f64, f64) {
    let x_axis = plane.u_dir;
    let y_axis = plane.v_dir;

    let mut total_area = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;

    for pts in polylines {
        let n = pts.len();
        if n < 3 {
            continue;
        }

        let pts_2d: Vec<(f64, f64)> = pts
            .iter()
            .map(|p| {
                let v = *p - plane.origin;
                (v.dot(x_axis), v.dot(y_axis))
            })
            .collect();

        let mut signed_area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            signed_area += pts_2d[i].0 * pts_2d[j].1 - pts_2d[j].0 * pts_2d[i].1;
        }
        signed_area *= 0.5;

        total_area += signed_area;

        if signed_area.abs() > TOLERANCE_LEN_MIN {
            for i in 0..n {
                let j = (i + 1) % n;
                let factor = pts_2d[i].0 * pts_2d[j].1 - pts_2d[j].0 * pts_2d[i].1;
                cx += (pts_2d[i].0 + pts_2d[j].0) * factor;
                cy += (pts_2d[i].1 + pts_2d[j].1) * factor;
            }
        }
    }

    if total_area.abs() < TOLERANCE_LEN_MIN {
        return (0.0, plane.origin, 0.0, 0.0, 0.0);
    }

    cx /= 6.0 * total_area;
    cy /= 6.0 * total_area;

    let centroid = plane.origin + x_axis * cx + y_axis * cy;

    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut ixy = 0.0;

    for pts in polylines {
        let n = pts.len();
        if n < 3 {
            continue;
        }

        let pts_2d: Vec<(f64, f64)> = pts
            .iter()
            .map(|p| {
                let v = *p - centroid;
                (v.dot(x_axis), v.dot(y_axis))
            })
            .collect();

        for i in 0..n {
            let j = (i + 1) % n;
            let x_i = pts_2d[i].0;
            let y_i = pts_2d[i].1;
            let x_j = pts_2d[j].0;
            let y_j = pts_2d[j].1;

            let factor = x_i * y_j - x_j * y_i;

            ixx += factor * (y_i * y_i + y_i * y_j + y_j * y_j);
            iyy += factor * (x_i * x_i + x_i * x_j + x_j * x_j);
            ixy += factor * (x_i * y_i + x_i * y_j + x_j * y_i + x_j * y_j);
        }
    }

    ixx /= 12.0;
    iyy /= 12.0;
    ixy /= 24.0;

    (total_area.abs(), centroid, ixx.abs(), iyy.abs(), ixy)
}

// = =  Multiple Section Support = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// Generate multiple sections at evenly spaced planes along an axis.
pub fn section_parallel_planes(
    brep: &BRep,
    origin: DVec3,
    direction: DVec3,
    spacing: f64,
    count: usize,
) -> Vec<SectionResult> {
    let dir = direction.normalize();
    let mut results = Vec::with_capacity(count);

    for i in 0..count {
        let plane_origin = origin + dir * (spacing * i as f64);
        let plane = Plane::new(plane_origin, dir);
        results.push(section_by_plane(brep, &plane));
    }

    results
}

/// Generate cross-sections along a path curve.
pub fn section_along_path(
    brep: &BRep,
    path: &Curve3,
    param_values: &[f64],
) -> Vec<SectionResult> {
    let mut results = Vec::with_capacity(param_values.len());

    for &t in param_values {
        let origin = path.point_at(t);
        let normal = path.tangent_at(t);
        let plane = Plane::new(origin, normal);
        results.push(section_by_plane(brep, &plane));
    }

    results
}

/// Cross-section generation along a path with automatic spacing.
pub fn cross_sections_along_path(brep: &BRep, path: &Curve3, count: usize) -> Vec<SectionResult> {
    let [t0, t1] = path.default_domain();

    let (t0, t1) = if !t0.is_finite() || !t1.is_finite() {
        let center = (t0 + t1) * 0.5;
        if center.is_finite() {
            (center - 50.0, center + 50.0)
        } else {
            (-50.0, 50.0)
        }
    } else {
        (t0, t1)
    };

    let param_values: Vec<f64> = (0..count)
        .map(|i| t0 + (t1 - t0) * i as f64 / (count - 1).max(1) as f64)
        .collect();

    section_along_path(brep, path, &param_values)
}

/// Stitch multiple section wires into a lofted solid.
pub fn stitch_sections_to_solid(sections: &[SectionResult], closed: bool) -> BRep {
    if sections.is_empty() {
        return BRep::new();
    }

    let mut result = BRep::new();
    let mut all_face_refs = Vec::new();

    let n = sections.len();
    let segments = if closed { n } else { n - 1 };

    for seg_idx in 0..segments {
        let curr_section = &sections[seg_idx];
        let next_section = &sections[(seg_idx + 1) % n];

        let curr_polylines = extract_polylines_from_section(curr_section);
        let next_polylines = extract_polylines_from_section(next_section);

        for (curr_pts, next_pts) in curr_polylines.iter().zip(next_polylines.iter()) {
            if let Some(face_ref) = create_ruled_face_tshape(&mut result, curr_pts, next_pts) {
                all_face_refs.push(face_ref);
            }
        }
    }

    if !all_face_refs.is_empty() {
        let shell = result.add_tshell(all_face_refs);
        result.add_tsolid(vec![shell]);
    }

    result
}

fn extract_polylines_from_section(section: &SectionResult) -> Vec<Vec<DVec3>> {
    section.curves.iter().map(|curve| curve.curve.sample_points(33)).collect()
}

fn create_ruled_face_tshape(brep: &mut BRep, pts1: &[DVec3], pts2: &[DVec3]) -> Option<Shape> {
    let n = pts1.len().min(pts2.len());
    if n < 2 {
        return None;
    }

    let resampled1 = resample_polyline(pts1, n);
    let resampled2 = resample_polyline(pts2, n);

    let mut wire_edges = Vec::new();

    for i in 0..n - 1 {
        let v00 = brep.add_tvertex(resampled1[i]);
        let v01 = brep.add_tvertex(resampled1[i + 1]);
        let v10 = brep.add_tvertex(resampled2[i]);
        let v11 = brep.add_tvertex(resampled2[i + 1]);

        let e1 = brep.add_tedge(None, v00.clone(), v10.clone(), [0.0, 1.0]);
        let e2 = brep.add_tedge(None, v10.clone(), v01.clone(), [0.0, 1.0]);
        let e3 = brep.add_tedge(None, v01.clone(), v00, [0.0, 1.0]);

        let _e4 = brep.add_tedge(None, v01.clone(), v10.clone(), [0.0, 1.0]);
        let _e5 = brep.add_tedge(None, v10, v11.clone(), [0.0, 1.0]);
        let _e6 = brep.add_tedge(None, v11, v01, [0.0, 1.0]);

        wire_edges.push(e1);
        wire_edges.push(e2);
        wire_edges.push(e3);
    }

    let wire = brep.add_twire(wire_edges);
    Some(brep.add_tface(None, wire, vec![], None, None, vec![], true))
}

fn resample_polyline(pts: &[DVec3], n: usize) -> Vec<DVec3> {
    if pts.len() == n {
        return pts.to_vec();
    }
    if pts.len() < 2 || n < 2 {
        return pts.to_vec();
    }

    let mut lengths = vec![0.0];
    let mut total = 0.0;
    for i in 1..pts.len() {
        total += (pts[i] - pts[i - 1]).length();
        lengths.push(total);
    }

    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let target = total * i as f64 / (n - 1) as f64;
        let seg = lengths.windows(2)
            .position(|w| target >= w[0] && target <= w[1])
            .unwrap_or(lengths.len() - 2);
        let seg_start = lengths[seg];
        let seg_end = lengths[seg + 1];
        let seg_len = seg_end - seg_start;
        let t = if seg_len > TOLERANCE_LEN_MIN { (target - seg_start) / seg_len } else { 0.0 };
        result.push(pts[seg].lerp(pts[seg + 1], t));
    }

    result
}

// = =  Analytic section curves = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

/// One result curve from [`section_curves`].
#[derive(Debug, Clone)]
pub enum SectionCurve {
    /// Exact analytic curve returned when the face has a recognized analytic surface.
    Analytic(Curve3),
    /// Polyline fallback for parametric surfaces (BSpline, Bezier, Offset, Torus, ...).
    Polyline(Vec<DVec3>),
}

/// Section a BRep with a plane, returning analytic curves where possible.
pub fn section_curves(brep: &BRep, plane: &Plane) -> Vec<SectionCurve> {
    use crate::inttools::{
        plane_cone::{PlaneConicalResult, intersect_plane_cone},
        plane_cylinder::{PlaneCylinderResult, intersect_plane_cylinder},
        plane_plane::{PlanePlaneResult, intersect_plane_plane},
        plane_sphere::{PlaneSphereResult, intersect_plane_sphere},
    };

    let mut results: Vec<SectionCurve> = Vec::new();

    if brep.tshapes.is_empty() {
        return results;
    }

    let merge_eps = plane_section_mesh_merge_eps(brep);

    let face_list: Vec<Shape> = {
        let mut list = Vec::new();
        for ts in &brep.tshapes {
            if let TShape::Solid(sd) = ts.as_ref() {
                for shell_sr in &sd.shells {
                    if let TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
                        for face_sr in &shd.faces {
                            list.push(face_sr.clone());
                        }
                    }
                }
            }
        }
        list
    };

    let mut face_global_idx = 0usize;
    for face_ref in &face_list {
        let face: Shape = face_ref.clone();
        let surf_opt = face_surface_from_ref(brep, face);

        if let Some(surface) = surf_opt {
            let analytic = match surface {
                Surface3::Plane(face_plane) => {
                    match intersect_plane_plane(plane, face_plane) {
                        PlanePlaneResult::Line(line) => Some(Curve3::Line(line)),
                        _ => None,
                    }
                }
                Surface3::Sphere(sph) => match intersect_plane_sphere(plane, sph) {
                    PlaneSphereResult::Circle(c) => Some(Curve3::Circle(c)),
                    PlaneSphereResult::TangentPoint(_) => None,
                    PlaneSphereResult::NoIntersection => None,
                },
                Surface3::Cylinder(cyl) => match intersect_plane_cylinder(plane, cyl) {
                    PlaneCylinderResult::Circle(c) => Some(Curve3::Circle(c)),
                    PlaneCylinderResult::Ellipse(e) => Some(Curve3::Ellipse(e)),
                    PlaneCylinderResult::TwoLines(l1, _l2) => Some(Curve3::Line(l1)),
                    PlaneCylinderResult::TangentLine(_) => None,
                    PlaneCylinderResult::NoIntersection => None,
                },
                Surface3::Cone(cone) => match intersect_plane_cone(plane, cone) {
                    PlaneConicalResult::Circle(c) => Some(Curve3::Circle(c)),
                    PlaneConicalResult::Ellipse(e) => Some(Curve3::Ellipse(e)),
                    PlaneConicalResult::Parabola(par) => Some(Curve3::Parabola(par)),
                    PlaneConicalResult::Hyperbola(hyp) => Some(Curve3::Hyperbola(hyp)),
                    PlaneConicalResult::SingleLine(l) => Some(Curve3::Line(l)),
                    PlaneConicalResult::TwoLines(l1, _l2) => Some(Curve3::Line(l1)),
                    PlaneConicalResult::Point(_) => None,
                    PlaneConicalResult::NoIntersection => None,
                },
                _ => {
                    let tris = collect_face_triangles_tshape(brep, face_ref.clone());
                    let segs: Vec<[DVec3; 2]> = tris
                        .into_iter()
                        .filter_map(|tri| triangle_section(plane, tri))
                        .collect();
                    if !segs.is_empty() {
                        let chains = chain_segments_eps(segs, merge_eps);
                        for chain in chains {
                            if chain.len() >= 2 {
                                results.push(SectionCurve::Polyline(chain));
                            }
                        }
                    }
                    face_global_idx += 1;
                    continue;
                }
            };

            if let Some(curve) = analytic {
                results.push(SectionCurve::Analytic(curve));
            }
        } else {
            let tris = collect_face_triangles_tshape(brep, face_ref.clone());
            let segs: Vec<[DVec3; 2]> = tris
                .into_iter()
                .filter_map(|tri| triangle_section(plane, tri))
                .collect();
            if !segs.is_empty() {
                let chains = chain_segments_eps(segs, merge_eps);
                for chain in chains {
                    if chain.len() >= 2 {
                        results.push(SectionCurve::Polyline(chain));
                    }
                }
            }
        }

        face_global_idx += 1;
    }

    results
}

// = =  Tests = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
