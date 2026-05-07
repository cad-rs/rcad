use super::{
    BuildError, basis_from_axis_ref, basis_from_x_y, do_mirror_brep, normalize_vector,
    transform_brep, translate_brep, validate_point, validate_positive,
};
use glam::{DMat3, DVec3};
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::WireEdge;
use rcad_kernel::{BRep, PrimitiveSolid};

pub fn box_primitive(width: f64, height: f64, depth: f64) -> Result<PrimitiveSolid, BuildError> {
    let width = validate_positive("width", width)?;
    let height = validate_positive("height", height)?;
    let depth = validate_positive("depth", depth)?;
    Ok(PrimitiveSolid::Box {
        width,
        height,
        depth,
    })
}

pub fn make_box_primitive(
    width: f64,
    height: f64,
    depth: f64,
) -> Result<PrimitiveSolid, BuildError> {
    box_primitive(width, height, depth)
}

pub fn box_brep(
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    width: f64,
    height: f64,
    depth: f64,
) -> Result<BRep, BuildError> {
    let origin = validate_point("origin", origin)?;
    let primitive = box_primitive(width, height, depth)?;
    let (x_axis, y_axis, z_axis) = basis_from_x_y(x_dir, y_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, origin, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_box_brep(
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    width: f64,
    height: f64,
    depth: f64,
) -> Result<BRep, BuildError> {
    box_brep(origin, x_dir, y_dir, width, height, depth)
}

pub fn sphere_primitive(radius: f64) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    Ok(PrimitiveSolid::Sphere { radius })
}

pub fn make_sphere_primitive(radius: f64) -> Result<PrimitiveSolid, BuildError> {
    sphere_primitive(radius)
}

pub fn sphere_brep(center: DVec3, radius: f64) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = sphere_primitive(radius)?;
    let mut brep = BRep::from_primitive(primitive);
    translate_brep(&mut brep, center);
    Ok(brep)
}

pub fn make_sphere_brep(center: DVec3, radius: f64) -> Result<BRep, BuildError> {
    sphere_brep(center, radius)
}

pub fn cylinder_primitive(radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    let radius = validate_positive("radius", radius)?;
    let height = validate_positive("height", height)?;
    Ok(PrimitiveSolid::Cylinder { radius, height })
}

pub fn make_cylinder_primitive(radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    cylinder_primitive(radius, height)
}

pub fn cylinder_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = cylinder_primitive(radius, height)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_cylinder_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    radius: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    cylinder_brep(center, axis, ref_dir, radius, height)
}

pub fn cone_primitive(base_radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    let base_radius = validate_positive("base_radius", base_radius)?;
    let height = validate_positive("height", height)?;
    Ok(PrimitiveSolid::Cone {
        base_radius,
        height,
    })
}

pub fn make_cone_primitive(base_radius: f64, height: f64) -> Result<PrimitiveSolid, BuildError> {
    cone_primitive(base_radius, height)
}

pub fn cone_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    base_radius: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = cone_primitive(base_radius, height)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_cone_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    base_radius: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    cone_brep(center, axis, ref_dir, base_radius, height)
}

/// Right conical frustum (truncated cone), matching OCCT `pcone` when both bottom and top radii are positive.
///
/// `center` is the **midpoint** between the circular face centers; `axis` points from the bottom face
/// toward the top face; `r_bottom` / `r_top` are radii in those end planes; `height` is the distance
/// between the planes. Built with [`super::ops::loft`] between two regular polygon approximations.
pub fn make_conical_frustum_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    r_bottom: f64,
    r_top: f64,
    height: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let height = validate_positive("height", height)?;
    let rb = validate_positive("r_bottom", r_bottom)?;
    let rt = validate_positive("r_top", r_top)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let half_h = height * 0.5;
    let bottom_c = center - y_axis * half_h;
    let top_c = center + y_axis * half_h;
    const N: usize = 32;
    let mut lo = Vec::with_capacity(N);
    let mut hi = Vec::with_capacity(N);
    use std::f64::consts::TAU;
    for i in 0..N {
        let ang = TAU * (i as f64) / (N as f64);
        let c = ang.cos();
        let s = ang.sin();
        lo.push(bottom_c + x_axis * (c * rb) + z_axis * (s * rb));
        hi.push(top_c + x_axis * (c * rt) + z_axis * (s * rt));
    }
    super::ops::loft(&[lo, hi])
}

pub fn torus_primitive(major_radius: f64, minor_radius: f64) -> Result<PrimitiveSolid, BuildError> {
    let major_radius = validate_positive("major_radius", major_radius)?;
    let minor_radius = validate_positive("minor_radius", minor_radius)?;
    Ok(PrimitiveSolid::Torus {
        major_radius,
        minor_radius,
    })
}

pub fn make_torus_primitive(
    major_radius: f64,
    minor_radius: f64,
) -> Result<PrimitiveSolid, BuildError> {
    torus_primitive(major_radius, minor_radius)
}

pub fn torus_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<BRep, BuildError> {
    let center = validate_point("center", center)?;
    let primitive = torus_primitive(major_radius, minor_radius)?;
    let (x_axis, y_axis, z_axis) = basis_from_axis_ref(axis, ref_dir)?;
    let mut brep = BRep::from_primitive(primitive);
    transform_brep(&mut brep, center, x_axis, y_axis, z_axis);
    Ok(brep)
}

pub fn make_torus_brep(
    center: DVec3,
    axis: DVec3,
    ref_dir: DVec3,
    major_radius: f64,
    minor_radius: f64,
) -> Result<BRep, BuildError> {
    torus_brep(center, axis, ref_dir, major_radius, minor_radius)
}

/// Create a BRep with a single rectangular planar face.
///
/// `origin` is the plane origin; `u_axis`/`v_axis` are the parametric axes
/// (the normal is `u_axis × v_axis`). The corners are `origin + u_axis*u + v_axis*v`.
pub fn make_planar_rect_brep(
    origin: DVec3,
    u_axis: DVec3,
    v_axis: DVec3,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> Result<BRep, BuildError> {
    let normal = u_axis.cross(v_axis).normalize();
    let surface = Surface3::Plane(Plane { origin, normal });

    let c0 = origin + u_axis * umin + v_axis * vmin;
    let c1 = origin + u_axis * umax + v_axis * vmin;
    let c2 = origin + u_axis * umax + v_axis * vmax;
    let c3 = origin + u_axis * umin + v_axis * vmax;

    let mut brep = BRep::default();
    let v0 = super::brep_builder::make_vertex(&mut brep, c0);
    let v1 = super::brep_builder::make_vertex(&mut brep, c1);
    let v2 = super::brep_builder::make_vertex(&mut brep, c2);
    let v3 = super::brep_builder::make_vertex(&mut brep, c3);

    let len01 = (c1 - c0).length();
    let len12 = (c2 - c1).length();
    let len23 = (c3 - c2).length();
    let len30 = (c0 - c3).length();

    let dir01 = (c1 - c0).normalize();
    let dir12 = (c2 - c1).normalize();
    let dir23 = (c3 - c2).normalize();
    let dir30 = (c0 - c3).normalize();

    let e0 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c0, direction: dir01 }),
        0.0, len01, v0, v1,
    )?;
    let e1 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c1, direction: dir12 }),
        0.0, len12, v1, v2,
    )?;
    let e2 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c2, direction: dir23 }),
        0.0, len23, v2, v3,
    )?;
    let e3 = super::brep_builder::make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: c3, direction: dir30 }),
        0.0, len30, v3, v0,
    )?;

    let w = super::brep_builder::make_wire(vec![
        WireEdge::new(e0, true),
        WireEdge::new(e1, true),
        WireEdge::new(e2, true),
        WireEdge::new(e3, true),
    ]);

    super::brep_builder::make_face(&mut brep, surface, w, vec![])?;
    Ok(brep)
}

/// Create a BRep with a single planar face bounded by a polygon, with optional inner holes.
///
/// `origin` and `normal` define the plane. `outer_polygon` gives the outer boundary
/// vertices (≥3 points in CCW order). `inner_polygons` are hole boundaries (each ≥3 points
/// in CW order).
pub fn make_planar_polygon_brep(
    origin: DVec3,
    normal: DVec3,
    outer_polygon: &[DVec3],
    inner_polygons: &[Vec<DVec3>],
) -> Result<BRep, BuildError> {
    let surface = Surface3::Plane(Plane { origin, normal });
    let mut brep = BRep::default();
    let outer_wire = make_polygon_wire(&mut brep, outer_polygon)?;
    let mut inner_wires = Vec::new();
    for poly in inner_polygons {
        inner_wires.push(make_polygon_wire(&mut brep, poly)?);
    }
    super::brep_builder::make_face(&mut brep, surface, outer_wire, inner_wires)?;
    Ok(brep)
}

fn make_polygon_wire(brep: &mut BRep, points: &[DVec3]) -> Result<rcad_kernel::topology::Wire, BuildError> {
    let n = points.len();
    if n < 3 {
        return Err(BuildError::DegenerateGeometry("polygon needs at least 3 points"));
    }
    let mut vertex_indices = Vec::with_capacity(n);
    for p in points {
        vertex_indices.push(super::brep_builder::make_vertex(brep, *p));
    }
    let mut wire_edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let vi = vertex_indices[i];
        let vj = vertex_indices[j];
        let dir = (points[j] - points[i]).normalize();
        let len = (points[j] - points[i]).length();
        let ei = super::brep_builder::make_edge(
            brep,
            Curve3::Line(Line3 { origin: points[i], direction: dir }),
            0.0, len, vi, vj,
        )?;
        wire_edges.push(WireEdge::new(ei, true));
    }
    Ok(super::brep_builder::make_wire(wire_edges))
}

/// Create a half-space solid bounded by a plane: a large box extending from `origin` in the
/// `-normal` direction (the "interior" side). Used by OCCT `mkvolume` translation to intersect
/// planar half-spaces.
///
/// `size` controls the extents: ±size in the plane directions, size in the -normal direction.
pub fn make_half_space_brep(
    origin: DVec3,
    normal: DVec3,
    size: f64,
) -> Result<BRep, BuildError> {
    let normal = normalize_vector("normal", normal)?;
    let origin = validate_point("origin", origin)?;
    let size = validate_positive("size", size)?;

    // Choose u_axis perpendicular to normal
    let abs = normal.abs();
    let candidate = if abs.x <= abs.y && abs.x <= abs.z {
        DVec3::X
    } else if abs.y <= abs.z {
        DVec3::Y
    } else {
        DVec3::Z
    };
    let u_axis = normal.cross(candidate).normalize();
    // v_axis = u_axis × normal gives u_axis × v_axis = -normal (the interior direction).
    let v_axis = u_axis.cross(normal).normalize();

    // box_brep creates a corner-box at `origin - u_axis*size - v_axis*size` that extends
    // +u_axis*2*size, +v_axis*2*size, +z_axis*size, giving ±size UV coverage around `origin`
    // and `size` units in the -normal (interior) direction.
    let corner = origin - u_axis * size - v_axis * size;
    box_brep(corner, u_axis, v_axis, 2.0 * size, 2.0 * size, size)
}

/// Build a convex polyhedron BRep directly from half-space plane equations (no boolean ops).
///
/// Each half-space is `(origin, normal)` where the interior is the side opposite the normal,
/// i.e. `normal · (p - origin) ≤ 0`. Input planes may contain coplanar duplicates — they
/// are merged into a single face by keeping the tightest constraint.
///
/// The algorithm computes all vertices by intersecting 3 planes at a time, keeps those
// satisfying all half-spaces, then builds faces from coplanar vertex sets.
pub fn make_convex_polyhedron_from_half_spaces(
    planes: &[(DVec3, DVec3)],
) -> Result<BRep, BuildError> {
    if planes.len() < 4 {
        return Err(BuildError::DegenerateGeometry(
            "need at least 4 planes for a convex polyhedron",
        ));
    }

    // Geometry-adaptive tolerance
    let extent: f64 = planes
        .iter()
        .map(|(o, _)| o.length_squared())
        .fold(0.0, f64::max)
        .sqrt()
        .max(1.0);
    let tol = extent * 1e-12 + 1e-10;

    // --- 1. Deduplicate coplanar planes, keep tightest constraint ---
    struct PlaneEq {
        n: DVec3,
        d: f64,
        origin: DVec3,
    }
    let mut eqs: Vec<PlaneEq> = Vec::new();

    for &(origin, normal) in planes {
        let n = normal.normalize();
        let d = n.dot(origin);
        let mut found = false;
        for eq in &mut eqs {
            let cos = n.dot(eq.n);
            if cos > 1.0 - tol * tol * tol {
                // Same direction – keep the tighter constraint (smaller d for n·p ≤ d)
                if d < eq.d {
                    eq.d = d;
                    eq.origin = origin;
                }
                found = true;
                break;
            }
            if cos < -1.0 + tol * tol * tol {
                // Opposite direction – one plane bounds +n, the other -n.
                // Both constraints are independent; keep both.
                // (But skip if one strictly dominates the other.)
                found = true;
                break;
            }
        }
        if !found {
            eqs.push(PlaneEq { n, d, origin });
        }
    }

    if eqs.len() < 4 {
        return Err(BuildError::DegenerateGeometry(
            "fewer than 4 unique planes after deduplication",
        ));
    }

    // --- 2. Compute all intersection points of 3 planes ---
    let mut verts: Vec<DVec3> = Vec::new();
    let np = eqs.len();

    for i in 0..np {
        for j in (i + 1)..np {
            for k in (j + 1)..np {
                let (ni, nj, nk) = (eqs[i].n, eqs[j].n, eqs[k].n);
                // Rows are plane normals (ni·p=di, …); columns store x/y/z components per plane.
                let m = DMat3::from_cols(
                    DVec3::new(ni.x, nj.x, nk.x),
                    DVec3::new(ni.y, nj.y, nk.y),
                    DVec3::new(ni.z, nj.z, nk.z),
                );
                if m.determinant().abs() < tol {
                    continue; // singular — planes meet in a line or are parallel
                }
                let p = m.inverse() * DVec3::new(eqs[i].d, eqs[j].d, eqs[k].d);

                // Check against all half-space constraints
                let valid = eqs.iter().all(|eq| eq.n.dot(p) - eq.d <= tol * tol);
                if !valid {
                    continue;
                }

                // Deduplicate
                let dup = verts.iter().any(|v| (v - p).length_squared() < tol * tol);
                if !dup {
                    verts.push(p);
                }
            }
        }
    }

    if verts.len() < 4 {
        return Err(BuildError::DegenerateGeometry(
            "fewer than 4 valid vertices — half-spaces may not form a bounded convex polyhedron",
        ));
    }

    // --- 3. Group vertices by plane ---
    let mut face_verts: Vec<Vec<usize>> = vec![Vec::new(); np];
    for (vi, &v) in verts.iter().enumerate() {
        for (ei, eq) in eqs.iter().enumerate() {
            if (eq.n.dot(v) - eq.d).abs() < tol * tol * tol {
                face_verts[ei].push(vi);
            }
        }
    }

    // --- 4. Build BRep ---
    use super::brep_builder::*;
    let mut brep = BRep::default();

    // Create all vertices upfront so edges on adjacent faces share indices
    let bv: Vec<usize> = verts.iter().map(|&v| make_vertex(&mut brep, v)).collect();

    for fi in 0..np {
        if face_verts[fi].len() < 3 {
            continue;
        }
        let n = eqs[fi].n;
        let origin = eqs[fi].origin;

        // Centroid of face vertices
        let centroid: DVec3 = face_verts[fi]
            .iter()
            .map(|&vi| verts[vi])
            .sum::<DVec3>()
            / face_verts[fi].len() as f64;

        // UV basis in the plane (u × v = n for CCW winding from +n)
        let abs = n.abs();
        let candidate = if abs.x <= abs.y && abs.x <= abs.z {
            DVec3::X
        } else if abs.y <= abs.z {
            DVec3::Y
        } else {
            DVec3::Z
        };
        let u_axis = n.cross(candidate).normalize();
        let v_axis = n.cross(u_axis).normalize(); // u × v = n

        // Sort vertices by angle around centroid in the UV plane
        let mut sorted: Vec<usize> = face_verts[fi].clone();
        sorted.sort_by(|&a, &b| {
            let da = verts[a] - centroid;
            let db = verts[b] - centroid;
            let aa = da.dot(u_axis).atan2(da.dot(v_axis));
            let ab = db.dot(u_axis).atan2(db.dot(v_axis));
            aa.partial_cmp(&ab).unwrap()
        });

        // Build wire
        let mut wire_edges = Vec::with_capacity(sorted.len());
        for w in 0..sorted.len() {
            let vi0 = bv[sorted[w]];
            let vi1 = bv[sorted[(w + 1) % sorted.len()]];
            let p0 = verts[sorted[w]];
            let p1 = verts[sorted[(w + 1) % sorted.len()]];
            let dir = (p1 - p0).normalize();
            let len = (p1 - p0).length();
            let e = make_edge(
                &mut brep,
                Curve3::Line(Line3 {
                    origin: p0,
                    direction: dir,
                }),
                0.0,
                len,
                vi0,
                vi1,
            )?;
            wire_edges.push(WireEdge::new(e, true));
        }

        let wire = make_wire(wire_edges);
        let surface = Surface3::Plane(Plane { origin, normal: n });
        make_face(&mut brep, surface, wire, vec![])?;
    }

    Ok(brep)
}

/// Mirror a BRep across a plane defined by `origin` and `normal`.
///
/// The mirrored BRep has inverted face normals and reversed wire orientations
/// to maintain consistent outward-facing normals.
pub fn mirror_brep(brep: &BRep, plane_origin: DVec3, plane_normal: DVec3) -> Result<BRep, BuildError> {
    let _ = validate_point("plane_origin", plane_origin)?;
    let n = normalize_vector("plane_normal", plane_normal)?;
    Ok(do_mirror_brep(brep, plane_origin, n))
}
