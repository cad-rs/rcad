
/// Build a tessellated BRep for `cone 閳?box` via Z-slice tessellation.
///
/// The cone is a Z-aligned conical frustum with center at `(cx, cy)` in XY,
/// extending from Z `cz_lo` to `cz_hi`, with bottom radius `cr_lo` and top
/// radius `cr_hi`. The box is axis-aligned `[bmin, bmax]`.
///
/// Builds only the overlap section (circle 閳?rect).
fn build_cone_box_intersection_tessellated(
    bmin: DVec3, bmax: DVec3,
    cx: f64, cy: f64,
    cz_lo: f64, cz_hi: f64,
    cr_lo: f64, cr_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cz_hi <= cz_lo + tol { return None; }
    if cr_lo < tol && cr_hi < tol { return None; }

    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let box_center = (bmin + bmax) * 0.5;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let cu = cx - box_center.x;
    let cv = cy - box_center.y;

    let n_slices = 64usize;
    let n_boundary = 256usize;
    let empty_wire = || Wire { edges: vec![] };

    // Overlap Z range
    let z0 = cz_lo.max(box_z_lo);
    let z1 = cz_hi.min(box_z_hi);
    if z1 <= z0 + tol { return None; }

    let mut verts: Vec<Vertex> = Vec::new();
    let mut add_v = |p: DVec3| -> usize {
        for (i, v) in verts.iter().enumerate() {
            if (v.point - p).length() < 1e-12 { return i; }
        }
        let idx = verts.len();
        verts.push(Vertex { point: p });
        idx
    };

    let mut faces: Vec<Face> = Vec::new();

    let to_world = |u: f64, v: f64, z: f64| -> DVec3 {
        DVec3::new(box_center.x + u, box_center.y + v, z)
    };

    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let ref_pt = DVec2::new(-eu, -ev);

    // Build wall faces via Z-slice tessellation
    let n = n_slices;
    let dz = (z1 - z0) / n as f64;

    // Pre-compute and remap all slice polygons (circle 閳?rect)
    let mut slices: Vec<Vec<DVec2>> = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let z = z0 + dz * i as f64;
        let r = (cr_lo + dr_dz * (z - cz_lo)).max(0.0);
        if r < tol {
            slices.push(vec![]);
        } else {
            let poly = build_circle_intersect_rect_polygon(cu, cv, r, eu, ev);
            if poly.len() >= 3 {
                slices.push(remap_polygon_arclength(&poly, n_boundary, ref_pt));
            } else {
                slices.push(vec![]);
            }
        }
    }

    // Build wall faces between adjacent remapped slices
    for i in 0..n {
        let bot = &slices[i];
        let top = &slices[i + 1];
        if bot.len() < 3 || top.len() < 3 { continue; }

        let n_pts = bot.len().min(top.len());
        let z_bot = z0 + dz * i as f64;
        let z_top = z0 + dz * (i + 1) as f64;

        let mut idx = Vec::with_capacity(2 * n_pts);
        for p in bot.iter() { idx.push(add_v(to_world(p.x, p.y, z_bot))); }
        for p in top.iter() { idx.push(add_v(to_world(p.x, p.y, z_top))); }

        let mut tris = Vec::with_capacity(n_pts * 2);
        for j in 0..n_pts {
            let k = (j + 1) % n_pts;
            tris.push([idx[j], idx[k], idx[n_pts + k]]);
            tris.push([idx[j], idx[n_pts + k], idx[n_pts + j]]);
        }
        faces.push(Face {
            outer_wire: empty_wire(), inner_wires: vec![],
            normal: DVec3::ZERO, triangles: tris,
            sample_point: None, mesh_dirty: false,
                surface_idx: None,
        });
    }

    // ---- Bottom cap ----
    let r_bot = (cr_lo + dr_dz * (z0 - cz_lo)).max(0.0);
    if r_bot > tol {
        let poly = build_circle_intersect_rect_polygon(cu, cv, r_bot, eu, ev);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, z0, -DVec3::Z, &to_world, &empty_wire);
        }
    }

    // ---- Top cap ----
    let r_top = (cr_lo + dr_dz * (z1 - cz_lo)).max(0.0);
    if r_top > tol {
        let poly = build_circle_intersect_rect_polygon(cu, cv, r_top, eu, ev);
        if poly.len() >= 3 {
            add_cap_face(&mut add_v, &mut faces, &poly, z1, DVec3::Z, &to_world, &empty_wire);
        }
    }

    if faces.is_empty() { return None; }

    let geom = GeomStore { face_internal_vertices: vec![],
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![],
    };

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Detect one direction of cone-box union.
fn try_union_cone_box_one_dir(cone_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(cone_brep)?;
    let [bmin, bmax] = try_as_axis_aligned_box(box_brep)?;

    let cx = center_xy.x;
    let cy = center_xy.y;
    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let tol = TOLERANCE_LEN_MIN;

    // Check Z overlap
    let inter_lo = cz_lo.max(box_z_lo);
    let inter_hi = cz_hi.min(box_z_hi);
    if inter_hi <= inter_lo + tol { return None; }

    // Compute radii at overlap boundaries
    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_inter_lo = (cr_lo + dr_dz * (inter_lo - cz_lo)).max(0.0);
    let r_at_inter_hi = (cr_lo + dr_dz * (inter_hi - cz_lo)).max(0.0);

    // Quick check: cone XY reach at overlap Z must reach the box XY
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    let min_r = r_at_inter_lo.min(r_at_inter_hi);
    if dist_center > box_half_diag + min_r + tol {
        return None; // Disjoint in XY
    }

    // Check containment: box entirely inside cone 閳?union is the cone
    let corners = [
        DVec2::new(bmin.x, bmin.y),
        DVec2::new(bmax.x, bmin.y),
        DVec2::new(bmax.x, bmax.y),
        DVec2::new(bmin.x, bmax.y),
    ];
    let all_corners_inside_at_lo = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_lo + tol).powi(2)
    });
    let all_corners_inside_at_hi = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_hi + tol).powi(2)
    });
    let box_z_inside_cone = box_z_lo >= cz_lo - tol && box_z_hi <= cz_hi + tol;
    if all_corners_inside_at_lo && all_corners_inside_at_hi && box_z_inside_cone {
        return Some(cone_brep.clone());
    }

    // Check containment: cone entirely inside box 閳?union is the box
    let cone_inside_box_xy_at_lo = cx - r_at_inter_lo >= bmin.x - tol
        && cx + r_at_inter_lo <= bmax.x + tol
        && cy - r_at_inter_lo >= bmin.y - tol
        && cy + r_at_inter_lo <= bmax.y + tol;
    let cone_inside_box_xy_at_hi = cx - r_at_inter_hi >= bmin.x - tol
        && cx + r_at_inter_hi <= bmax.x + tol
        && cy - r_at_inter_hi >= bmin.y - tol
        && cy + r_at_inter_hi <= bmax.y + tol;
    let cone_z_inside_box = cz_lo >= box_z_lo - tol && cz_hi <= box_z_hi + tol;
    if cone_inside_box_xy_at_lo && cone_inside_box_xy_at_hi && cone_z_inside_box {
        return Some(box_brep.clone());
    }

    // Need the tessellated path: verify the cross-section is non-empty
    let test_r = r_at_inter_lo.max(r_at_inter_hi);
    let test_poly = build_circle_union_rect_polygon(
        cx - box_center_xy.x, cy - box_center_xy.y,
        test_r, eu, ev,
    );
    if test_poly.len() < 3 { return None; }

    // Try analytic builder first
    if let Some(result) = crate::cone_box_analytic::build_cone_box_union_analytic(cone_brep, box_brep) {
        return Some(result);
    }

    build_cone_box_union_tessellated(bmin, bmax, cx, cy, cz_lo, cz_hi, cr_lo, cr_hi)
}

/// Fast path: cone-box Union via Z-slice tessellation.
pub fn try_union_cone_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_cone_box_one_dir(a, b).or_else(|| try_union_cone_box_one_dir(b, a))
}

/// Detect one direction of cone-box intersection.
fn try_intersection_cone_box_one_dir(cone_brep: &BRep, box_brep: &BRep) -> Option<BRep> {
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(cone_brep)?;
    let [bmin, bmax] = try_as_axis_aligned_box(box_brep)?;

    let cx = center_xy.x;
    let cy = center_xy.y;
    let box_z_lo = bmin.z;
    let box_z_hi = bmax.z;
    let eu = (bmax.x - bmin.x) * 0.5;
    let ev = (bmax.y - bmin.y) * 0.5;
    let tol = TOLERANCE_LEN_MIN;

    // Check Z overlap
    let inter_lo = cz_lo.max(box_z_lo);
    let inter_hi = cz_hi.min(box_z_hi);
    if inter_hi <= inter_lo + tol { return None; }

    // Compute radii at overlap boundaries
    let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_inter_lo = (cr_lo + dr_dz * (inter_lo - cz_lo)).max(0.0);
    let r_at_inter_hi = (cr_lo + dr_dz * (inter_hi - cz_lo)).max(0.0);

    // Quick check: cone XY reach at overlap Z must reach the box XY
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    let min_r = r_at_inter_lo.min(r_at_inter_hi);
    if dist_center > box_half_diag + min_r + tol {
        return None; // Disjoint in XY
    }

    // Check containment: box entirely inside cone 閳?intersection is the box
    let corners = [
        DVec2::new(bmin.x, bmin.y),
        DVec2::new(bmax.x, bmin.y),
        DVec2::new(bmax.x, bmax.y),
        DVec2::new(bmin.x, bmax.y),
    ];
    let all_corners_inside_at_lo = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_lo + tol).powi(2)
    });
    let all_corners_inside_at_hi = corners.iter().all(|c| {
        (c.x - cx).powi(2) + (c.y - cy).powi(2) <= (r_at_inter_hi + tol).powi(2)
    });
    let box_z_inside_cone = box_z_lo >= cz_lo - tol && box_z_hi <= cz_hi + tol;
    if all_corners_inside_at_lo && all_corners_inside_at_hi && box_z_inside_cone {
        return Some(box_brep.clone());
    }

    // Check containment: cone entirely inside box 閳?intersection is the cone
    let cone_inside_box_xy_at_lo = cx - r_at_inter_lo >= bmin.x - tol
        && cx + r_at_inter_lo <= bmax.x + tol
        && cy - r_at_inter_lo >= bmin.y - tol
        && cy + r_at_inter_lo <= bmax.y + tol;
    let cone_inside_box_xy_at_hi = cx - r_at_inter_hi >= bmin.x - tol
        && cx + r_at_inter_hi <= bmax.x + tol
        && cy - r_at_inter_hi >= bmin.y - tol
        && cy + r_at_inter_hi <= bmax.y + tol;
    let cone_z_inside_box = cz_lo >= box_z_lo - tol && cz_hi <= box_z_hi + tol;
    if cone_inside_box_xy_at_lo && cone_inside_box_xy_at_hi && cone_z_inside_box {
        return Some(cone_brep.clone());
    }

    // Need the tessellated path: verify the cross-section is non-empty
    let test_r = r_at_inter_lo.max(r_at_inter_hi);
    let test_poly = build_circle_intersect_rect_polygon(
        cx - box_center_xy.x, cy - box_center_xy.y,
        test_r, eu, ev,
    );
    if test_poly.len() < 3 { return None; }

    // Try analytic builder first
    if let Some(result) = crate::cone_box_analytic::build_cone_box_intersection_analytic(cone_brep, box_brep) {
        return Some(result);
    }

    build_cone_box_intersection_tessellated(bmin, bmax, cx, cy, cz_lo, cz_hi, cr_lo, cr_hi)
}

/// Fast path: cone-box Intersection via Z-slice tessellation.
pub fn try_intersection_cone_box(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_cone_box_one_dir(a, b).or_else(|| try_intersection_cone_box_one_dir(b, a))
}

/// Fast path: union of a box-with-cavity solid and its matching cavity shell.
///
/// When `boolean_op(Difference, big_box, small_box)` produces a hollow box, the
/// Pave-Filler's fuse step can't merge the cavity walls back into the solid.
/// Detect this case and return the outer box directly.
///
/// Conditions:
/// - `a` is a solid whose outer shell is a box (6 planar faces)
/// - `a` has at least one inner shell (a cavity)
/// - `b` is a simple box (passes `try_as_box`)
/// - `b` fits entirely inside `a`'s bounding box
pub fn try_union_fill_box_cavity(a: &BRep, b: &BRep) -> Option<BRep> {
    if a.solids.len() != 1 {
        return None;
    }
    let a_solid = &a.solids[0];
    if a_solid.shells.is_empty() {
        return None;
    }

    // Count total faces across all shells.
    let total_faces: usize = a_solid.shells.iter().map(|sh| sh.faces.len()).sum();
    if total_faces < 6 {
        return None;
    }

    // Must have evidence of a cavity: either multiple shells, or >6 faces in total
    // (the latter occurs when extract_solids collapsed multiple shells into one).
    if a_solid.shells.len() < 2 && total_faces <= 6 {
        return None;
    }

    // The outer shell's faces must all be planar (forms a box boundary).
    let outer_face_count = a_solid.shells[0].faces.len();
    for fi in 0..outer_face_count {
        let si_slot = a.geom.face_surface.get(fi);
        let si = *si_slot?.as_ref()?;
        match a.geom.surfaces.get(si)? {
            Surface3::Plane(_) => {}
            _ => return None,
        }
    }

    // B must be inside A's bounding box (filling the cavity).
    // If B is empty (no faces), skip this check 閳?A has internal faces that need
    // cleaning (e.g., box-box difference with collapsed shells).
    let [amin, amax] = a.bounding_box()?;
    let tol = 1e-6;
    if crate::brep_algo::total_surface_area(b) > 0.0 {
        let [bmin, bmax] = b.bounding_box()?;
        if bmin.x < amin.x - tol
            || bmin.y < amin.y - tol
            || bmin.z < amin.z - tol
            || bmax.x > amax.x + tol
            || bmax.y > amax.y + tol
            || bmax.z > amax.z + tol
        {
            return None;
        }
    }

    let w = amax.x - amin.x;
    let h = amax.y - amin.y;
    let d = amax.z - amin.z;
    if w <= tol || h <= tol || d <= tol {
        return None;
    }

    make_box_brep(amin, DVec3::X, DVec3::Y, w, h, d).ok()
}

/// Fast path for `cone 閳?box` boolean difference.
///
/// Detects a Z-aligned conical frustum (possibly Z-rotated and translated in XY)
/// minus an axis-aligned box.  Builds the result via Z-slice tessellation (the
/// inverse of [`try_difference_box_cone`]: the kept shape is the cone with a
/// box-shaped channel removed).
pub fn try_difference_cone_box(a: &BRep, b: &BRep) -> Option<BRep> {
    // Detect Z-aligned cone frustum (a)
    let (center_xy, cz_lo, cz_hi, cr_lo, cr_hi) = detect_z_axis_cone_frustum(a)?;

    // Detect axis-aligned box (b)
    let [bmin, bmax] = try_as_axis_aligned_box(b)?;

    // Compute Z overlap
    let z_lo = cz_lo.max(bmin.z);
    let z_hi = cz_hi.min(bmax.z);
    if z_hi <= z_lo + TOLERANCE_LEN_MIN {
        return Some(a.clone());
    }

    let cx = center_xy.x;
    let cy = center_xy.y;

    // Clamp radii to non-negative
    let dr = (cr_hi - cr_lo) / (cz_hi - cz_lo);
    let r_at_zlo = cr_lo + dr * (z_lo - cz_lo);
    let r_at_zhi = cr_lo + dr * (z_hi - cz_lo);
    let r_lo = r_at_zlo.max(TOLERANCE_COORD_SUB);
    let r_hi = r_at_zhi.max(TOLERANCE_COORD_SUB);

    // Quick check: if the cone doesn't reach the box XY at the overlap Z,
    // return cone unchanged.
    let min_r = r_lo.min(r_hi);
    let box_half_diag = ((bmax.x - bmin.x).powi(2) + (bmax.y - bmin.y).powi(2)).sqrt() * 0.5;
    let box_center_xy = DVec2::new((bmin.x + bmax.x) * 0.5, (bmin.y + bmax.y) * 0.5);
    let dist_center = (box_center_xy - DVec2::new(cx, cy)).length();
    if dist_center > box_half_diag + min_r + TOLERANCE_LEN_MIN {
        // Box is entirely outside the cone's XY reach 閳?no intersection
        return Some(a.clone());
    }

    // Phase 3: Check if cone cross-section circle is fully inside the box XY rect
    // at every Z level in the overlap range. Since radius varies linearly with Z and
    // containment is convex, checking at both Z boundaries is sufficient.
    if cx - r_lo >= bmin.x - TOLERANCE_LEN_MIN
        && cx + r_lo <= bmax.x + TOLERANCE_LEN_MIN
        && cy - r_lo >= bmin.y - TOLERANCE_LEN_MIN
        && cy + r_lo <= bmax.y + TOLERANCE_LEN_MIN
        && cx - r_hi >= bmin.x - TOLERANCE_LEN_MIN
        && cx + r_hi <= bmax.x + TOLERANCE_LEN_MIN
        && cy - r_hi >= bmin.y - TOLERANCE_LEN_MIN
        && cy + r_hi <= bmax.y + TOLERANCE_LEN_MIN
    {
        use rcad_modeling::make_conical_frustum_brep;
        let dr_dz = (cr_hi - cr_lo) / (cz_hi - cz_lo);
        let mut parts: Vec<BRep> = Vec::new();

        // Cone portion below the box Z range
        if cz_lo < bmin.z - TOLERANCE_LEN_MIN {
            let z_to = bmin.z.min(cz_hi);
            let h = z_to - cz_lo;
            if h > TOLERANCE_LEN_MIN {
                let r_at_z_to = cr_lo + dr_dz * (z_to - cz_lo);
                if let Ok(f) = make_conical_frustum_brep(
                    DVec3::new(cx, cy, (cz_lo + z_to) * 0.5),
                    DVec3::Z,
                    DVec3::X,
                    cr_lo,
                    r_at_z_to,
                    h,
                ) {
                    parts.push(f);
                }
            }
        }

        // Cone portion above the box Z range
        if cz_hi > bmax.z + TOLERANCE_LEN_MIN {
            let z_from = bmax.z.max(cz_lo);
            let h = cz_hi - z_from;
            if h > TOLERANCE_LEN_MIN {
                let r_at_z_from = cr_lo + dr_dz * (z_from - cz_lo);
                if let Ok(f) = make_conical_frustum_brep(
                    DVec3::new(cx, cy, (z_from + cz_hi) * 0.5),
                    DVec3::Z,
                    DVec3::X,
                    r_at_z_from,
                    cr_hi,
                    h,
                ) {
                    parts.push(f);
                }
            }
        }

        if parts.is_empty() {
            // Cone is entirely inside the box 閳?empty result
            return Some(BRep::default());
        }
        if parts.len() == 1 {
            return Some(parts.swap_remove(0));
        }
        // Multiple portions (both below and above) 閳?fall through to tessellated path
        // which handles the compound case.
    }

    // Phase 4: analytic fast path (extracted to cone_box_analytic module)
    if let Some(result) = crate::cone_box_analytic::build_cone_minus_box_analytic(a, b) {
        return Some(result);
    }

    build_cone_minus_box_tessellated(bmin, bmax, cx, cy, z_lo, z_hi, r_lo, r_hi)
}

/// Fast path: cylinder 閳?box boolean difference.
///
/// When the box fully contains the cylinder cross-section in XY (at the Z-range
/// of intersection), the result is simply the cylinder clipped to Z ranges not
/// covered by the box 閳?avoiding the Pave-Filler which can produce hundreds
/// of unnecessary faces for symmetric configurations (e.g. bopcut_simple V5/V6).
pub fn try_difference_box_minus_brep_with_hole(a: &BRep, b: &BRep) -> Option<BRep> {
    let _bx = try_as_box(a)?;

    // b must have at least one inner wire (hole) on a planar face.
    let shell = b.solids.first()?.shells.first()?;
    let has_hole = shell.faces.iter().any(|f| !f.inner_wires.is_empty());
    if !has_hole {
        return None;
    }

    // Find the cylindrical surface in b.
    let cyl_si = b.geom.surfaces.iter().position(|s| matches!(s, Surface3::Cylinder(_)))?;
    let Surface3::Cylinder(cyl) = &b.geom.surfaces[cyl_si] else { return None; };

    // Find the face that uses this cylindrical surface 閳?get V range for height.
    let cyl_fi = b.geom.face_surface.iter().position(|fs| fs.map_or(false, |si| si == cyl_si))?;
    let v_range = b.geom.face_surface_range.get(cyl_fi)
        .and_then(|r| *r)
        .map(|r| (r[2], r[3]))?;
    let height = v_range.1 - v_range.0;
    if height <= 1e-12 {
        return None;
    }

    // Build a cylinder BRep matching the extracted parameters.
    let cyl_center = cyl.origin + cyl.axis * (height / 2.0);
    let cyl_brep = make_cylinder_brep(cyl_center, cyl.axis, cyl.ref_dir, cyl.radius, height).ok()?;

    // Compute box 閳?cylinder via the existing intersection fast path.
    // try_intersection_cylinder_box tries both (a, cylinder) and (cylinder, a).
    try_intersection_cylinder_box(a, &cyl_brep)
}

/// When the cylinder exits through exactly one box face, build the difference
/// (box 閳?cylinder) analytically: clone the box, remove the exit face, and add
/// a half-cylinder wall face to seal the opening.
///
/// The cylinder axis must be Z-aligned (enforced by the caller).  Currently
/// only handles the case where the cylinder Z-range matches the box Z-range
/// (z_lo 閳?閳姀w, z_hi 閳?ew) 閳?otherwise returns `None`.
fn build_box_minus_cylinder_one_v_face_exit(
    box_brep: &BRep,
    bx: &BoxInfo,
    _u_idx: usize, _v_idx: usize, z_idx: usize,
    _cu: f64, _cv: f64, _cz: f64,
    eu: f64, ev: f64, ew: f64,
    bc: DVec3, u_ax: DVec3, v_ax: DVec3,
    cyl_origin: DVec3,
    cyl_r: f64, cyl_height: f64,
    exit_side: i32,
    is_u: bool,
) -> Option<BRep> {
    // Cylinder must span the full box Z-range so the wall connects directly
    // to the Z- and Z+ box faces (no cap faces needed).
    let z_ax = bx.axes[z_idx];
    let cyl_center = cyl_origin + z_ax * (cyl_height / 2.0);
    let cz = (cyl_center - bc).dot(z_ax);
    let z_lo = cz - cyl_height / 2.0;
    let z_hi = cz + cyl_height / 2.0;
    if (z_lo - (-ew)).abs() > 1e-6 || (z_hi - ew).abs() > 1e-6 {
        return None;
    }

    // 1. Clone the box and remove the exit face.
    let mut brep = box_brep.clone();
    let exit_normal = if is_u {
        exit_side as f64 * u_ax
    } else {
        exit_side as f64 * v_ax
    };
    // Box faces use exact world-space normals 閳?no tolerance needed.
    let exit_fi = brep.solids[0].shells[0].faces.iter().position(|f| f.normal == exit_normal)?;
    brep.solids[0].shells[0].faces.remove(exit_fi);

    // 2. Identify the 4 vertices on the exit face (they form the wall boundary).
    // Collect unique vertex indices from the exit face's outer wire edges.
    // Box face normalised axes: u_ax=X, v_ax=Y, z_ax=Z 閳?standard vertex layout.
    let cyl_cx = cyl_origin.x; // centre X of cylinder (surface origin)
    let cyl_cy = cyl_origin.y; // centre Y of cylinder (surface origin)

    // Vertices of the exit face in wire-traversal order.
    // For the standard box (u=X, v=Y, z=Z):
    //   Y- face (F2): v0(x_min,y_min,z_min), v1(x_max,y_min,z_min),
    //                 v5(x_max,y_min,z_max), v4(x_min,y_min,z_max)
    //   Y+ face (F3): v3(x_min,y_max,z_min), v2(x_max,y_max,z_min),
    //                 v6(x_max,y_max,z_max), v7(x_min,y_max,z_max)
    //   X- face (F4): v0(x_min,y_min,z_min), v3(x_min,y_max,z_min),
    //                 v7(x_min,y_max,z_max), v4(x_min,y_min,z_max)
    //   X+ face (F5): v1(x_max,y_min,z_min), v2(x_max,y_max,z_min),
    //                 v6(x_max,y_max,z_max), v5(x_max,y_min,z_max)
    //
    // We map the 4 vertices to cylinder-wall corners:
    //   va = "right/down"  @ z_min (鑳?閳?0 mod 2锜?
    //   vb = "left/up"     @ z_min (鑳?閳?锜?
    //   vc = "left/up"     @ z_max (鑳?閳?锜?
    //   vd = "right/down"  @ z_max (鑳?閳?0 mod 2锜?
    // Wire order: va 閳?vb (bottom arc), vb 閳?vc (left gen),
    //             vc 閳?vd (top arc),   vd 閳?va (right gen reversed)

    fn sign_f64(v: f64) -> f64 { if v < 0.0 { -1.0 } else { 1.0 } }

    let (va, vb, vc, vd) = if !is_u {
        // Exit through v-direction face.
        let v_exit = exit_side as f64 * ev;
        let z = |zi: f64, off: f64| bc + zi * eu * u_ax + v_exit * v_ax + off * ew * z_ax;
        // Find vertex indices by matching positions (box vertices are exact).
        let va_pos = z(1.0, -1.0); // (+eu, exit_v, 閳姀w) 閳?right-bottom, 鑳?0
        let vb_pos = z(-1.0, -1.0); // (閳姀u, exit_v, 閳姀w) 閳?left-bottom, 鑳?锜?
        let vc_pos = z(-1.0, 1.0); // (閳姀u, exit_v, +ew) 閳?left-top,  鑳?锜?
        let vd_pos = z(1.0, 1.0); // (+eu, exit_v, +ew) 閳?right-top, 鑳?0
        (va_pos, vb_pos, vc_pos, vd_pos)
    } else {
        // Exit through u-direction face.
        let u_exit = exit_side as f64 * eu;
        let z = |vi: f64, off: f64| bc + u_exit * u_ax + vi * ev * v_ax + off * ew * z_ax;
        let va_pos = z(-1.0, -1.0); // (exit_u, 閳姀v, 閳姀w) 閳?bottom at v_min, 鑳?0
        let vb_pos = z(1.0, -1.0); // (exit_u, +ev, 閳姀w) 閳?bottom at v_max, 鑳?锜?
        let vc_pos = z(1.0, 1.0); // (exit_u, +ev, +ew) 閳?top at v_max,  鑳?锜?
        let vd_pos = z(-1.0, 1.0); // (exit_u, 閳姀v, +ew) 閳?top at v_min,  鑳?0
        (va_pos, vb_pos, vc_pos, vd_pos)
    };

    // Find box vertex indices by position.
    let tol_vi = 1e-10;
    let find_vi = |p: DVec3| -> Option<usize> {
        brep.vertices.iter().position(|v| v.point.distance_squared(p) < tol_vi)
    };
    let vi_va = find_vi(va)?;
    let vi_vb = find_vi(vb)?;
    let vi_vc = find_vi(vc)?;
    let vi_vd = find_vi(vd)?;

    // 3. Create edges for the cylinder wall.
    let h = 2.0 * ew; // wall height = box Z-extent

    // Bottom arc at z = 閳姀w: from va 閳?vb along y閳? (鑳? 0閳崻鈧?on Circle3 normal=+Z)
    let pi = std::f64::consts::PI;
    let bot_center = DVec3::new(cyl_cx, cyl_cy, -ew);
    let bot_arc = make_edge(
        &mut brep,
        Curve3::Circle(Circle3::new(bot_center, DVec3::Z, cyl_r )),
        0.0, pi, vi_va, vi_vb,
    ).ok()?;

    // Left generator: vb 閳?vc (up along z at 鑳?锜?
    let left_gen = make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: vb, direction: z_ax }),
        0.0, h, vi_vb, vi_vc,
    ).ok()?;

    // Top arc at z = +ew: from vc 閳?vd along y閳? (鑳? 锜洪埆? on Circle3 normal=+Z)
    let top_center = DVec3::new(cyl_cx, cyl_cy, ew);
    let top_arc = make_edge(
        &mut brep,
        Curve3::Circle(Circle3::new(top_center, DVec3::Z, cyl_r )),
        pi, 0.0, vi_vc, vi_vd,
    ).ok()?;

    // Right generator: vd 閳?va (down along z at 鑳?0, stored as va閳姵d, rev 閳?vd閳姵a)
    let right_gen = make_edge(
        &mut brep,
        Curve3::Line(Line3 { origin: va, direction: z_ax }),
        0.0, h, vi_va, vi_vd,
    ).ok()?;

    // Push empty edge-pcurve vecs for each new edge (make_edge does NOT
    // push to edge_pcurves).
    for _ in 0..4 {
        brep.geom.edge_pcurves.push(Vec::new());
    }

    // 4. Create CylindricalSurface for the wall.
    let surf = Surface3::Cylinder(CylindricalSurface {
        origin: DVec3::new(cyl_cx, cyl_cy, -ew),
        axis: DVec3::Z,
        radius: cyl_r,
        ref_dir: DVec3::X,
    });
    let si = brep.geom.surfaces.len();
    brep.geom.surfaces.push(surf);

    // 5. Create the wall face with wire.
    // Wire: va 閳?vb (bot_arc_fwd), vb 閳?vc (left_gen_fwd),
    //       vc 閳?vd (top_arc_fwd), vd 閳?va (right_gen_rev)
    let wall_wire = Wire {
        edges: vec![
            WireEdge::fwd(bot_arc),
            WireEdge::fwd(left_gen),
            WireEdge::fwd(top_arc),
            WireEdge::rev(right_gen),
        ],
    };
    let fi = brep.solids[0].shells[0].faces.len();
    brep.solids[0].shells[0].faces.push(Face {
        outer_wire: wall_wire,
        inner_wires: Vec::new(),
        normal: DVec3::ZERO,
        triangles: Vec::new(),
        sample_point: None,
        mesh_dirty: true,
                surface_idx: None,
    });
    while brep.geom.face_surface.len() <= fi {
        brep.geom.face_surface.push(None);
    }
    brep.geom.face_surface[fi] = Some(si);
    while brep.geom.face_surface_range.len() <= fi {
        brep.geom.face_surface_range.push(None);
    }
    // UV range: u 閳?[锜? 2锜篯 (y閳? half of CylindricalSurface ref_dir=X),
    //           v 閳?[0, h]
    brep.geom.face_surface_range[fi] = Some([pi, 2.0 * pi, 0.0, h]);

    Some(brep)
}

/// Fast path: cylinder 閳?box boolean intersection.
///
/// When the box fully contains the cylinder cross-section in XY (at the Z-range
/// of intersection), the result is simply the cylinder clipped to the box's
/// Z range 閳?avoiding the Pave-Filler which produces hundreds of unnecessary
/// faces for symmetric configurations (e.g. bopcommon_simple V5/V6).
pub fn try_intersection_cylinder_box(a: &BRep, b: &BRep) -> Option<BRep> {
    // Try both orderings.
    try_intersect_cylinder_box_one_dir(a, b)
        .or_else(|| try_intersect_cylinder_box_one_dir(b, a))
}

/// Helper: find which theta values in [0, 2锜? satisfy all clip-plane constraints.
/// For clip plane (inward_normal n, cut_dist d): valid where cos(鑳冮埈鎹? 閳?閳妿/r.
/// Returns sorted disjoint intervals within [0, 2锜?.
pub(crate) fn compute_valid_theta_ranges(r: f64, clip_planes: &[(DVec3, f64)]) -> Vec<(f64, f64)> {
    let pi = std::f64::consts::PI;
    let two_pi = 2.0 * pi;

    let mut valid = vec![(0.0, two_pi)];

    for &(n, d) in clip_planes {
        let cd = (-d / r).clamp(-1.0, 1.0);
        let alpha = cd.acos();
        if alpha >= pi - 1e-12 {
            continue; // full circle, no constraint
        }
        let phi = n.y.atan2(n.x);

        // Constraint: 鑳?閳?[锠侀埈鎹? 锠?浼猐 mod 2锜?
        let lo = phi - alpha;
        let hi = phi + alpha;

        let constraint = {
            let lo_norm = lo.rem_euclid(two_pi);
            let hi_norm = hi.rem_euclid(two_pi);
            if lo_norm <= hi_norm {
                vec![(lo_norm, hi_norm)]
            } else {
                vec![(lo_norm, two_pi), (0.0, hi_norm)]
            }
        };

        let mut next = Vec::new();
        for &(vl, vr) in &valid {
            for &(cl, cr) in &constraint {
                let l = f64::max(vl, cl);
                let r = f64::min(vr, cr);
                if r > l + 1e-12 {
                    next.push((l, r));
                }
            }
        }
        valid = next;
        if valid.is_empty() {
            break;
        }
    }
    valid
}

fn split_theta_ranges(intervals: &[(f64, f64)], split_thetas: &[f64]) -> Vec<(f64, f64)> {
    let two_pi = 2.0 * std::f64::consts::PI;
    if intervals.is_empty() {
        return Vec::new();
    }

    let mut splits: Vec<f64> = split_thetas.iter()
        .map(|theta| {
            let t = theta.rem_euclid(two_pi);
            if t >= two_pi - 1e-12 { 0.0 } else { t }
        })
        .collect();
    splits.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    splits.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

    let mut result = Vec::new();
    for &(lo, hi) in intervals {
        let mut cuts = vec![lo];
        for &split in &splits {
            if split > lo + 1e-12 && split < hi - 1e-12 {
                cuts.push(split);
            }
        }
        cuts.push(hi);
        for pair in cuts.windows(2) {
            if pair[1] > pair[0] + 1e-12 {
                result.push((pair[0], pair[1]));
            }
        }
    }
    result
}

fn cylinder_box_tangent_split_thetas(
    full_u: bool,
    full_v: bool,
    cu: f64,
    cv: f64,
    u_ax: DVec3,
    v_ax: DVec3,
    eu: f64,
    ev: f64,
    cyl_r: f64,
    tol: f64,
) -> Vec<f64> {
    let mut split_thetas = Vec::new();
    for &(full_axis, cp, ax, ext) in &[(full_u, cu, u_ax, eu), (full_v, cv, v_ax, ev)] {
        let _ = full_axis;
        if (cp - cyl_r + ext).abs() <= tol {
            split_thetas.push((-ax).y.atan2((-ax).x));
        }
        if (cp + cyl_r - ext).abs() <= tol {
            split_thetas.push(ax.y.atan2(ax.x));
        }
    }
    split_thetas
}

/// Compute the complement (inverse) of a set of 鑳?intervals within [0, 2锜?.
/// Each interval (lo, hi) is assumed sorted and non-overlapping.
/// The result covers the parts of [0, 2锜? not in any input interval.
pub(crate) fn compute_complement_theta_ranges(intervals: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let two_pi = 2.0 * std::f64::consts::PI;
    if intervals.is_empty() {
        return vec![(0.0, two_pi)];
    }
    // Sort by lo first 閳?compute_valid_theta_ranges may return intervals in
    // insertion order (depends on the order clip planes were pushed), not
    // sorted by angle. Without sorting, the sequential sweep produces wrong
    // complement when a later interval has a smaller lo.
    let mut sorted: Vec<(f64, f64)> = intervals.to_vec();
    sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut result = Vec::new();
    let mut prev = 0.0;
    for &(lo, hi) in &sorted {
        if lo > prev + 1e-12 {
            result.push((prev, lo));
        }
        prev = prev.max(hi);
    }
    if two_pi - prev > 1e-12 {
        result.push((prev, two_pi));
    }
    result
}

/// On which clip-plane boundary (index into `info`) does `theta` lie?
fn find_plane_for_theta(theta: f64, info: &[(f64, f64, DVec3, f64)], r: f64) -> Option<usize> {
    for (i, &(phi, alpha, _n, d)) in info.iter().enumerate() {
        if alpha >= std::f64::consts::PI - 1e-12 {
            continue;
        }
        let diff = (theta - phi).cos() + d / r;
        if diff.abs() < 1e-8 {
            return Some(i);
        }
    }
    None
}

/// Solve for the XY corner where two clip planes intersect.
/// Each plane's line: n璺?p閳妽enter.xy) = 閳妿  閳? n璺痯 = n璺痗enter.xy 閳?d
fn corner_of_planes(
    n1: DVec3, d1: f64,
    n2: DVec3, d2: f64,
    center: DVec3,
) -> DVec3 {
    // 2鑴? system: [n1.x n1.y; n2.x n2.y] * p = [n1璺痗_xy 閳?d1; n2璺痗_xy 閳?d2]
    let cx = center.x;
    let cy = center.y;
    let a = n1.x; let b = n1.y;
    let c = n2.x; let d = n2.y;
    let rhs1 = cx * n1.x + cy * n1.y - d1;
    let rhs2 = cx * n2.x + cy * n2.y - d2;
    let det = a * d - b * c;
    if det.abs() < 1e-15 {
        // Fallback: average the two plane origins
        let o1 = DVec3::new(cx - n1.x * d1, cy - n1.y * d1, 0.0);
        let o2 = DVec3::new(cx - n2.x * d2, cy - n2.y * d2, 0.0);
        return (o1 + o2) * 0.5;
    }
    let inv = 1.0 / det;
    DVec3::new(
        (rhs1 * d - rhs2 * b) * inv,
        (a * rhs2 - c * rhs1) * inv,
        0.0,
    )
}

/// Check if a point in XY satisfies all clip-plane constraints (within tolerance).
fn point_satisfies_all(p: DVec3, clip_planes: &[(DVec3, f64)], center: DVec3) -> bool {
    let tol = 1e-8;
    for &(n, d) in clip_planes {
        if n.dot(p - center) < -d - tol {
            return false;
        }
    }
    true
}

/// Build the shortest sequence of plane indices from `p_from` to `p_to`
/// where each consecutive pair has a valid (non-parallel, constraint-satisfying) corner.
///
/// The chain is used to route the gap boundary path through intermediate clip planes.
fn build_plane_chain(
    p_from: usize, p_to: usize,
    clip_planes: &[(DVec3, f64)],
    center: DVec3,
) -> Vec<usize> {
    if p_from == p_to {
        return vec![p_from];
    }

    // Build chain through intermediate planes sorted by inward normal angle.
    let n_planes = clip_planes.len();
    let mut indices: Vec<usize> = (0..n_planes).collect();
    indices.sort_by(|&i, &j| {
        let ai = clip_planes[i].0.y.atan2(clip_planes[i].0.x);
        let aj = clip_planes[j].0.y.atan2(clip_planes[j].0.x);
        ai.partial_cmp(&aj).unwrap()
    });

    let pos_from = indices.iter().position(|&i| i == p_from).unwrap();
    let _pos_to = indices.iter().position(|&i| i == p_to).unwrap();

    // Build forward chain (increasing angle, wraps around)
    let mut fwd: Vec<usize> = Vec::new();
    let mut i = pos_from;
    loop {
        fwd.push(indices[i]);
        if indices[i] == p_to { break; }
        i = (i + 1) % n_planes;
        if i == pos_from { break; } // full circle
    }

    // Build backward chain (decreasing angle, wraps around)
    let mut bwd: Vec<usize> = Vec::new();
    let mut i = pos_from;
    loop {
        bwd.push(indices[i]);
        if indices[i] == p_to { break; }
        i = if i == 0 { n_planes - 1 } else { i - 1 };
        if i == pos_from { break; }
    }

    // Validate a chain: every consecutive pair must have a non-parallel
    // corner that satisfies all constraints.
    let valid_chain = |chain: &[usize]| -> bool {
        for w in chain.windows(2) {
            let (pa, pb) = (w[0], w[1]);
            let (na, da) = (clip_planes[pa].0, clip_planes[pa].1);
            let (nb, db) = (clip_planes[pb].0, clip_planes[pb].1);
            let det = na.x * nb.y - na.y * nb.x;
            if det.abs() < 1e-12 { return false; }
            let c = corner_of_planes(na, da, nb, db, center);
            if !point_satisfies_all(c, clip_planes, center) { return false; }
        }
        true
    };

    let fwd_ok = valid_chain(&fwd);
    let bwd_ok = valid_chain(&bwd);

    if fwd_ok && bwd_ok {
        // Both chains are valid.  The fwd chain follows increasing 锠?order and
        // includes all intermediate clip planes (always correct for box faces).
        // The bwd chain may skip intermediate planes, producing a shortcut
        // chord that cuts through the interior of the valid region instead of
        // tracing the full boundary.  Always prefer fwd when both are valid.
        fwd
    } else if fwd_ok {
        fwd
    } else if bwd_ok {
        bwd
    } else {
        // No valid chain found 閳?fall back to direct (may produce wrong geometry,
        // but prevents a crash / empty gap).
        vec![p_from, p_to]
    }
}

/// Build a BRep for the intersection of a Z-aligned cylinder with up to 4
/// vertical clip planes (box faces parallel to the cylinder axis).
///
/// Each clip plane is `(inward_normal, cut_dist)`.  The result is the portion
/// of the cylinder that satisfies ALL clip-plane constraints simultaneously.
///
/// The result may have multiple cylindrical-wall faces (one per valid thetas run),
/// top/bottom planar caps with wire boundaries that alternate between circle
/// arcs and clip-plane chord segments, and one rectangular side face per
/// clip plane.
fn build_cylinder_box_intersection_brep(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
) -> BRep {
    build_cylinder_box_intersection_brep_with_splits(center, r, h, clip_planes, &[])
}

fn build_cylinder_box_intersection_brep_with_splits(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
    split_thetas: &[f64],
) -> BRep {
    let intervals = split_theta_ranges(&compute_valid_theta_ranges(r, clip_planes), split_thetas);
    build_cylinder_box_clipped_brep(center, r, h, &intervals, clip_planes, false, false, true)
}

/// Build the middle Z-slice of C \ B (cylinder minus box) for partial XY overlap.
///
/// The result is the portion of the cylinder at the box's Z-range, with the
/// box's XY cross-section removed.  The wall consists of the complement of the
/// valid (inside-box) theta intervals.  Caps alternate between complement arcs
/// and chords on the clip planes.  Side faces are generated on each clip plane
/// from gap segments.
///
/// When the valid intervals are empty (box entirely within the cylinder
/// cross-section, e.g. bopcut_simple V1 where the box is the inscribed square),
/// the complement covers [0, 2pi) and the result is a full cylinder wall with
/// donut caps (full circle minus box polygon) and side faces on each clip plane.
fn build_cylinder_box_difference_middle(
    center: DVec3,
    r: f64,
    h: f64,
    clip_planes: &[(DVec3, f64)],
    skip_bottom_cap: bool,
    skip_top_cap: bool,
) -> BRep {
    let intersection_intervals = compute_valid_theta_ranges(r, clip_planes);
    if intersection_intervals.is_empty() {
        // Full-wall case: box entirely within cylinder cross-section.
        // Build full cylinder wall + donut caps + side faces on clip planes.
        return build_cylinder_box_difference_full_wall(center, r, h, clip_planes);
    }
    // Check if all clip planes are parallel (e.g. cylinder center inside box
    // in one axis).  The downstream gap routing cannot handle parallel-only
    // planes because `build_plane_chain` requires non-parallel corners.
    let all_parallel = clip_planes.len() >= 2
        && clip_planes.iter().skip(1).all(|(n, _)| {
            let (n0, _) = clip_planes[0];
            (n0.x * n.y - n0.y * n.x).abs() < 1e-12
        });
    if all_parallel {
        return build_cylinder_box_difference_parallel_only_skip(center, r, h, clip_planes, skip_bottom_cap, skip_top_cap);
    }
    let complement = compute_complement_theta_ranges(&intersection_intervals);
    build_cylinder_box_clipped_brep(center, r, h, &complement, clip_planes, skip_bottom_cap, skip_top_cap, false)
}