
/// Build a tessellated BRep for coaxial `cylinder Èñ?torus` where the torus
/// sits on the cylinder wall (same Z axis, same XY center).
///
/// Cylinder: [cyl_z_lo, cyl_z_hi], radius cyl_r.
/// Torus: centered at (0,0,tor_z), major radius R, minor radius r_m.
///
/// Builds via Z-slice tessellation: at each Z the cross-section is a circle
/// with radius = max(cyl_r, R + sqrt(r_mÈì?Èñ?(z Èñ?tor_z)Èì?).
fn build_cylinder_torus_union_tessellated(
    center_xy: DVec2,
    cyl_z_lo: f64, cyl_z_hi: f64, cyl_r: f64,
    tor_z: f64, R: f64, r_m: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    if cyl_r < tol || cyl_z_hi <= cyl_z_lo + tol || r_m < tol { return None; }

    let n_slices = 64usize;
    let n_slices_circ = 16usize;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

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
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    let circle_poly = |r: f64| -> Vec<DVec2> {
        let mut poly = Vec::with_capacity(n_arc + 1);
        for k in 0..=n_arc {
            let ang = tau * k as f64 / n_arc as f64;
            let (s, c) = ang.sin_cos();
            poly.push(DVec2::new(r * c, r * s));
        }
        poly
    };

    let radius_at = |z: f64| -> f64 {
        let dz = (z - tor_z).abs();
        if dz <= r_m {
            let bulge = R + (r_m * r_m - dz * dz).sqrt();
            cyl_r.max(bulge)
        } else {
            cyl_r
        }
    };

    // Bulge Z-range (clamped to cylinder)
    let bulge_z_lo = (tor_z - r_m).max(cyl_z_lo);
    let bulge_z_hi = (tor_z + r_m).min(cyl_z_hi);
    let has_bulge = bulge_z_hi > bulge_z_lo + tol;

    // Constant-radius circle (cylinder wall)
    let cyl_poly = circle_poly(cyl_r);

    // 1. Below bulge (constant radius)
    if bulge_z_lo > cyl_z_lo + tol {
        add_wall_section(&mut add_v, &mut faces, &cyl_poly, cyl_z_lo, bulge_z_lo, n_slices_circ, &to_world, &empty_wire);
    }

    // 2. Bulge section (varying radius via Z-slices)
    if has_bulge {
        let dz = (bulge_z_hi - bulge_z_lo) / n_slices as f64;
        for i in 0..n_slices {
            let za = bulge_z_lo + dz * i as f64;
            let zb = bulge_z_lo + dz * (i + 1) as f64;
            let ra = radius_at(za).max(0.0);
            let rb = radius_at(zb).max(0.0);
            if ra < tol && rb < tol { continue; }

            let nn = n_arc;
            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(ra * c, ra * s, za)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(rb * c, rb * s, zb)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        }
    }

    // 3. Above bulge (constant radius)
    if cyl_z_hi > bulge_z_hi + tol {
        add_wall_section(&mut add_v, &mut faces, &cyl_poly, bulge_z_hi, cyl_z_hi, n_slices_circ, &to_world, &empty_wire);
    }

    // Caps at cylinder ends (radius may differ from cyl_r if bulge extends to end)
    let bottom_r = radius_at(cyl_z_lo);
    let top_r = radius_at(cyl_z_hi);
    let bottom_poly = if (bottom_r - cyl_r).abs() < tol { cyl_poly.clone() } else { circle_poly(bottom_r) };
    let top_poly = if (top_r - cyl_r).abs() < tol { cyl_poly.clone() } else { circle_poly(top_r) };

    add_cap_face(&mut add_v, &mut faces, &bottom_poly, cyl_z_lo, -DVec3::Z, &to_world, &empty_wire);
    add_cap_face(&mut add_v, &mut faces, &top_poly, cyl_z_hi, DVec3::Z, &to_world, &empty_wire);

    if faces.is_empty() { return None; }

    let geom = GeomStore { face_internal_vertices: vec![],
        curves: vec![], surfaces: vec![], curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; faces.len()],
        edge_pcurves: vec![], edge_curve_range: vec![],
        edge_degenerated: vec![], vertex_tolerance: vec![],
        edge_tolerance: vec![], face_tolerance: vec![],
        curve2d_range: vec![], face_surface_range: vec![None; faces.len()],
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Fast path: coaxial Z-aligned cylinder + torus Union.
///
/// Detects a Z-aligned cylinder and a Z-aligned torus sharing the same XY
/// center, where the torus major radius Èñ?cylinder radius (torus protrudes
/// from or is flush with the cylinder wall).
fn try_union_cylinder_torus_one_dir(cyl_brep: &BRep, torus_brep: &BRep) -> Option<BRep> {
    let (cyl_bottom, cyl_axis, cyl_r, cyl_h) = try_cylinder_center_axis_radius_height(cyl_brep)?;
    let (tor_center, tor_axis, R, r_m) = torus_info(torus_brep)?;

    // Both must be Z-aligned
    if cyl_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN { return None; }
    if tor_axis.dot(DVec3::Z).abs() < 1.0 - TOLERANCE_AXIS_ALIGN { return None; }

    // Coaxial: same XY center
    let tol = TOLERANCE_LEN_MIN;
    if (cyl_bottom.x - tor_center.x).abs() > tol || (cyl_bottom.y - tor_center.y).abs() > tol {
        return None;
    }

    let cyl_z_lo = cyl_bottom.z;
    let cyl_z_hi = cyl_bottom.z + cyl_h;
    let tor_z = tor_center.z;

    // Torus must overlap with cylinder Z-range
    if tor_z + r_m < cyl_z_lo + tol || tor_z - r_m > cyl_z_hi - tol {
        return None;
    }

    // Torus major radius must be at least approximately the cylinder radius
    // (so the torus protrudes from or is flush with the wall)
    if R < cyl_r - tol {
        return None;
    }

    build_cylinder_torus_union_tessellated(
        DVec2::new(tor_center.x, tor_center.y),
        cyl_z_lo, cyl_z_hi, cyl_r,
        tor_z, R, r_m,
    )
}

/// Bidirectional wrapper for cylinder + torus Union.
pub fn try_union_cylinder_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_cylinder_torus_one_dir(a, b)
        .or_else(|| try_union_cylinder_torus_one_dir(b, a))
}

/// Fast path: Union of two same-center tori (same R, r).
///
/// Detects two single-face torus primitives at the same center with matching
/// radii, and returns a BRep containing both full torus faces AS-IS (no
/// trimming).  The actual union trimming is deferred to a subsequent boolean
/// step, where the Pave-Filler processes each torus face against the third
/// operand, and the A-A overlap check in the BooleanBuilder removes sub-faces
/// that are inside the other original torus.
fn try_union_torus_torus_one_dir(a: &BRep, b: &BRep) -> Option<BRep> {
    // Only fast-path single-face tori.  Multi-face operands (e.g. the output of
    // an earlier fast-path union) must go through the Pave-Filler for correct
    // trimming and surface-area computation.
    let total_faces = |brep: &BRep| -> usize {
        brep.solids.iter().flat_map(|s| s.shells.iter()).flat_map(|sh| &sh.faces).count()
    };
    if total_faces(a) != 1 || total_faces(b) != 1 {
        return None;
    }

    let (center_a, _axis_a, R_a, r_a) = torus_info(a)?;
    let (center_b, _axis_b, R_b, r_b) = torus_info(b)?;

    // Same center (within tolerance)
    if (center_a - center_b).length() > TOLERANCE_LEN_MIN {
        return None;
    }
    // Matching major and minor radii
    if (R_a - R_b).abs() > TOLERANCE_LEN_MIN || (r_a - r_b).abs() > TOLERANCE_LEN_MIN {
        return None;
    }

    // Run through the full boolean pipeline to produce a properly merged single-solid
    // BRep.  A multi-solid result would cause the Pave-Filler to miss A-A face pairs
    // when this result is used as an operand in a subsequent boolean operation.
    let mut result = a.clone();
    result.append_disjoint_brep(b);
    Some(result)
}

/// Bidirectional wrapper for torus + torus Union.
pub fn try_union_torus_torus(a: &BRep, b: &BRep) -> Option<BRep> {
    let mut tori = collect_same_center_tori(a)?;
    tori.extend(collect_same_center_tori(b)?);
    if tori.len() >= 3 {
        if let Some(mesh) = build_same_center_tori_union_mesh(&tori) {
            return Some(mesh);
        }
    }

    try_union_torus_torus_one_dir(a, b)
        .or_else(|| try_union_torus_torus_one_dir(b, a))
}

fn collect_same_center_tori(brep: &BRep) -> Option<Vec<ToroidalSurface>> {
    let mut tori = Vec::new();
    let mut flat_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for _face in &shell.faces {
                let surf_idx = brep.geom.face_surface.get(flat_idx).and_then(|o| *o)?;
                let Surface3::Torus(torus) = brep.geom.surfaces.get(surf_idx)? else {
                    return None;
                };
                tori.push(*torus);
                flat_idx += 1;
            }
        }
    }

    if tori.is_empty() {
        return None;
    }

    let first = tori[0];
    if !tori.iter().all(|t| {
        (t.center - first.center).length() <= TOLERANCE_LEN_MIN
            && (t.major_radius - first.major_radius).abs() <= TOLERANCE_LEN_MIN
            && (t.minor_radius - first.minor_radius).abs() <= TOLERANCE_LEN_MIN
    }) {
        return None;
    }

    Some(tori)
}

fn point_inside_torus_solid(torus: &ToroidalSurface, p: DVec3, tol: f64) -> bool {
    let axis = torus.axis.normalize_or(DVec3::Z);
    let local = p - torus.center;
    let axial = local.dot(axis);
    let radial = local - axis * axial;
    let tube_dist_sq = (radial.length() - torus.major_radius).powi(2) + axial * axial;
    tube_dist_sq <= (torus.minor_radius + tol).powi(2)
}

fn build_same_center_tori_union_mesh(tori: &[ToroidalSurface]) -> Option<BRep> {
    if tori.len() < 3 {
        return None;
    }

    let n_u = 192usize;
    let n_v = 96usize;
    let tau = std::f64::consts::TAU;
    let mut vertices: Vec<Vertex> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    let tol = TOLERANCE_ABS * tori[0].major_radius.max(tori[0].minor_radius).max(1.0);

    let push_vertex = |p: DVec3, vertices: &mut Vec<Vertex>| -> usize {
        let idx = vertices.len();
        vertices.push(Vertex { point: p });
        idx
    };

    for (ti, torus) in tori.iter().enumerate() {
        for iu in 0..n_u {
            let u0 = tau * iu as f64 / n_u as f64;
            let u1 = tau * (iu + 1) as f64 / n_u as f64;
            let uc = 0.5 * (u0 + u1);
            for iv in 0..n_v {
                let v0 = tau * iv as f64 / n_v as f64;
                let v1 = tau * (iv + 1) as f64 / n_v as f64;
                let vc = 0.5 * (v0 + v1);
                let sample = torus.point_at(uc, vc);
                if tori.iter().enumerate().any(|(oi, other)| {
                    oi != ti && point_inside_torus_solid(other, sample, tol)
                }) {
                    continue;
                }

                let p00 = torus.point_at(u0, v0);
                let p10 = torus.point_at(u1, v0);
                let p11 = torus.point_at(u1, v1);
                let p01 = torus.point_at(u0, v1);
                let i00 = push_vertex(p00, &mut vertices);
                let i10 = push_vertex(p10, &mut vertices);
                let i11 = push_vertex(p11, &mut vertices);
                let i01 = push_vertex(p01, &mut vertices);
                triangles.push([i00, i10, i11]);
                triangles.push([i00, i11, i01]);
            }
        }
    }

    if triangles.is_empty() {
        return None;
    }

    let face = Face {
        outer_wire: Wire { edges: vec![] },
        inner_wires: vec![],
        normal: DVec3::Z,
        triangles,
        sample_point: None,
        mesh_dirty: false,
                surface_idx: None,
    };

    let mut brep = BRep {
        vertices,
        edges: vec![],
        solids: vec![Solid {
            shells: vec![Shell { faces: vec![face] }],
        }],
        geom: GeomStore::default(),
        compound: None,
        compsolid: None,
    };
    rcad_kernel::resize_tolerance_arrays(&mut brep);
    Some(brep)
}

/// Build a tessellated BRep for the union of two coaxial Z-aligned cones.
///
/// Each cone i spans [z_i_lo, z_i_hi] with radius r_i(z) = r_i_lo +
/// (r_i_hi Èñ?r_i_lo)Áí?z Èñ?z_i_lo)/(z_i_hi Èñ?z_i_lo).  The outer envelope
/// at each Z is max(rÈñ?z), rÈñ?z), 0).  Ring faces are added at boundaries
/// where the envelope radius changes discontinuously.
fn build_coaxial_cones_union_tessellated(
    center_xy: DVec2,
    z1_lo: f64, z1_hi: f64, r1_lo: f64, r1_hi: f64,
    z2_lo: f64, z2_hi: f64, r2_lo: f64, r2_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let n_arc = 128usize;
    let tau = std::f64::consts::TAU;
    let empty_wire = || Wire { edges: vec![] };

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
        DVec3::new(center_xy.x + u, center_xy.y + v, z)
    };

    // Radius function for each cone, returns 0 outside its Z range
    let cone_r = |z: f64, z_lo: f64, z_hi: f64, r_lo: f64, r_hi: f64| -> f64 {
        if z < z_lo - tol || z > z_hi + tol { return 0.0; }
        let zc = z.clamp(z_lo, z_hi);
        let t = if z_hi > z_lo { (zc - z_lo) / (z_hi - z_lo) } else { 0.0 };
        (r_lo + (r_hi - r_lo) * t).max(0.0)
    };

    // Outer envelope radius at Z
    let env_r = |z: f64| -> f64 {
        cone_r(z, z1_lo, z1_hi, r1_lo, r1_hi)
            .max(cone_r(z, z2_lo, z2_hi, r2_lo, r2_hi))
    };

    let z_union_lo = z1_lo.min(z2_lo);
    let z_union_hi = z1_hi.max(z2_hi);
    if z_union_hi <= z_union_lo + tol { return None; }

    // Collect all unique Z boundaries
    let mut bounds: Vec<f64> = vec![z_union_lo, z_union_hi];
    for z in &[z1_lo, z1_hi, z2_lo, z2_hi] {
        if *z > z_union_lo + tol && *z < z_union_hi - tol {
            bounds.push(*z);
        }
    }
    bounds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    bounds.dedup();

    let nn = n_arc;
    // Build wall for each segment, then ring at each interior boundary
    for bi in 0..bounds.len() - 1 {
        let za = bounds[bi];
        let zb = bounds[bi + 1];
        if zb - za < tol { continue; }

        let ra_mid = env_r((za + zb) * 0.5);
        if ra_mid < tol { continue; } // no material in this segment

        // Build varying-radius wall for this segment
        let n_slices = 16usize.max(1);
        let dz = (zb - za) / n_slices as f64;
        for i in 0..n_slices {
            let z0 = za + dz * i as f64;
            let z1 = za + dz * (i + 1) as f64;
            let r0 = env_r(z0).max(0.0);
            let r1 = env_r(z1).max(0.0);
            if r0 < tol && r1 < tol { continue; }

            let mut idx = Vec::with_capacity(2 * (nn + 1));
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(r0 * c, r0 * s, z0)));
            }
            for k in 0..=nn {
                let ang = tau * k as f64 / nn as f64;
                let (s, c) = ang.sin_cos();
                idx.push(add_v(to_world(r1 * c, r1 * s, z1)));
            }

            let mut tris = Vec::with_capacity(nn * 2);
            for j in 0..nn {
                tris.push([idx[j], idx[j + 1], idx[nn + 1 + j + 1]]);
                tris.push([idx[j], idx[nn + 1 + j + 1], idx[nn + 1 + j]]);
            }
            faces.push(Face {
                outer_wire: empty_wire(), inner_wires: vec![],
                normal: DVec3::ZERO, triangles: tris,
                sample_point: None, mesh_dirty: false,
                surface_idx: None,
            });
        }

        // Ring face at the end of this segment (check radius discontinuity)
        if bi + 1 < bounds.len() - 1 {
            let z_bound = zb;
            let r_left = env_r(z_bound - tol * 0.5);
            let r_right = env_r(z_bound + tol * 0.5);
            let r_min = r_left.min(r_right);
            let r_max = r_left.max(r_right);
            if r_max - r_min > tol && r_min > tol {
                // Add annular ring
                let r_outer = r_max;
                let r_inner = r_min;
                let mut outer_idx: Vec<usize> = (0..nn).map(|i| {
                    let ang = tau * i as f64 / nn as f64;
                    let (s, c) = ang.sin_cos();
                    add_v(to_world(r_outer * c, r_outer * s, z_bound))
                }).collect();
                let inner_idx: Vec<usize> = (0..nn).map(|i| {
                    let ang = tau * i as f64 / nn as f64;
                    let (s, c) = ang.sin_cos();
                    add_v(to_world(r_inner * c, r_inner * s, z_bound))
                }).collect();
                outer_idx.extend(&inner_idx);
                let mut tris = Vec::with_capacity(nn * 2);
                for i in 0..nn {
                    let k = (i + 1) % nn;
                    tris.push([outer_idx[i], outer_idx[k], inner_idx[k]]);
                    tris.push([outer_idx[i], inner_idx[k], inner_idx[i]]);
                }
                faces.push(Face {
                    outer_wire: empty_wire(), inner_wires: vec![],
                    normal: DVec3::ZERO, triangles: tris,
                    sample_point: None, mesh_dirty: false,
                surface_idx: None,
                });
            }
        }
    }

    // Caps
    let r_bot = env_r(z_union_lo);
    if r_bot > tol {
        let bot_poly: Vec<DVec2> = (0..=nn).map(|i| {
            let ang = tau * i as f64 / nn as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(r_bot * c, r_bot * s)
        }).collect();
        add_cap_face(&mut add_v, &mut faces, &bot_poly, z_union_lo, -DVec3::Z, &to_world, &empty_wire);
    }
    let r_top = env_r(z_union_hi);
    if r_top > tol {
        let top_poly: Vec<DVec2> = (0..=nn).map(|i| {
            let ang = tau * i as f64 / nn as f64;
            let (s, c) = ang.sin_cos();
            DVec2::new(r_top * c, r_top * s)
        }).collect();
        add_cap_face(&mut add_v, &mut faces, &top_poly, z_union_hi, DVec3::Z, &to_world, &empty_wire);
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
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Fast path: coaxial Z-aligned cone frustums Union.
///
/// Detects two Z-aligned conical frustums sharing the same XY center and
/// builds the union via Z-slice tessellation of the outer radius envelope.
fn try_union_coaxial_cones_one_dir(a: &BRep, b: &BRep) -> Option<BRep> {
    let (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi) = detect_z_axis_cone(a)?;
    let (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi) = detect_z_axis_cone(b)?;

    // Coaxial
    let tol = TOLERANCE_LEN_MIN;
    if (c1_xy.x - c2_xy.x).abs() > tol || (c1_xy.y - c2_xy.y).abs() > tol {
        return None;
    }

    // Z ranges must overlap or be adjacent
    if c1_z_hi < c2_z_lo - tol && (c2_z_lo - c1_z_hi) > tol * 100.0 { return None; }
    if c2_z_hi < c1_z_lo - tol && (c1_z_lo - c2_z_hi) > tol * 100.0 { return None; }

    build_coaxial_cones_union_tessellated(
        c1_xy,
        c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi,
        c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi,
    )
}

/// Bidirectional wrapper for coaxial cone + cone Union.
pub fn try_union_coaxial_cones(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_coaxial_cones_one_dir(a, b)
        .or_else(|| try_union_coaxial_cones_one_dir(b, a))
}

/// Return the CCW boundary polygon of the union of two circles.
///
/// `c1`, `c2` Èñ?centers, `r1`, `r2` Èñ?radii.
/// `n1` Èñ?sample count for C1's boundary arc, `n2` Èñ?for C2's arc.
///
/// Cases handled: concentric (larger circle), contained (larger circle),
/// disjoint (both full circles), and overlapping (two arcs joined).
fn two_circle_union_ccw_pts(
    c1: DVec2, r1: f64, c2: DVec2, r2: f64, n1: usize, n2: usize,
) -> Vec<DVec2> {
    let d = (c2 - c1).length();
    let tol = 1e-12;

    if d < tol || d + r1.min(r2) <= r1.max(r2) + tol {
        // Concentric or one contains the other Èñ?return the larger circle full.
        let (_c, r) = if r1 >= r2 { (c1, r1) } else { (c2, r2) };
        let tau = std::f64::consts::TAU;
        let n = n1 + n2;
        return (0..=n).map(|i| {
            let ang = tau * i as f64 / n as f64;
            let (s, c) = ang.sin_cos();
            c + DVec2::new(r * c, r * s)
        }).collect();
    }

    if d >= r1 + r2 - tol {
        // Disjoint Èñ?both full circles.
        let tau = std::f64::consts::TAU;
        let mut pts: Vec<DVec2> = (0..n1).map(|i| {
            let ang = tau * i as f64 / n1 as f64;
            let (s, c) = ang.sin_cos();
            c1 + DVec2::new(r1 * c, r1 * s)
        }).collect();
        let _last = *pts.last().unwrap_or(&c1);
        pts.extend((0..=n2).map(|i| {
            let ang = tau * i as f64 / n2 as f64;
            let (s, c) = ang.sin_cos();
            c2 + DVec2::new(r2 * c, r2 * s)
        }));
        // Drop duplicate start point
        if pts.len() > 3 && (pts.last().unwrap() - pts[0]).length_squared() < tol * tol {
            pts.pop();
        }
        return pts;
    }

    // Overlapping Èñ?trace the outer envelope arcs.
    let dir = (c2 - c1) / d;          // unit vector C1 Èñ?C2
    let perp = DVec2::new(-dir.y, dir.x); // 90Èé?CCW

    let ix = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let _iy = (r1 * r1 - ix * ix).max(0.0).sqrt();

    let cos_t1 = (ix / r1).clamp(-1.0, 1.0);
    let theta1 = cos_t1.acos();                  // C1 arc half-angle
    let C = ((r1 * r1 - d * d - r2 * r2) / (2.0 * d * r2)).clamp(-1.0, 1.0);
    let theta2 = C.acos();                       // C2 arc half-angle

    let tau = std::f64::consts::TAU;

    // C1 outer arc: from theta1 Èñ?2Èî?Èñ?theta1 (through Èî? the left/away side).
    let a1_start = theta1;
    let a1_end = tau - theta1;
    let sweep1 = ((a1_end - a1_start) % tau + tau) % tau; // positive

    let mut pts = Vec::with_capacity(n1 + n2 + 2);
    for k in 0..=n1 {
        let frac = k as f64 / n1 as f64;
        let ang = a1_start + sweep1 * frac;
        let (s, c) = ang.sin_cos();
        pts.push(c1 + DVec2::new(r1 * c, r1 * s));
    }

    // C2 outer arc: from Èñ≥ÓÖüÂ¥? Èñ?+Èë? (through 0, the right/away side).
    // In C2's local frame (+X = dir = away from C1, +Y = perp = CCW).
    let a2_start = -theta2;
    let a2_end = theta2;
    let sweep2 = a2_end - a2_start; // = 2 * theta2, always positive

    for k in 1..n2 {
        let frac = k as f64 / n2 as f64;
        let ang = a2_start + sweep2 * frac;
        let (s, c) = ang.sin_cos();
        pts.push(c2 + r2 * (c * dir + s * perp));
    }

    // Close the polygon (remove trailing duplicate start-point).
    if pts.len() > 3 && (pts.last().unwrap() - pts[0]).length_squared() < tol * tol {
        pts.pop();
    }
    pts
}

/// Build a tessellated BRep for the union of two Z-aligned cone frustums with
/// offset XY centers.
///
/// The first cone is assumed to be the one with the larger Z span (the one that
/// provides the "above-overlap" walls).
fn build_cone_cone_union_tessellated(
    c1_xy: DVec2, c1_z_lo: f64, c1_z_hi: f64, c1_r_lo: f64, c1_r_hi: f64,
    c2_xy: DVec2, c2_z_lo: f64, c2_z_hi: f64, c2_r_lo: f64, c2_r_hi: f64,
) -> Option<BRep> {
    let tol = TOLERANCE_LEN_MIN;
    let overlap_lo = c1_z_lo.max(c2_z_lo);
    let overlap_hi = c1_z_hi.min(c2_z_hi);
    let has_overlap = overlap_hi > overlap_lo + tol;
    let has_top = c1_z_hi > overlap_hi + tol;

    // Cone 1 radius at any Z (the tall cone that continues above overlap).
    let c1_dr_dz = (c1_r_hi - c1_r_lo) / (c1_z_hi - c1_z_lo);
    let c1_r = |z: f64| (c1_r_lo + c1_dr_dz * (z - c1_z_lo)).max(0.0);

    // Cone 2 radius at any Z.
    let c2_dr_dz = (c2_r_hi - c2_r_lo) / (c2_z_hi - c2_z_lo);
    let c2_r = |z: f64| (c2_r_lo + c2_dr_dz * (z - c2_z_lo)).max(0.0);

    const N_ARC: usize = 48;   // sample points per circle boundary arc
    const N_SLICE: usize = 32; // Z slices in the overlap region

    let empty_wire = || Wire { edges: vec![] };
    // Simple world transform: points are already in world XY.
    let to_world = |x: f64, y: f64, z: f64| -> DVec3 { DVec3::new(x, y, z) };

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

    // Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì Overlap region [overlap_lo, overlap_hi] (two-circle union) Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì
    if has_overlap {
        // Pre-compute slice polygons at each Z level.
        let n_sl = N_SLICE.max(1);
        let mut polys: Vec<Vec<DVec2>> = Vec::with_capacity(n_sl + 1);
        for i in 0..=n_sl {
            let z = overlap_lo + (overlap_hi - overlap_lo) * i as f64 / n_sl as f64;
            let r_a = c1_r(z);
            let r_b = c2_r(z);
            polys.push(two_circle_union_ccw_pts(c1_xy, r_a, c2_xy, r_b, N_ARC, N_ARC));
        }

        // Extrude walls between consecutive slices.
        for i in 0..n_sl {
            let z0 = overlap_lo + (overlap_hi - overlap_lo) * i as f64 / n_sl as f64;
            let z1 = overlap_lo + (overlap_hi - overlap_lo) * (i + 1) as f64 / n_sl as f64;
            let n = polys[i].len();
            if n < 3 { continue; }

            let mut idx = Vec::with_capacity(2 * n);
            for p in &polys[i] { idx.push(add_v(to_world(p.x, p.y, z0))); }
            for p in &polys[i] { idx.push(add_v(to_world(p.x, p.y, z1))); }
            let mut tris = Vec::with_capacity(n * 2);
            for j in 0..n {
                let k = (j + 1) % n;
                tris.push([idx[j], idx[k], idx[n + k]]);
                tris.push([idx[j], idx[n + k], idx[n + j]]);
            }
            faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
        }

        // Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì Interface face at overlap_hi: two-circle union minus the circle of cone 1 Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì
        if has_top {
            let r_at = c1_r(overlap_hi);
            let poly = two_circle_union_ccw_pts(c1_xy, r_at, c2_xy, c2_r(overlap_hi), N_ARC, N_ARC);
            add_interface_face(
                &mut add_v, &mut faces,
                &poly, overlap_hi, -DVec3::Z,
                c1_xy, r_at,
                &to_world, &empty_wire,
            );
        }
    }

    // Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì Top region [overlap_hi, c1_z_hi] (single cone) Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì
    if has_top {
        let n_sl = N_SLICE.max(1);
        let z0 = overlap_hi;
        let z1 = c1_z_hi;
        // Build wall using the single-circle polygon.
        let r_lo = c1_r(z0);
        let _r_hi = c1_r(z1);
        if r_lo > tol {
            let bot_poly = two_circle_union_ccw_pts(c1_xy, r_lo, c2_xy, c2_r(z0), N_ARC, N_ARC);
            // Remap to same vertex count as interface uses.
            let poly_single = build_circle_polygon(c1_xy.x, c1_xy.y, r_lo);
            let n_pts = bot_poly.len();
            let poly = if poly_single.len() != n_pts {
                let ref_pt = c1_xy + DVec2::new(r_lo, 0.0);
                remap_polygon_arclength(&poly_single, n_pts, ref_pt)
            } else {
                poly_single
            };

            // Bottom cap at overlap_hi if no interface was created.
            if !has_overlap {
                add_cap_face(&mut add_v, &mut faces, &poly, z0, -DVec3::Z, &to_world, &empty_wire);
            }

            // Build wall from z0 to z1.
            let dz = (z1 - z0) / n_sl as f64;
            // We need to handle the transition here: use poly at z0 and circle at z1.
            // But since both use the same vertex count (after remapping), the wall
            // connects them directly.
            let top_r_at = |z: f64| {
                let r = c1_r(z);
                let mut p = build_circle_polygon(c1_xy.x, c1_xy.y, r);
                if p.len() != n_pts {
                    let ref_pt = c1_xy + DVec2::new(r, 0.0);
                    p = remap_polygon_arclength(&p, n_pts, ref_pt);
                }
                p
            };

            for i in 0..n_sl {
                let za = z0 + dz * i as f64;
                let zb = z0 + dz * (i + 1) as f64;
                let pts_a = top_r_at(za);
                let pts_b = top_r_at(zb);
                let n = pts_a.len();
                if n < 3 { continue; }
                let mut idx = Vec::with_capacity(2 * n);
                for p in &pts_a { idx.push(add_v(to_world(p.x, p.y, za))); }
                for p in &pts_b { idx.push(add_v(to_world(p.x, p.y, zb))); }
                let mut tris = Vec::with_capacity(n * 2);
                for j in 0..n {
                    let k = (j + 1) % n;
                    tris.push([idx[j], idx[k], idx[n + k]]);
                    tris.push([idx[j], idx[n + k], idx[n + j]]);
                }
                faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
            }

            // Top cap at c1_z_hi.
            let r_top = c1_r(z1);
            if r_top > tol {
                let mut top_poly = build_circle_polygon(c1_xy.x, c1_xy.y, r_top);
                if top_poly.len() != n_pts {
                    let ref_pt = c1_xy + DVec2::new(r_top, 0.0);
                    top_poly = remap_polygon_arclength(&top_poly, n_pts, ref_pt);
                }
                add_cap_face(&mut add_v, &mut faces, &top_poly, z1, DVec3::Z, &to_world, &empty_wire);
            }
        }
    }

    // Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì Bottom cap Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì
    if !has_overlap || c1_z_lo < overlap_lo - tol {
        let z = c1_z_lo;
        let r = c1_r(z);
        if r > tol {
            let poly = two_circle_union_ccw_pts(c1_xy, r, c2_xy, c2_r(z), N_ARC, N_ARC);
            add_cap_face(&mut add_v, &mut faces, &poly, z, -DVec3::Z, &to_world, &empty_wire);
        }
    } else {
        let z = overlap_lo;
        let r_a = c1_r(z);
        let r_b = c2_r(z);
        if r_a > tol || r_b > tol {
            let poly = two_circle_union_ccw_pts(c1_xy, r_a, c2_xy, r_b, N_ARC, N_ARC);
            add_cap_face(&mut add_v, &mut faces, &poly, z, -DVec3::Z, &to_world, &empty_wire);
        }
    }

    // Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì Top cap (only if no top region was built above) Èñ≥ÂÖâÂÅìÈñ≥ÂÖâÂÅì
    if !has_top {
        let z = c1_z_hi;
        let r = c1_r(z);
        if r > tol {
            let poly = if c1_z_hi > overlap_hi + tol || !has_overlap {
                build_circle_polygon(c1_xy.x, c1_xy.y, r)
            } else {
                two_circle_union_ccw_pts(c1_xy, r, c2_xy, c2_r(z), N_ARC, N_ARC)
            };
            add_cap_face(&mut add_v, &mut faces, &poly, z, DVec3::Z, &to_world, &empty_wire);
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
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}

/// Detect cone-cone union with offset centers (different XY, same Z-axis).
fn try_union_offset_cones_one_dir(a: &BRep, b: &BRep) -> Option<BRep> {
    let (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi) = detect_z_axis_cone_frustum(a)?;
    let (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi) = detect_z_axis_cone_frustum(b)?;

    // Ensure c1 is the one with larger Z span (it continues above overlap).
    let (s_a, s_b) = if (c1_z_hi - c1_z_lo) >= (c2_z_hi - c2_z_lo) {
        ((c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi),
         (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi))
    } else {
        ((c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi),
         (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi))
    };

    let (c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi) = s_a;
    let (c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi) = s_b;

    let tol = TOLERANCE_LEN_MIN;

    // Not offset Èñ?delegate to coaxial path.
    if (c1_xy.x - c2_xy.x).abs() < tol && (c1_xy.y - c2_xy.y).abs() < tol {
        return None;
    }

    // Z ranges must overlap.
    if c1_z_hi < c2_z_lo - tol || c2_z_hi < c1_z_lo - tol {
        return None;
    }

    build_cone_cone_union_tessellated(
        c1_xy, c1_z_lo, c1_z_hi, c1_r_lo, c1_r_hi,
        c2_xy, c2_z_lo, c2_z_hi, c2_r_lo, c2_r_hi,
    )
}

/// Bidirectional wrapper for offset cone-cone union.
pub fn try_union_offset_cones(a: &BRep, b: &BRep) -> Option<BRep> {
    try_union_offset_cones_one_dir(a, b)
        .or_else(|| try_union_offset_cones_one_dir(b, a))
}

/// Build BRep for outer conical frustum minus inner conical frustum (coaxial, Z-aligned).
/// Result is a hollow conical frustum: outer lateral + bottom cap + top annulus + inner lateral + cavity floor.
///
/// All faces are triangulated (no analytic surfaces) in a single shell.
fn build_conical_frustum_minus_frustum_brep(
    zo_lo: f64, zo_hi: f64, ro_lo: f64, ro_hi: f64,
    zi_lo: f64, zi_hi: f64, ri_lo: f64, ri_hi: f64,
) -> Option<BRep> {
    use std::f64::consts::TAU;
    const N: usize = 48; // Circumferential divisions

    let empty_wire = || Wire { edges: vec![] };
    let tol = TOLERANCE_LEN_MIN;

    // Overlap Z range (inner clamped to outer)
    let z_olap_lo = zi_lo.max(zo_lo);
    let z_olap_hi = zi_hi.min(zo_hi);
    let has_overlap = z_olap_hi > z_olap_lo + tol;

    // Inner cone radius at overlap boundaries
    let dri = (ri_hi - ri_lo) / (zi_hi - zi_lo);
    let ri_at_olap_lo = ri_lo + dri * (z_olap_lo - zi_lo);
    let ri_at_olap_hi = ri_lo + dri * (z_olap_hi - zi_lo);

    let ring_pts = |z: f64, r: f64| -> Vec<DVec3> {
        (0..N).map(|i| {
            let ang = TAU * i as f64 / N as f64;
            let (c, s) = ang.sin_cos();
            DVec3::new(r * c, r * s, z)
        }).collect()
    };

    let outer_bot = ring_pts(zo_lo, ro_lo);
    let outer_top = ring_pts(zo_hi, ro_hi);

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

    // 1. Outer lateral (z=zo_lo to z=zo_hi)
    {
        let mut tris = Vec::new();
        wall_grid(&mut add_v, &outer_bot, &outer_top, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    }

    // 2. Bottom at z=zo_lo
    //    - If inner starts at or below outer bottom: annulus
    //    - Otherwise: full disk
    if zi_lo <= zo_lo + tol && has_overlap && ri_at_olap_lo < ro_lo - tol {
        // Hole goes through the bottom
        let mut tris = Vec::new();
        annulus_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_lo), ri_at_olap_lo, ro_lo, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    } else {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_lo), ro_lo, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    }

    // 3. Top at z=zo_hi
    //    - If inner ends at or above outer top: annulus
    //    - Otherwise: full disk
    if zi_hi >= zo_hi - tol && has_overlap && ri_at_olap_hi < ro_hi - tol {
        let mut tris = Vec::new();
        annulus_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_hi), ri_at_olap_hi, ro_hi, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    } else {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, zo_hi), ro_hi, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    }

    // 4. Inner lateral Èñ?cavity wall (overlap region only)
    if has_overlap {
        let inner_bot = ring_pts(z_olap_lo, ri_at_olap_lo);
        let inner_top = ring_pts(z_olap_hi, ri_at_olap_hi);
        let mut tris = Vec::new();
        wall_grid(&mut add_v, &inner_bot, &inner_top, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    }

    // 5. Cavity floor Èñ?where hole starts (if inner starts above outer bottom)
    if zi_lo > zo_lo + tol && has_overlap && ri_at_olap_lo > tol {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, z_olap_lo), ri_at_olap_lo, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
    }

    // 6. Cavity ceiling Èñ?where hole ends (if inner ends below outer top)
    if zi_hi < zo_hi - tol && has_overlap && ri_at_olap_hi > tol {
        let mut tris = Vec::new();
        disk_tri_fan(&mut add_v, DVec3::new(0.0, 0.0, z_olap_hi), ri_at_olap_hi, N, &mut tris);
        faces.push(Face { outer_wire: empty_wire(), inner_wires: vec![], normal: DVec3::ZERO, triangles: tris, sample_point: None, mesh_dirty: false, surface_idx: None });
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
        edge_same_parameter: vec![], edge_same_range: vec![], edge_vertex_params: vec![]};

    Some(BRep {
        vertices: verts, edges: vec![],
        solids: vec![Solid { shells: vec![Shell { faces }] }],
        geom, compound: None, compsolid: None,
    })
}
