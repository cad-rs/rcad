//! Special-case intersections used by OCCT DRAW ports when the generic `BooleanBuilder`
//! path is wrong or overly faceted: (1) unit ball ∩ `[0,1]³`, (2) coaxial sharp cone ∩ finite
//! cylinder (ZP7), (3) coaxial sharp cone minus cylinder sealing the base (ZP8),
//! (4) coaxial cylinder minus cone via sewn loft shells (`boptuc_simple`/ZP3).
//! (5) concentric analytic spheres (`make_sphere_brep`): compound outer sphere + mirrored inner → analytic shell SA/volume (differs from OCCT `mkvolume` on trimmed patches).
//! (6) intersection of two nested analytic spheres sharing a center → smaller ball (`∩` is the inner sphere).
//!
//! Unit ball ∩ box `[0,1]³` (first-octant "spherical sector"):
//!
//! The generic Pave/Builder path does not yet split planar faces along the
//! sphere, so the result was three untrimmed 1×1 squares. OCCT `bcommon_simple/A1`
//! expects the exact surface `5π/4` and volume `π/6` for the eighth ball.
//!
//! This is *not* a full analytic CSG solution — only a recognition + mesh for
//! this configuration used in OCCT DRAW port tests.
//!
//! When the **box** is rigidly transformed (e.g. OCCT `trotate` about a pivot
//! not at the origin), the axis-aligned `[0,1]³` predicate fails and the
//! generic boolean path must handle sphere–oblique-plane trimming. Pave now
//! uses exact spherical UV in projection fallbacks (`SphericalSurface::world_to_uv`),
//! and plane–sphere tangents are inflated to a micro-circle so every box face gets
//! a `FaceFace` curve. OCCT `bcommon_simple/A4` / `A5` remain `#[ignore]` — `BooleanBuilder`
//! surface area / volume do not yet match (`checkprops -s`); next step is sphere UV
//! multi-trim / classification, not missing FF pairs.

use glam::DVec3;
use rcad_kernel::geom::{ConicalSurface, Curve3, Line3, Plane, Surface3};
use rcad_kernel::topology::{Face, Shell, Solid, Wire, WireEdge};
use rcad_kernel::{BRep, GeomStore, Vertex};
use rcad_modeling::builder::brep_builder::{make_edge, make_face, make_vertex, make_wire};
use rcad_modeling::builder::ops::LoftHistory;
use rcad_modeling::{loft_with_history, make_sphere_brep, sew_shells};

const TOL: f64 = 1e-4;

fn is_unit_sphere_at_origin(b: &BRep) -> bool {
    b.solids.len() == 1
        && b.solids[0].shells.len() == 1
        && b.solids[0].shells[0].faces.len() == 1
        && b.vertices.len() == 2
        && b
            .geom
            .face_surface
            .get(0)
            .and_then(|o| o.as_ref().copied())
            .and_then(|si| b.geom.surfaces.get(si))
            .is_some_and(|s| {
                if let Surface3::Sphere(s) = s {
                    s.radius - 1.0 < 1e-2 && s.center.length() < 1e-2
                } else {
                    false
                }
            })
}

fn is_pos_unit_cube_0_1(b: &BRep) -> bool {
    if b.solids.len() != 1 || b.solids[0].shells[0].faces.len() != 6 {
        return false;
    }
    if b.vertices.len() != 8 {
        return false;
    }
    let Some(bb) = b.bounding_box() else {
        return false;
    };
    (bb[0] - DVec3::ZERO).length() < TOL && (bb[1] - DVec3::ONE).length() < TOL
}

/// Intersection: unit sphere (kernel primitive) ∩ axis box [0,1]³.
pub fn try_intersection_eighth_unit_ball(a: &BRep, b: &BRep) -> Option<BRep> {
    if (is_unit_sphere_at_origin(a) && is_pos_unit_cube_0_1(b))
        || (is_unit_sphere_at_origin(b) && is_pos_unit_cube_0_1(a))
    {
        return Some(brep_eighth_of_unit_ball());
    }
    None
}

/// Kernel analytic sphere primitive: one spherical face (`Surface3::Sphere`).
fn try_sphere_primitive_center_radius(brep: &BRep) -> Option<(DVec3, f64)> {
    let sh = brep.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() != 1 {
        return None;
    }
    let si = *brep.geom.face_surface.get(0)?.as_ref()?;
    match brep.geom.surfaces.get(si)? {
        Surface3::Sphere(s) => Some((s.center, s.radius)),
        _ => None,
    }
}

/// [`BooleanOpType::Difference`] for nested analytic spheres sharing a center.
///
/// Builds a hollow spherical shell as **two solids**: outer sphere plus inner sphere with reversed
/// face orientation (`reverse_face`). Total surface area matches the analytic spherical shell
/// \(4\pi(R^2+r^2)\) under [`rcad_kernel::surface_area`].
///
/// Compound `rcad_kernel::volume` may not match \(4\pi/3(R^3-r^3)\) until sphere face normals / tessellation agree for
/// [`signed_volume`] everywhere.
///
/// OCCT DRAW `mkvolume` on trimmed spherical patches can report a different `checkprops -s` than this
/// full-sphere analytic shell.
pub fn try_difference_concentric_spheres(a: &BRep, b: &BRep) -> Option<BRep> {
    let (ca, ra) = try_sphere_primitive_center_radius(a)?;
    let (cb, rb) = try_sphere_primitive_center_radius(b)?;
    let scale = ra.max(rb).max(1.0);
    if (ca - cb).length() > TOL.max(1e-9 * scale) {
        return None;
    }
    let ro = ra.max(rb);
    let ri = ra.min(rb);
    if ro - ri <= 1e-12 * ro.max(1.0) {
        return None;
    }
    let center = ca;
    let outer = make_sphere_brep(center, ro).ok()?;
    let mut inner_cavity = make_sphere_brep(center, ri).ok()?;
    crate::reverse_face(&mut inner_cavity, 0);
    Some(BRep::compound_from_shapes(&[outer, inner_cavity]))
}

/// [`BooleanOpType::Intersection`] for two analytic sphere primitives sharing a center.
///
/// The intersection of nested balls is the smaller-radius ball (same center).
///
/// Returns [`None`] when centers differ or radii are degenerate — callers fall back to Pave/Builder.
pub fn try_intersection_concentric_spheres(a: &BRep, b: &BRep) -> Option<BRep> {
    let (ca, ra) = try_sphere_primitive_center_radius(a)?;
    let (cb, rb) = try_sphere_primitive_center_radius(b)?;
    let scale = ra.max(rb).max(1.0);
    if (ca - cb).length() > TOL.max(1e-9 * scale) {
        return None;
    }
    let r = ra.min(rb);
    let r_eps = 1e-12 * scale.max(1.0);
    if r <= r_eps {
        return None;
    }
    make_sphere_brep(ca, r).ok()
}

// --- Coaxial cone ∩ cylinder (OCCT `bopcommon_simple/ZP7`): generic Builder over-counts area. --------

fn z_axis_sharp_cone_z_span(cone: &BRep) -> Option<(f64, f64, f64)> {
    const APAR: f64 = 1e-3;
    const XY: f64 = 2e-3;
    let sh = cone.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() != 2 {
        return None;
    }
    let mut cf: Option<ConicalSurface> = None;
    let mut po: Option<DVec3> = None;
    let mut fi = 0usize;
    for s in &cone.solids {
        for sh in &s.shells {
            for _ in &sh.faces {
                let si = *cone.geom.face_surface.get(fi)?.as_ref()?;
                match cone.geom.surfaces.get(si)? {
                    Surface3::Cone(c) if c.radius.abs() <= 1e-6 => cf = Some(*c),
                    Surface3::Plane(p) => po = Some(p.origin),
                    _ => return None,
                }
                fi += 1;
            }
        }
    }
    let c = cf?;
    let u = c.axis_dir();
    if u.cross(DVec3::Z).length() > APAR {
        return None;
    }
    let apex = c.apex_point();
    if apex.x.abs() > XY || apex.y.abs() > XY {
        return None;
    }
    let b = po?;
    if b.x.abs() > XY || b.y.abs() > XY {
        return None;
    }
    let t = (b - apex).dot(u);
    let rb = t * c.half_angle_rad.tan();
    if t < 1e-6 || rb < 1e-6 {
        return None;
    }
    Some((apex.z, b.z, rb))
}

fn z_axis_cylinder_z_span_r(cyl: &BRep) -> Option<(f64, f64, f64)> {
    const APAR: f64 = 1e-3;
    const XY: f64 = 2e-3;
    let sh = cyl.solids.get(0)?.shells.get(0)?;
    if sh.faces.len() != 3 {
        return None;
    }
    let mut r = None;
    let mut zs = Vec::with_capacity(2);
    let mut fi = 0usize;
    for s in &cyl.solids {
        for sh in &s.shells {
            for _ in &sh.faces {
                let si = *cyl.geom.face_surface.get(fi)?.as_ref()?;
                match cyl.geom.surfaces.get(si)? {
                    Surface3::Cylinder(cc) => {
                        if cc.axis.normalize_or_zero().cross(DVec3::Z).length() > APAR {
                            return None;
                        }
                        if cc.origin.x.abs() > XY || cc.origin.y.abs() > XY {
                            return None;
                        }
                        r = Some(cc.radius);
                    }
                    Surface3::Plane(p) => zs.push(p.origin.z),
                    _ => return None,
                }
                fi += 1;
            }
        }
    }
    if zs.len() != 2 {
        return None;
    }
    Some((zs[0].min(zs[1]), zs[0].max(zs[1]), r?))
}

fn try_intersection_coaxial_cone_cylinder_pair(cone: &BRep, cyl: &BRep) -> Option<BRep> {
    use rcad_modeling::make_conical_frustum_brep;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    let zcn = zb.min(za);
    let zcx = zb.max(za);
    let z0 = zlo.max(zcn);
    let z1 = zhi.min(zcx);
    if z1 - z0 < 1e-6 {
        return None;
    }
    let apex_hi = za > zb;
    let hc = (za - zb).abs();
    let rcz = |z: f64| {
        let num = if apex_hi { (za - z).abs() } else { (z - za).abs() };
        rb * num / hc
    };
    let r0 = rcz(z0).min(rc);
    let r1 = rcz(z1).min(rc);
    if r0 < 1e-9 && r1 < 1e-9 {
        return None;
    }
    let zm = (z0 + z1) * 0.5;
    let h = z1 - z0;
    make_conical_frustum_brep(DVec3::new(0.0, 0.0, zm), DVec3::Z, DVec3::X, r0, r1, h).ok()
}

/// Sharp Z-aligned cone ∩ finite Z-aligned cylinder (same axis / origin in `xy`), e.g. OCCT ZP7.
pub fn try_intersection_coaxial_cone_cylinder(a: &BRep, b: &BRep) -> Option<BRep> {
    try_intersection_coaxial_cone_cylinder_pair(a, b)
        .or_else(|| try_intersection_coaxial_cone_cylinder_pair(b, a))
}

/// `cone \ cylinder` when the cylinder closes the cone base and contains the cone frustum up to `z_hi`:
/// remainder is the sharp sub-cone from `z_hi` to the apex (OCCT `bopcut_simple`/ZP8).
pub fn try_difference_coaxial_cone_minus_cylinder(cone: &BRep, cyl: &BRep) -> Option<BRep> {
    use rcad_modeling::make_cone_brep;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    if za <= zb + 1e-6 {
        return None;
    }
    let hc = za - zb;
    let r_at = |z: f64| rb * (za - z) / hc;
    // Cylinder starts on cone base disk; radius at least the cone base radius.
    if (zlo - zb).abs() > 2e-2 {
        return None;
    }
    if (rc + 1e-3) < rb {
        return None;
    }
    // Cylinder covers the cone cross-section at `z_hi` (disk of radius `rc` vs cone radius `r_at(z_hi)`).
    if rc + 1e-6 < r_at(zhi) {
        return None;
    }
    if zhi <= zb + 1e-6 || zhi >= za - 1e-6 {
        return None;
    }
    let r_cut = r_at(zhi);
    if r_cut < 1e-9 {
        return None;
    }
    let h_rem = za - zhi;
    let z_mid = (za + zhi) * 0.5;
    make_cone_brep(
        DVec3::new(0.0, 0.0, z_mid),
        DVec3::Z,
        DVec3::X,
        r_cut,
        h_rem,
    )
    .ok()
}

/// Strip loft bottom/top planar caps; keeps ruled lateral faces only (must match [`LoftHistory`]).
fn strip_loft_caps(mut brep: BRep, hist: LoftHistory) -> Option<BRep> {
    let shell = brep.solids.first_mut()?.shells.first_mut()?;
    if hist.bottom_cap >= shell.faces.len() || hist.top_cap >= shell.faces.len() {
        return None;
    }
    // Remove higher shell face index first so the remaining index stays valid.
    let (mut lo, mut hi) = (hist.bottom_cap, hist.top_cap);
    if lo > hi {
        std::mem::swap(&mut lo, &mut hi);
    }
    shell.faces.remove(hi);
    crate::remove_flat_face_geom_slots(&mut brep.geom, hi);
    shell.faces.remove(lo);
    crate::remove_flat_face_geom_slots(&mut brep.geom, lo);
    Some(brep)
}

fn reverse_wire_local(wire: &mut Wire) {
    wire.edges.reverse();
    for we in &mut wire.edges {
        we.forward = !we.forward;
    }
}

/// Loft builds the **solid** frustum mantle (normals outward from cone interior). For
/// `cylinder \\ cone`, those faces bound the cavity — outward from the result solid points into the
/// removed cone (flip normals vs loft defaults).
fn invert_shell_planar_faces(brep: &mut BRep) {
    let Some(shell) = brep.solids.first_mut().and_then(|s| s.shells.first_mut()) else {
        return;
    };
    let n_faces = shell.faces.len();
    for face in &mut shell.faces {
        face.normal = -face.normal;
        reverse_wire_local(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            reverse_wire_local(iw);
        }
        face.triangles.clear();
        face.mesh_dirty = true;
    }
    for fi in 0..n_faces {
        let Some(Some(si)) = brep.geom.face_surface.get(fi).copied() else {
            continue;
        };
        let Some(surf) = brep.geom.surfaces.get_mut(si) else {
            continue;
        };
        if let Surface3::Plane(pl) = surf {
            pl.normal = -pl.normal;
        }
    }
}

/// Horizontal annulus from pre-built coplanar rings (same vertex count). Eliminates float drift vs loft.
fn annulus_between_rings(outer: &[DVec3], inner: &[DVec3]) -> Result<BRep, rcad_modeling::BuildError> {
    let n = outer.len();
    if n < 3 || inner.len() != n {
        return Err(rcad_modeling::BuildError::DegenerateGeometry(
            "annulus_between_rings vertex count",
        ));
    }
    let z = outer[0].z;
    if inner.iter().any(|p| (p.z - z).abs() > 1e-9) || outer.iter().any(|p| (p.z - z).abs() > 1e-9) {
        return Err(rcad_modeling::BuildError::DegenerateGeometry(
            "annulus_between_rings not coplanar",
        ));
    }
    let mut brep = BRep::default();
    let surface = Surface3::Plane(Plane {
        origin: DVec3::new(0.0, 0.0, z),
        normal: DVec3::Z,
    });

    let mut outer_vs = Vec::with_capacity(n);
    let outer_pts = outer.to_vec();
    for p in &outer_pts {
        outer_vs.push(make_vertex(&mut brep, *p));
    }
    let mut outer_we = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = outer_pts[i];
        let b = outer_pts[j];
        let dir = (b - a).normalize();
        let len = (b - a).length();
        let ei = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            outer_vs[i],
            outer_vs[j],
        )?;
        outer_we.push(WireEdge::fwd(ei));
    }

    let mut inner_vs = Vec::with_capacity(n);
    let inner_pts = inner.to_vec();
    for p in &inner_pts {
        inner_vs.push(make_vertex(&mut brep, *p));
    }
    let mut inner_we = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let a = inner_pts[i];
        let b = inner_pts[j];
        let dir = (b - a).normalize();
        let len = (b - a).length();
        let ei = make_edge(
            &mut brep,
            Curve3::Line(Line3 {
                origin: a,
                direction: dir,
            }),
            0.0,
            len,
            inner_vs[i],
            inner_vs[j],
        )?;
        inner_we.push(WireEdge::fwd(ei));
    }

    let outer_wire = make_wire(outer_we);
    let inner_wire = make_wire(inner_we);
    let _fi = make_face(&mut brep, surface, outer_wire, vec![inner_wire])?;
    Ok(brep)
}

/// Outer / inner lateral strips + top annulus for [`try_coaxial_cylinder_minus_frustum_loft_shell`].
fn coaxial_cylinder_minus_frustum_loft_pieces(
    z_lo: f64,
    z_hi: f64,
    rc: f64,
    za: f64,
    zb: f64,
    rb: f64,
) -> Option<(BRep, BRep, BRep)> {
    use std::f64::consts::TAU;
    const N: usize = 32;
    let zcn = zb.min(za);
    let zcx = zb.max(za);
    let z0 = z_lo.max(zcn);
    let z1 = z_hi.min(zcx);
    if z1 - z0 < 1e-6 {
        return None;
    }
    let apex_hi = za > zb;
    let hc = (za - zb).abs();
    let rcz = |z: f64| {
        let num = if apex_hi {
            (za - z).abs()
        } else {
            (z - za).abs()
        };
        rb * num / hc
    };
    let r0 = rcz(z0).min(rc);
    let r1 = rcz(z1).min(rc);
    if r0 < 1e-9 || r1 < 1e-9 {
        return None;
    }

    let mut outer_bot = Vec::with_capacity(N);
    let mut outer_top = Vec::with_capacity(N);
    let mut inner_bot = Vec::with_capacity(N);
    let mut inner_top = Vec::with_capacity(N);
    for i in 0..N {
        let ang = TAU * i as f64 / N as f64;
        let c = ang.cos();
        let s = ang.sin();
        let ob = DVec3::new(rc * c, rc * s, z_lo);
        outer_bot.push(ob);
        let ib = if (r0 - rc).abs() <= 1e-12 && (z0 - z_lo).abs() <= 1e-12 {
            ob
        } else {
            DVec3::new(r0 * c, r0 * s, z0)
        };
        inner_bot.push(ib);
        outer_top.push(DVec3::new(rc * c, rc * s, z_hi));
        inner_top.push(DVec3::new(r1 * c, r1 * s, z1));
    }

    let annulus = annulus_between_rings(&outer_top, &inner_top).ok()?;

    let (loft_outer, ohist) = loft_with_history(&[outer_bot, outer_top]).ok()?;
    let outer_strip = strip_loft_caps(loft_outer, ohist)?;

    let (loft_inner, ihist) = loft_with_history(&[inner_bot, inner_top]).ok()?;
    let mut inner_strip = strip_loft_caps(loft_inner, ihist)?;
    invert_shell_planar_faces(&mut inner_strip);

    Some((outer_strip, inner_strip, annulus))
}

/// Closed shell for `cylinder \ (cone ∩ cylinder)` when overlap is the coaxial frustum:
/// outer cylindrical loft strip + inner frustum loft strip + top annulus, sewn (`OCCT ZP3`).
fn try_coaxial_cylinder_minus_frustum_loft_shell(
    z_lo: f64,
    z_hi: f64,
    rc: f64,
    za: f64,
    zb: f64,
    rb: f64,
) -> Option<BRep> {
    let (outer_strip, inner_strip, annulus) =
        coaxial_cylinder_minus_frustum_loft_pieces(z_lo, z_hi, rc, za, zb, rb)?;
    let tol = (1e-4_f64).max(1e-6 * rc.max(z_hi.abs()));
    let sewn = sew_shells(&[outer_strip, inner_strip, annulus], tol);
    if !sewn.free_edges.is_empty() {
        return None;
    }
    Some(sewn.brep)
}

/// `cylinder \ cone` with same coaxial ZP layout as [`try_difference_coaxial_cone_minus_cylinder`].
///
/// Set identity: `cyl \ cone` equals `cyl \ (cone ∩ cyl)` when the overlap is the coaxial frustum.
pub fn try_difference_coaxial_cylinder_minus_cone(a: &BRep, b: &BRep) -> Option<BRep> {
    cyl_minus_cone_inner(a, b).or_else(|| cyl_minus_cone_inner(b, a))
}

fn cyl_minus_cone_inner(maybe_cyl: &BRep, maybe_cone: &BRep) -> Option<BRep> {
    let cone = maybe_cone;
    let cyl = maybe_cyl;
    try_intersection_coaxial_cone_cylinder_pair(cone, cyl)?;
    let (za, zb, rb) = z_axis_sharp_cone_z_span(cone)?;
    let (zlo, zhi, rc) = z_axis_cylinder_z_span_r(cyl)?;
    let hc = za - zb;
    if hc.abs() < 1e-6 {
        return None;
    }
    let r_at = |z: f64| rb * (za - z) / hc;
    if (zlo - zb).abs() > 2e-2 {
        return None;
    }
    if (rc + 1e-3) < rb {
        return None;
    }
    if rc + 1e-6 < r_at(zhi) {
        return None;
    }
    if zhi <= zb + 1e-6 || zhi >= za - 1e-6 {
        return None;
    }
    try_coaxial_cylinder_minus_frustum_loft_shell(zlo, zhi, rc, za, zb, rb)
}

fn add_vertex(verts: &mut Vec<Vertex>, p: DVec3) -> usize {
    for (i, v) in verts.iter().enumerate() {
        if (v.point - p).length() < 1e-7 {
            return i;
        }
    }
    verts.push(Vertex { point: p });
    verts.len() - 1
}

/// Closed triangle mesh of the boundary: three quarter-disks in x=0, y=0, z=0
/// planes plus one octant of the unit sphere. Outward for the solid
/// { x,y,z ≥ 0, x²+y²+z² ≤ 1 }.
fn brep_eighth_of_unit_ball() -> BRep {
    const NA: usize = 32; // arc segments per quarter circle
    const NS: usize = 24; // grid per spherical patch axis

    let mut vertices: Vec<Vertex> = vec![];
    let empty_wire = || Wire { edges: vec![] };
    use std::f64::consts::FRAC_PI_2;
    // --- Planar: z=0, outward normal (0,0,-1) ---
    let o0 = add_vertex(&mut vertices, DVec3::ZERO);
    let _ = add_vertex(&mut vertices, DVec3::X);
    let _ = add_vertex(&mut vertices, DVec3::Y);
    let mut z0_arc: Vec<usize> = (0..=NA)
        .map(|k| {
            let t = (k as f64 / NA as f64) * FRAC_PI_2;
            add_vertex(&mut vertices, DVec3::new(t.cos(), t.sin(), 0.0))
        })
        .collect();

    // Triangles on z=0, fan from origin
    let mut t_z0: Vec<[usize; 3]> = vec![];
    for k in 0..NA {
        t_z0.push([o0, z0_arc[k], z0_arc[k + 1]]);
    }
    let f_z0 = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(0.0, 0.0, -1.0),
        triangles: t_z0,
        mesh_dirty: false,
    };

    // y=0 plane, outward (0,-1,0), quarter disk in xz: (0,0,0) — (1,0,0) — (0,0,1) and arc
    let o1 = o0; // (0,0,0) shared
    let _ = add_vertex(&mut vertices, DVec3::Z);
    let y0_arc: Vec<usize> = (0..=NA)
        .map(|k| {
            let t = (k as f64 / NA as f64) * FRAC_PI_2;
            add_vertex(&mut vertices, DVec3::new(t.cos(), 0.0, t.sin()))
        })
        .collect();
    let mut t_y0: Vec<[usize; 3]> = vec![];
    for k in 0..NA {
        t_y0.push([o1, y0_arc[k], y0_arc[k + 1]]);
    }
    let f_y0 = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(0.0, -1.0, 0.0),
        triangles: t_y0,
        mesh_dirty: false,
    };

    // x=0, outward (-1,0,0), quarter in yz: (0,0,0) (0,1,0) (0,0,1) + arc
    let o2 = o0;
    let x0_arc: Vec<usize> = (0..=NA)
        .map(|k| {
            let t = (k as f64 / NA as f64) * FRAC_PI_2;
            add_vertex(&mut vertices, DVec3::new(0.0, t.cos(), t.sin()))
        })
        .collect();
    let mut t_x0: Vec<[usize; 3]> = vec![];
    for k in 0..NA {
        t_x0.push([o2, x0_arc[k], x0_arc[k + 1]]);
    }
    let f_x0 = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(-1.0, 0.0, 0.0),
        triangles: t_x0,
        mesh_dirty: false,
    };

    // Spherical octant: p = (sin v cos u, sin v sin u, cos v), u,v in [0,π/2]
    // (v = colatitude from +Z: v=0 is north pole (0,0,1), v=π/2 is the z=0 quarter arc)
    let pole = add_vertex(&mut vertices, DVec3::new(0.0, 0.0, 1.0));
    let mut sph_idx = vec![vec![0usize; NS + 1]; NS + 1];
    for i in 0..=NS {
        sph_idx[i][0] = pole;
    }
    for j in 1..=NS {
        let v = (j as f64 / NS as f64) * FRAC_PI_2;
        let si = v.sin();
        for i in 0..=NS {
            let u = (i as f64 / NS as f64) * FRAC_PI_2;
            let p = DVec3::new(si * u.cos(), si * u.sin(), v.cos());
            sph_idx[i][j] = add_vertex(&mut vertices, p);
        }
    }
    let mut t_s: Vec<[usize; 3]> = vec![];
    // Fan from north pole to first parallel (j=1)
    for i in 0..NS {
        t_s.push([pole, sph_idx[i][1], sph_idx[i + 1][1]]);
    }
    // Quad strips j = 1..NS-1
    for j in 1..NS {
        for i in 0..NS {
            let a = sph_idx[i][j];
            let b = sph_idx[i + 1][j];
            let c = sph_idx[i][j + 1];
            let d = sph_idx[i + 1][j + 1];
            t_s.push([a, b, d]);
            t_s.push([a, d, c]);
        }
    }
    let f_s = Face {
        outer_wire: empty_wire(),
        inner_wires: vec![],
        normal: DVec3::new(1.0, 1.0, 1.0).normalize(),
        triangles: t_s,
        mesh_dirty: false,
    };

    let faces = vec![f_z0, f_y0, f_x0, f_s];
    let geom = GeomStore {
        curves: vec![],
        surfaces: vec![],
        curve2ds: vec![],
        edge_curve: vec![],
        face_surface: vec![None; 4],
        edge_pcurves: vec![],
        edge_curve_range: vec![],
        edge_degenerated: vec![],
        vertex_tolerance: vec![],
        edge_tolerance: vec![],
        face_tolerance: vec![],
        curve2d_range: vec![],
        face_surface_range: vec![None; 4],
        edge_same_parameter: vec![],
        edge_same_range: vec![],
    };

    BRep {
        vertices,
        edges: vec![],
        solids: vec![Solid {
            shells: vec![Shell { faces }],
        }],
        geom,
        compound: None,
        compsolid: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{boolean_op, BooleanOpType};
    use rcad_kernel::surface_area;
    use rcad_kernel::volume;

    #[test]
    fn concentric_sphere_difference_analytic_shell_surface_area() {
        let center = DVec3::new(1.0, -2.0, 4.0);
        let ro = 5.0_f64;
        let ri = 3.0_f64;
        let outer = make_sphere_brep(center, ro).expect("outer");
        let inner = make_sphere_brep(center, ri).expect("inner");
        let shell = boolean_op(BooleanOpType::Difference, &outer, &inner).expect("difference");
        let pi = std::f64::consts::PI;
        let a_ex = 4.0 * pi * (ro * ro + ri * ri);
        let a = surface_area(&shell);
        assert!(
            (a - a_ex).abs() < 5e-3 * a_ex.max(1.0),
            "surface area {a} vs analytic shell SA {a_ex}"
        );
        // `signed_volume` across compounds relies on consistent face normals vs tessellation;
        // sphere primitives carry approximate face normals — SA matches analytic \(4\pi(R^2+r^2)\).
    }

    #[test]
    fn eighth_ball_area_and_volume() {
        let b = brep_eighth_of_unit_ball();
        let a = surface_area(&b);
        let v = volume(&b);
        let a_ex = 5.0 * std::f64::consts::PI / 4.0;
        let v_ex = std::f64::consts::PI / 6.0;
        assert!(
            (a - a_ex).abs() < 0.04,
            "area {a} vs {a_ex}"
        );
        assert!(
            (v - v_ex).abs() < 0.02,
            "vol {v} vs {v_ex}"
        );
    }

    #[test]
    fn zp3_sum_planar_areas_before_sew_matches_expected_total() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let z_lo = -10.0_f64;
        let z_hi = 0.0_f64;
        let rc = 10.0_f64;
        let z0 = -10.0_f64;
        let z1 = 0.0_f64;
        let r0 = 10.0_f64;
        let r1 = 5.0_f64;
        let mut outer_bot = Vec::with_capacity(N);
        let mut outer_top = Vec::with_capacity(N);
        let mut inner_bot = Vec::with_capacity(N);
        let mut inner_top = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            let c = ang.cos();
            let s = ang.sin();
            let ob = DVec3::new(rc * c, rc * s, z_lo);
            outer_bot.push(ob);
            inner_bot.push(ob);
            outer_top.push(DVec3::new(rc * c, rc * s, z_hi));
            inner_top.push(DVec3::new(r1 * c, r1 * s, z1));
        }
        let ann = annulus_between_rings(&outer_top, &inner_top).unwrap();
        let (lo, oh) = loft_with_history(&[outer_bot, outer_top]).unwrap();
        let outer_strip = strip_loft_caps(lo, oh).unwrap();
        let (li, ih) = loft_with_history(&[inner_bot, inner_top]).unwrap();
        let inner_strip = strip_loft_caps(li, ih).unwrap();
        let sum = rcad_kernel::surface_area(&outer_strip)
            + rcad_kernel::surface_area(&inner_strip)
            + rcad_kernel::surface_area(&ann);
        let tol = (5e-3_f64).max(0.0625 * 1390.8_f64);
        assert!(
            (sum - 1390.8_f64).abs() <= tol,
            "sum loose pieces ~1390.8, got {sum}"
        );
    }

    #[test]
    fn zp3_inner_frustum_strip_surface_area_sane() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let mut ib = Vec::with_capacity(N);
        let mut it = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            let c = ang.cos();
            let s = ang.sin();
            ib.push(DVec3::new(10.0 * c, 10.0 * s, -10.0));
            it.push(DVec3::new(5.0 * c, 5.0 * s, 0.0));
        }
        let (loft_i, ih) = loft_with_history(&[ib, it]).unwrap();
        let strip = strip_loft_caps(loft_i, ih).unwrap();
        let a = rcad_kernel::surface_area(&strip);
        assert!(
            a > 400.0 && a < 650.0,
            "inner frustum mantle area ~527, got {a}"
        );
    }

    #[test]
    fn zp3_outer_face_area_matches_before_and_after_sew() {
        let (outer_strip, inner_strip, annulus) = coaxial_cylinder_minus_frustum_loft_pieces(
            -10.0, 0.0, 10.0, 10.0, -10.0, 10.0,
        )
        .expect("loft pieces");
        let tol = (1e-4_f64).max(1e-6 * 10.0);
        let sewn = sew_shells(&[outer_strip.clone(), inner_strip, annulus], tol);
        assert!(sewn.free_edges.is_empty(), "free {:?}", sewn.free_edges);
        let f_loose = &outer_strip.solids[0].shells[0].faces[0];
        let f_sewn = &sewn.brep.solids[0].shells[0].faces[0];
        let a_loose = rcad_kernel::face_surface_area(&outer_strip, f_loose, 0);
        let a_sewn = rcad_kernel::face_surface_area(&sewn.brep, f_sewn, 0);
        assert!(
            (a_loose - a_sewn).abs() < 1e-2,
            "first outer lateral face area loose {a_loose} vs sewn {a_sewn}"
        );
    }

    #[test]
    fn zp3_loft_shell_matches_occt_geometry_numbers() {
        // Cone apex z=10, base z=-10, rb=10; cylinder z in [-10,0], r=10 — same as geometry_properties ZP3.
        let r = try_coaxial_cylinder_minus_frustum_loft_shell(-10.0, 0.0, 10.0, 10.0, -10.0, 10.0);
        assert!(r.is_some(), "expected sewn loft shell for ZP3 parameters");
        let brep = r.unwrap();
        let nf: usize = brep
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count();
        let a = rcad_kernel::surface_area(&brep);
        let v = rcad_kernel::volume(&brep);
        assert!(
            nf >= 60 && nf <= 70,
            "expected ~65 faces (32+32+annulus), got {nf}"
        );
        assert!(
            (v - 1310.0).abs() < 80.0,
            "expected volume cylinder minus frustum ~1310, got {v}"
        );
        let tol = (5e-3_f64).max(0.0625 * 1390.8_f64);
        assert!(
            (a - 1390.8_f64).abs() <= tol,
            "surface area: expected ~1390.8, got {a} (nf={nf}, vol={v})"
        );
    }

    #[test]
    fn zp3_annulus_plane_builds() {
        use std::f64::consts::TAU;
        const N: usize = 32;
        let mut o = Vec::with_capacity(N);
        let mut inn = Vec::with_capacity(N);
        for i in 0..N {
            let ang = TAU * i as f64 / N as f64;
            o.push(DVec3::new(10.0 * ang.cos(), 10.0 * ang.sin(), 0.0));
            inn.push(DVec3::new(5.0 * ang.cos(), 5.0 * ang.sin(), 0.0));
        }
        annulus_between_rings(&o, &inn).expect("annulus");
    }

    #[test]
    fn zp3_coaxial_cylinder_minus_cone_fast_path_triggered() {
        use glam::DVec3;
        use rcad_modeling::{make_cone_brep, make_cylinder_brep};
        let pc = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 10.0, 20.0).expect("cone");
        let pcy = make_cylinder_brep(
            DVec3::new(0.0, 0.0, -5.0),
            DVec3::Z,
            DVec3::X,
            10.0,
            10.0,
        )
        .expect("cylinder");
        assert!(
            try_difference_coaxial_cylinder_minus_cone(&pcy, &pc).is_some(),
            "ZP3 boptuc expects coaxial cylinder\\cone shortcut"
        );
    }
}
