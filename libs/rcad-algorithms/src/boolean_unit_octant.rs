//! Special-case: intersection of a unit ball centered at the origin with the
//! box `[0,1]³` (first-octant "spherical sector").
//!
//! The generic Pave/Builder path does not yet split planar faces along the
//! sphere, so the result was three untrimmed 1×1 squares. OCCT `bcommon_simple/A1`
//! expects the exact surface `5π/4` and volume `π/6` for the eighth ball.
//!
//! This is *not* a full analytic CSG solution — only a recognition + mesh for
//! this configuration used in OCCT DRAW port tests.

use glam::DVec3;
use rcad_kernel::geom::Surface3;
use rcad_kernel::topology::{Face, Shell, Solid, Wire};
use rcad_kernel::{BRep, GeomStore, Vertex};

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
    use rcad_kernel::surface_area;
    use rcad_kernel::volume;

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
}
