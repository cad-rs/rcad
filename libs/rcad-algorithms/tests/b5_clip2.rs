use glam::DVec3;
use rcad_algorithms::inttools::edge_face::clip_line_to_polygon_with_tol;
use rcad_kernel::geom::{Line3, Plane};

#[test]
fn test_clip_exact() {
    // Use the EXACT planes and line from the PaveFiller debug:
    // p1=(0.000,1.000,1.000) n1=(1.000,0.000,0.000)  [A-face[4] promoted BSpline]
    // p2=(0.000,0.250,0.000) n2=(-0.000,-1.000,-0.000)  [B-face[8] Plane]
    // line o=(-0.000,0.250,0.000) d=(0.000,0.000,-1.000)

    let plane_a = Plane { origin: DVec3::new(0.0, 1.0, 1.0), normal: DVec3::new(1.0, 0.0, 0.0) };
    let plane_b = Plane { origin: DVec3::new(0.0, 0.25, 0.0), normal: DVec3::new(0.0, -1.0, 0.0) };
    let line = Line3 { origin: DVec3::new(0.0, 0.25, 0.0), direction: DVec3::new(0.0, 0.0, -1.0) };

    // For A-face[4] boundary (x=0 face)
    let verts_a = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 1.0),
        DVec3::new(0.0, 0.0, 1.0),
    ];

    // For B-face[8] boundary (y=0.25 face)
    let verts_b = vec![
        DVec3::new(0.0, 0.25, 0.0),
        DVec3::new(1.0, 0.25, 0.0),
        DVec3::new(1.0, 0.25, 1.0),
        DVec3::new(0.0, 0.25, 1.0),
    ];

    let tol = 1e-7;
    let ra = clip_line_to_polygon_with_tol(&line, &plane_a, &verts_a, tol);
    let rb = clip_line_to_polygon_with_tol(&line, &plane_b, &verts_b, tol);
    eprintln!("EXACT test: ranges_a={:?} ranges_b={:?}", ra, rb);
}
