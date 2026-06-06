use glam::DVec3;
use rcad_algorithms::inttools::edge_face::clip_line_to_polygon_with_tol;
use rcad_kernel::geom::{Line3, Plane};

#[test]
fn test_clip() {
    // B-face[8] at y=0.25
    let plane = Plane { origin: DVec3::new(0.0, 0.25, 0.0), normal: DVec3::new(0.0, -1.0, 0.0) };
    
    // Line: intersection of x=0 and y=0.25
    let line = Line3 { origin: DVec3::new(0.0, 0.25, 0.0), direction: DVec3::new(0.0, 0.0, 1.0) };
    
    // Face boundary vertices
    let verts = vec![
        DVec3::new(0.0, 0.25, 0.0),
        DVec3::new(1.0, 0.25, 0.0),
        DVec3::new(1.0, 0.25, 1.0),
        DVec3::new(0.0, 0.25, 1.0),
    ];
    
    let clip_tol = 1e-7;
    let ranges = clip_line_to_polygon_with_tol(&line, &plane, &verts, clip_tol);
    eprintln!("ranges for B-face[8]: {:?}", ranges);

    // Also test for A-face[4] at x=0
    let plane_a = Plane { origin: DVec3::new(0.0, -0.5, 1.5), normal: DVec3::new(-1.0, 0.0, 0.0) };
    let verts_a = vec![
        DVec3::new(0.0, 0.0, 0.0),
        DVec3::new(0.0, 1.0, 0.0),
        DVec3::new(0.0, 1.0, 1.0),
        DVec3::new(0.0, 0.0, 1.0),
    ];
    let ranges_a = clip_line_to_polygon_with_tol(&line, &plane_a, &verts_a, clip_tol);
    eprintln!("ranges for A-face[4]: {:?}", ranges_a);
}
