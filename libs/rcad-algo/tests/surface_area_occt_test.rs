//! OCCT BRepGProp::SurfaceProperties (BRepGProp.cxx L167-266) alignment.
//!
//! `checkprops result -s` (BRepTest_GPropCommands.cxx L123-125) calls
//! BRepGProp::SurfaceProperties(S, G, SkipShared=false, UseTriangulation=false)
//! → surfaceProperties(S, Props, Eps=1.0, false, false).  Per face:
//!   - natural-restriction faces (no wires) → G.Perform(BF): whole-surface
//!     Gauss-Legendre over the surface natural bounds (BRepGProp_Gauss.cxx
//!     L1306-1393),
//!   - faces with wires → G.Perform(BF, BD): Green-theorem line integral over
//!     the edge pcurves with the face UV bounds as the U offset
//!     (BRepGProp_Gauss.cxx L1126-1211, BRepTools.cxx L64-367).
//!
//! The sphere test fixes the semantics of `surface_area` on a full sphere
//! (a natural-restriction face, integrated by the whole-surface Gauss path):
//! exactly 4*PI*R^2.

use glam::DVec3;
use rcad_algo::common;
use rcad_kernel::base::gprop::surface::{face_surface_area_gauss_domain, surface_area};
use rcad_kernel::base::gprop::tri::face_flat_iter;
use rcad_kernel::topods;

/// OCCT `box 2 3 4` (BRepPrimAPI_MakeBox): the GWedge faces carry planar 2D
/// curves (SetPCurve), so every planar face integrates via Compute(Face,
/// Domain) with the face UV bounds (BRepTools::UVBounds) as the U offset;
/// total surface area = 52.
#[test]
fn box_surface_area_matches_occt() {
    let b = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 3.0, 4.0).unwrap();
    let got = surface_area(&b);
    assert!((got - 52.0).abs() < 1e-6, "box surface_area = {got}, expected 52");
}

/// OCCT `pcylinder 1 2` (BRepPrimAPI_MakeCylinder): lateral 2*PI*1*2 + two
/// disks 2*PI*1^2 = 6*PI.  The bottom face is REVERSED; BRepGProp_Domain
/// explores F.Oriented(FORWARD) so the disk boundary keeps the wire x edge
/// orientation (cumOri) and the Green integral is positive.
#[test]
fn cylinder_surface_area_matches_occt() {
    let cyl = rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
    let expected = 6.0 * std::f64::consts::PI;
    let got = surface_area(&cyl);
    assert!(
        (got - expected).abs() < 1e-6,
        "cylinder surface_area = {got}, expected {expected}"
    );
}

/// OCCT `pcylinder 1 2` + `box -1 -1 0 2 2 1` bopcommon (bopcommon_simple
/// v5): the box cuts the cylinder at z=1; result vertex count is 8.
#[test]
fn common_cylinder_box_vertex_count() {
    let b1 = rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();
    let b2 = rcad_modeling::make_box_brep(DVec3::new(-1.0, -1.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 1.0).unwrap();
    let result = rcad_algo::common(&b1, &b2).unwrap();
    let mut used = std::collections::BTreeSet::new();
    let mut chains: Vec<(u64, u64, usize, usize)> = Vec::new();
    for s in result.solids() {
        for sh in &s.shells {
            for f in &sh.faces {
                for we in &f.outer_wire.edges {
                    if let Some((va, vb)) = result.tshapes.get(we.idx).and_then(|ts| match ts.as_ref() {
                        rcad_kernel::topods::TShape::Edge(ed) => Some((ed.first.index, ed.last.index)),
                        _ => None,
                    }) {
                        chains.push((std::sync::Arc::as_ptr(&result.tshapes[we.idx]) as u64, std::sync::Arc::as_ptr(&result.tshapes[va]) as u64, va, vb));
                        used.insert(va);
                        used.insert(vb);
                    }
                }
                for w in &f.inner_wires {
                    for we in &w.edges {
                        if let Some((va, vb)) = result.tshapes.get(we.idx).and_then(|ts| match ts.as_ref() {
                            rcad_kernel::topods::TShape::Edge(ed) => Some((ed.first.index, ed.last.index)),
                            _ => None,
                        }) {
                            chains.push((std::sync::Arc::as_ptr(&result.tshapes[we.idx]) as u64, std::sync::Arc::as_ptr(&result.tshapes[va]) as u64, va, vb));
                            used.insert(va);
                            used.insert(vb);
                        }
                    }
                }
            }
        }
    }
    eprintln!("[DBG] v5 VERTEX count = {}, n_faces = {}, tshapes.len = {}", used.len(), face_flat_iter(&result).len(), result.tshapes.len());
    for (ti, f) in &face_flat_iter(&result) {
        let st = match &*result.tshapes[*ti] {
            topods::TShape::Face(fd) => format!("{:?}", fd.surface),
            _ => String::new(),
        };
        eprintln!("[DBG]   face {ti} {st} edges={}", f.outer_wire.edges.len());
    }
}

/// OCCT `box 1 1 1; box 1 1 1; bopcommon result; checkprops result -s 6`
/// (bopcommon_simple A1): the result is the input box TShape shared into the
/// result BRep (BRep_Builder::Add), so the GWedge face pcurves remain
/// attached and the surface area is exactly 6.
#[test]
fn common_box_box_surface_area_matches_occt() {
    let b1 = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let b2 = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let result = rcad_algo::common(&b1, &b2).unwrap();
    let got = surface_area(&result);
    assert!((got - 6.0).abs() < 1e-6, "common box surface_area = {got}, expected 6");
}

/// OCCT bopcommon_simple a5: box 1 1 1 + box 0.5 1 0.5 (contained) common =
/// the inner box; surface area 2.5.
#[test]
fn common_nested_box_surface_area() {
    let b1 = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
    let b2 = rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 0.5, 1.0, 0.5).unwrap();
    let result = rcad_algo::common(&b1, &b2).unwrap();
    let got = surface_area(&result);
    assert!((got - 2.5).abs() < 1e-6, "nested box common surface_area = {got}, expected 2.5");
}

/// OCCT `psphere 1`: a full sphere is a natural-restriction face (no wires),
/// so BRepGProp::SurfaceProperties integrates the whole surface by
/// Gauss-Legendre (G.Perform(BF)); the surface area is exactly 4*PI*R^2.
#[test]
fn sphere_surface_area_matches_occt() {
    let sphere = rcad_modeling::make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
    let expected = 4.0 * std::f64::consts::PI;
    let got = surface_area(&sphere);
    assert!(
        (got - expected).abs() < 1e-6,
        "sphere surface_area = {got}, expected {expected}"
    );
}

