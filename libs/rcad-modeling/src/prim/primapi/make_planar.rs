//! Planar face construction helpers (legacy API used by generated DRAW tests).
//!
//! OCCT `mkplane` equivalent: build a single rectangular planar face.

use glam::{DVec2, DVec3};
use rcad_kernel::BRep;
use rcad_kernel::geom::{Curve2d, Curve3, Line2d, Line3, Plane, Surface3};
use rcad_kernel::topods::Shape;

/// Build a rectangular planar face with the given UV extents.
///
/// Mirrors the removed legacy `make_planar_rect_brep` (originally built on the
/// old builder API) using the current `topods::BRep` pool. The result is a
/// single plane face (OCCT `mkplane` semantics).
pub fn make_planar_rect_brep(
    origin: DVec3,
    u_axis: DVec3,
    v_axis: DVec3,
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
) -> Result<BRep, crate::BuildError> {
    let normal = u_axis.cross(v_axis).normalize();
    let surface = Surface3::Plane(Plane::new(origin, normal));

    let c0 = origin + u_axis * umin + v_axis * vmin;
    let c1 = origin + u_axis * umax + v_axis * vmin;
    let c2 = origin + u_axis * umax + v_axis * vmax;
    let c3 = origin + u_axis * umin + v_axis * vmax;

    let mut brep = BRep::new();
    let v0 = brep.add_tvertex(c0);
    let v1 = brep.add_tvertex(c1);
    let v2 = brep.add_tvertex(c2);
    let v3 = brep.add_tvertex(c3);

    let line = |p0: DVec3, p1: DVec3| -> Curve3 {
        let d = p1 - p0;
        Curve3::Line(Line3 {
            origin: p0,
            direction: d.normalize(),
        })
    };
    let e0 = brep.add_tedge(Some(line(c0, c1)), v0.clone(), v1.clone(), [0.0, (c1 - c0).length()]);
    let e1 = brep.add_tedge(Some(line(c1, c2)), v1.clone(), v2.clone(), [0.0, (c2 - c1).length()]);
    let e2 = brep.add_tedge(Some(line(c2, c3)), v2.clone(), v3.clone(), [0.0, (c3 - c2).length()]);
    let e3 = brep.add_tedge(Some(line(c3, c0)), v3.clone(), v0.clone(), [0.0, (c0 - c3).length()]);

    let wire = brep.add_twire(vec![e0, e1, e2, e3]);
    brep.add_tface(Some(surface), wire, vec![], None, None, vec![], false);
    Ok(brep)
}

/// OCCT `BRepLib_MakeFace(W, OnlyPlane=true)` translation
/// (BRepLib_MakeFace.cxx L180-258, Add L873-899, CheckInside L903-925).
///
/// The surface is the `gp_Pln(gp_Pnt, gp_Dir)` plane built from the caller
/// origin and normal (BRepLib_FindSurface.cxx L546-575).  The outer polygon
/// and each hole become a closed wire added to the face; every edge receives
/// its planar 2D curve (BRepLib::SameParameter, BRepLib_MakeFace.cxx L250).
pub fn make_planar_polygon_brep(
    origin: DVec3,
    normal: DVec3,
    outer: &[DVec3],
    holes: &[Vec<DVec3>],
) -> Result<BRep, crate::BuildError> {
    if outer.len() < 3 {
        return Err(crate::BuildError::DegenerateGeometry(
            "polygon needs at least 3 points",
        ));
    }
    let plane = Plane::new(origin, normal);
    let mut brep = BRep::new();
    let (outer_wire, outer_edges) = make_polygon_wire(&mut brep, outer)?;
    let mut inner_wires = Vec::new();
    let mut inner_edges = Vec::new();
    let mut hole_polys: Vec<Vec<DVec3>> = Vec::new();
    for poly in holes {
        // A valid planar face orients the hole wires opposite to the outer
        // wire (the hole interior stays OUT of the bounded face — OCCT
        // BRepGProp_Gauss boundary integral subtracts a CW hole).  Reverse the
        // point chain so the hole boundary runs clockwise.
        let reversed: Vec<DVec3> = poly.iter().rev().copied().collect();
        let (w, edges) = make_polygon_wire(&mut brep, &reversed)?;
        inner_wires.push(w);
        inner_edges.push(edges);
        hole_polys.push(reversed);
    }
    let face = brep.add_tface(
        Some(Surface3::Plane(plane)),
        outer_wire,
        inner_wires,
        None,
        None,
        vec![],
        false,
    );
    let face_key = (face.ptr_id(), face.location);
    let project = |p: DVec3| -> DVec2 {
        DVec2::new((p - origin).dot(plane.u_dir), (p - origin).dot(plane.v_dir))
    };
    // OCCT BRepLib::SameParameter: attach the planar 2D curves to the edges.
    let polygons = std::iter::once(outer).chain(hole_polys.iter().map(|v| v.as_slice()));
    let mut edge_groups = vec![outer_edges];
    edge_groups.extend(inner_edges);
    for (poly, edges) in polygons.zip(edge_groups) {
        for (k, e) in edges.into_iter().enumerate() {
            let p0 = project(poly[k]);
            let p1 = project(poly[(k + 1) % poly.len()]);
            let d = p1 - p0;
            let dir = if d.length_squared() > 0.0 {
                d.normalize()
            } else {
                DVec2::X
            };
            let len = (poly[(k + 1) % poly.len()] - poly[k]).length();
            brep.edge_mut_inplace(e).pcurves.insert(
                face_key,
                (Curve2d::Line(Line2d::new(p0, dir)), 0.0, len),
            );
        }
    }
    // OCCT CheckInside (BRepLib_MakeFace.cxx L903-925): a bounded planar wire
    // classifies the infinite 2D point OUT, so no wire reversal is produced.
    Ok(brep)
}

fn make_polygon_wire(
    brep: &mut BRep,
    points: &[DVec3],
) -> Result<(Shape, Vec<Shape>), crate::BuildError> {
    let n = points.len();
    if n < 3 {
        return Err(crate::BuildError::DegenerateGeometry(
            "polygon needs at least 3 points",
        ));
    }
    let mut vertices = Vec::with_capacity(n);
    for p in points {
        if !p.is_finite() {
            return Err(crate::BuildError::NonFiniteValue("polygon vertex"));
        }
        vertices.push(brep.add_tvertex(*p));
    }
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        let j = (i + 1) % n;
        let delta = points[j] - points[i];
        let len = delta.length();
        if len < 1e-12 {
            return Err(crate::BuildError::DegenerateGeometry(
                "zero-length polygon edge",
            ));
        }
        let e = brep.add_tedge(
            Some(Curve3::Line(Line3::new(points[i], delta / len))),
            vertices[i].clone(),
            vertices[j].clone(),
            [0.0, len],
        );
        edges.push(e);
    }
    let wire = brep.add_twire(edges.clone());
    Ok((wire, edges))
}
