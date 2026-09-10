//! OCCT BRepBuilderAPI_MakeFace (TKTopAlgo) — face builders.
//!
//! OCCT BRepBuilderAPI_MakeFace(gp_Pln) builds a natural-restriction plane
//! face (no wires); the (Geom_Plane | Geom_CylindricalSurface, UMin, UMax,
//! VMin, VMax, TolDegen) forms build a bounded face whose outer wire is the
//! parameter-domain rectangle; MakeFace(wire) builds a face on the plane
//! fitted from the wire (BRepLib_MakeFace + BRepLib_FindSurface).
//!
//! Each bounded face attaches planar 2D curves (pcurves) to its edges
//! (BRepLib::SameParameter, BRepLib_MakeFace.cxx L250) so the
//! BRepGProp_Sinert boundary integral (base/gprop/surface) computes the area.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Circle3, Curve2d, Curve3, CylindricalSurface, Line2d, Line3, Plane, Surface3,
};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{BRep, Orientation};

use crate::BuildError;

/// OCCT BRepBuilderAPI_MakeFace(gp_Pln) — natural-restriction plane face.
pub fn make_face_plane_brep(plane: &Plane) -> Result<BRep, BuildError> {
    let mut brep = BRep::new();
    let w = brep.add_twire(vec![]);
    brep.add_tface(
        Some(Surface3::Plane(*plane)),
        w,
        vec![],
        None,
        None,
        vec![],
        true,
    );
    Ok(brep)
}

/// OCCT BRepBuilderAPI_MakeFace(Geom_Plane, UMin, UMax, VMin, VMax, TolDegen):
/// a bounded rectangular plane face over the parameter domain [u1,u2]x[v1,v2].
/// Area = (u2-u1)*(v2-v1).
pub fn make_face_plane_bounds_brep(
    plane: &Plane,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<BRep, BuildError> {
    let p = |u: f64, v: f64| plane.origin + plane.u_dir * u + plane.v_dir * v;
    let corners = [p(u1, v1), p(u2, v1), p(u2, v2), p(u1, v2)];
    let uv = [
        DVec2::new(u1, v1),
        DVec2::new(u2, v1),
        DVec2::new(u2, v2),
        DVec2::new(u1, v2),
    ];
    let curves = [
        Curve3::Line(Line3::new(corners[0], (corners[1] - corners[0]).normalize())),
        Curve3::Line(Line3::new(corners[1], (corners[2] - corners[1]).normalize())),
        Curve3::Line(Line3::new(corners[2], (corners[3] - corners[2]).normalize())),
        Curve3::Line(Line3::new(corners[3], (corners[0] - corners[3]).normalize())),
    ];
    let ranges = [
        [0.0, (corners[1] - corners[0]).length()],
        [0.0, (corners[2] - corners[1]).length()],
        [0.0, (corners[3] - corners[2]).length()],
        [0.0, (corners[0] - corners[3]).length()],
    ];
    bounded_face(Surface3::Plane(*plane), &corners, &uv, &curves, &ranges)
}

/// OCCT BRepBuilderAPI_MakeFace(Geom_CylindricalSurface, UMin, UMax, VMin,
/// VMax, TolDegen): a bounded cylindrical face over [u1,u2]x[v1,v2].
/// Area = radius * (u2-u1) * (v2-v1).
pub fn make_face_cylinder_bounds_brep(
    cyl: &CylindricalSurface,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Result<BRep, BuildError> {
    let x = cyl.ref_dir;
    let y = cyl.y_axis();
    let p = |u: f64, v: f64| {
        cyl.origin + (x * cyl.radius * u.cos() + y * cyl.radius * u.sin()) + cyl.axis * v
    };
    let corners = [p(u1, v1), p(u2, v1), p(u2, v2), p(u1, v2)];
    let uv = [
        DVec2::new(u1, v1),
        DVec2::new(u2, v1),
        DVec2::new(u2, v2),
        DVec2::new(u1, v2),
    ];
    // E1/E3 run along v=const (parameter u, circle); E2/E4 along u=const
    // (parameter v, straight line).  E3/E4 run in the reversed parameter
    // direction so the wire closes (E3: u2->u1, E4: v2->v1).
    let curves = [
        Curve3::Circle(Circle3 {
            center: cyl.origin + cyl.axis * v1,
            normal: cyl.axis,
            x_dir: x,
            y_dir: y,
            radius: cyl.radius,
        }),
        Curve3::Line(Line3::new(p(u2, v1), cyl.axis)),
        Curve3::Circle(Circle3 {
            center: cyl.origin + cyl.axis * v2,
            normal: cyl.axis,
            x_dir: x,
            y_dir: y,
            radius: cyl.radius,
        }),
        Curve3::Line(Line3::new(p(u1, v2), -cyl.axis)),
    ];
    let ranges = [[u1, u2], [0.0, v2 - v1], [u2, u1], [0.0, v2 - v1]];
    bounded_face(Surface3::Cylinder(*cyl), &corners, &uv, &curves, &ranges)
}

/// OCCT BRepBuilderAPI_MakeFace(wire) — face on the plane fitted from the
/// wire (BRepLib_MakeFace + BRepLib_FindSurface).  The wire must be planar;
/// the plane is fitted from the first non-collinear vertex triple and each
/// edge receives its planar 2D curve (BRepLib::SameParameter).
pub fn make_face_from_wire_brep(brep: &mut BRep, wire: Shape) -> Result<Shape, BuildError> {
    // Collect the distinct vertex positions of the wire edges.
    let wd = brep.wire(wire.clone()).clone();
    let mut pts: Vec<DVec3> = Vec::new();
    for e in &wd.edges {
        let ed = brep.edge(e.clone());
        for sv in [&ed.first, &ed.last] {
            let p = brep.vertex(sv.clone()).point;
            if !pts.iter().any(|q| (q - p).length_squared() < 1e-18) {
                pts.push(p);
            }
        }
    }
    if pts.len() < 3 {
        return Err(BuildError::DegenerateGeometry(
            "wire needs at least 3 points",
        ));
    }
    // Fit a plane from the first non-collinear triple.
    let p0 = pts[0];
    let mut plane = None;
    'outer: for i in 1..pts.len() {
        for j in (i + 1)..pts.len() {
            let n = (pts[i] - p0).cross(pts[j] - p0);
            if n.length_squared() > 1e-24 {
                plane = Some(Plane::new(p0, n));
                break 'outer;
            }
        }
    }
    let plane = plane.ok_or(BuildError::DegenerateGeometry("collinear wire"))?;
    let face = brep.add_tface(
        Some(Surface3::Plane(plane)),
        wire.clone(),
        vec![],
        None,
        None,
        vec![],
        false,
    );
    // BRepLib::SameParameter: attach the planar 2D curves to the edges.
    let face_key = (face.ptr_id(), face.location);
    for e in &wd.edges {
        let ed = brep.edge(e.clone());
        let a = brep.vertex(ed.first.clone()).point;
        let b = brep.vertex(ed.last.clone()).point;
        let to_uv = |p: DVec3| {
            DVec2::new(
                (p - plane.origin).dot(plane.u_dir),
                (p - plane.origin).dot(plane.v_dir),
            )
        };
        let uv_a = to_uv(a);
        let uv_b = to_uv(b);
        let d = uv_b - uv_a;
        let pc = Curve2d::Line(Line2d::new(
            uv_a,
            if d.length_squared() > 1e-30 {
                d.normalize()
            } else {
                DVec2::X
            },
        ));
        brep.edge_mut_inplace(e.clone())
            .pcurves
            .insert(face_key, (pc, 0.0, (b - a).length()));
    }
    Ok(face)
}

/// Build the bounded face from corner points, 2D pcurve endpoints and 3D edge
/// curves (each FORWARD over its range).  Attaches the planar 2D curves
/// (BRepLib::SameParameter) so the Gauss boundary integral can run.
fn bounded_face(
    surface: Surface3,
    corners: &[DVec3; 4],
    uv: &[DVec2; 4],
    curves: &[Curve3; 4],
    ranges: &[[f64; 2]; 4],
) -> Result<BRep, BuildError> {
    // OCCT BRepLib_MakeEdge::Init vertex orientation contract: the first
    // endpoint is stored FORWARD, the second REVERSED (TopExp::Vertices
    // identifies Vfirst/Vlast by orientation; see make_planar_rect_brep).
    let rev = |sr: Shape| Shape {
        orientation: Orientation::Reversed,
        ..sr
    };
    let mut brep = BRep::new();
    let vs: Vec<Shape> = corners.iter().map(|&c| brep.add_tvertex(c)).collect();
    let mut edges = Vec::with_capacity(4);
    for i in 0..4 {
        let j = (i + 1) % 4;
        let e = brep.add_tedge(
            Some(curves[i].clone()),
            vs[i].clone(),
            rev(vs[j].clone()),
            ranges[i],
        );
        edges.push(e);
    }
    let wire = brep.add_twire(edges.clone());
    let face = brep.add_tface(Some(surface), wire, vec![], None, None, vec![], false);
    let face_key = (face.ptr_id(), face.location);
    for (i, e) in edges.into_iter().enumerate() {
        let a = uv[i];
        let b = uv[(i + 1) % 4];
        let d = b - a;
        let len = d.length();
        let pc = Curve2d::Line(Line2d::new(a, if len > 1e-30 { d / len } else { DVec2::X }));
        brep.edge_mut_inplace(e).pcurves.insert(face_key, (pc, 0.0, len));
    }
    Ok(brep)
}
