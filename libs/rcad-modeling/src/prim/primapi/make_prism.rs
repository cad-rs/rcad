// OCCT BRepPrimAPI_MakePrism(face, gp_Vec(0,0,h)) 1:1 translation.
//
// The OCCT reference prism (pavefiller_stage_dump.cpp "PRISM") is built as a
// rectangular face in the XOY plane extruded by (0,0,h).  The 4 lateral faces
// of the sweep carry Geom_SurfaceOfLinearExtrusion surfaces (BRepSweep_Prism),
// so they are NOT planes: the FF classifies them as parametric and routes them
// through IntPatch_ImpPrm (walking).  This builder reproduces that structure:
// 2 planar caps + 4 linear-extrusion lateral faces.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, LinearExtrusionSurface, Plane, Surface3};
use rcad_kernel::topods::{self, Orientation, Shape};
use rcad_kernel::BRep;

pub fn make_prism_brep(
    width: f64,
    depth: f64,
    height: f64,
) -> Result<BRep, crate::BuildError> {
    let w = width.abs();
    let d = depth.abs();
    let h = height.abs();
    let mut t = BRep::new();
    let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };

    // Vertices: bottom 0..3 (z=0), top 4..7 (z=h).
    let v = [
        t.add_tvertex(DVec3::new(0.0, 0.0, 0.0)),
        t.add_tvertex(DVec3::new(w, 0.0, 0.0)),
        t.add_tvertex(DVec3::new(w, d, 0.0)),
        t.add_tvertex(DVec3::new(0.0, d, 0.0)),
        t.add_tvertex(DVec3::new(0.0, 0.0, h)),
        t.add_tvertex(DVec3::new(w, 0.0, h)),
        t.add_tvertex(DVec3::new(w, d, h)),
        t.add_tvertex(DVec3::new(0.0, d, h)),
    ];
    let ln = |a: DVec3, b: DVec3| Curve3::Line(Line3::new(a, b - a));
    // Bottom edges (z=0).
    let b_ed = [
        t.add_tedge(Some(ln(DVec3::new(0.0, 0.0, 0.0), DVec3::new(w, 0.0, 0.0))), v[0].clone(), v[1].clone(), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(w, 0.0, 0.0), DVec3::new(w, d, 0.0))), v[1].clone(), v[2].clone(), [0.0, d]),
        t.add_tedge(Some(ln(DVec3::new(w, d, 0.0), DVec3::new(0.0, d, 0.0))), v[2].clone(), v[3].clone(), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(0.0, d, 0.0), DVec3::new(0.0, 0.0, 0.0))), v[3].clone(), v[0].clone(), [0.0, d]),
    ];
    // Vertical edges.
    let e_ver = [
        t.add_tedge(Some(ln(DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 0.0, h))), v[0].clone(), v[4].clone(), [0.0, h]),
        t.add_tedge(Some(ln(DVec3::new(w, 0.0, 0.0), DVec3::new(w, 0.0, h))), v[1].clone(), v[5].clone(), [0.0, h]),
        t.add_tedge(Some(ln(DVec3::new(w, d, 0.0), DVec3::new(w, d, h))), v[2].clone(), v[6].clone(), [0.0, h]),
        t.add_tedge(Some(ln(DVec3::new(0.0, d, 0.0), DVec3::new(0.0, d, h))), v[3].clone(), v[7].clone(), [0.0, h]),
    ];
    // Top edges (z=h).
    let t_ed = [
        t.add_tedge(Some(ln(DVec3::new(0.0, 0.0, h), DVec3::new(w, 0.0, h))), v[4].clone(), v[5].clone(), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(w, 0.0, h), DVec3::new(w, d, h))), v[5].clone(), v[6].clone(), [0.0, d]),
        t.add_tedge(Some(ln(DVec3::new(w, d, h), DVec3::new(0.0, d, h))), v[6].clone(), v[7].clone(), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(0.0, d, h), DVec3::new(0.0, 0.0, h))), v[7].clone(), v[4].clone(), [0.0, d]),
    ];

    // Wires.  Bottom/top follow BRepPrim_GWedge; the lateral wires follow the
    // box (make_box.rs) lateral face topology.
    let wires = [
        t.add_twire(vec![b_ed[0].clone(), b_ed[1].clone(), b_ed[2].clone(), b_ed[3].clone()]),
        t.add_twire(vec![t_ed[0].clone(), rev(t_ed[3].clone()), rev(t_ed[2].clone()), rev(t_ed[1].clone())]),
        t.add_twire(vec![b_ed[0].clone(), e_ver[1].clone(), rev(t_ed[0].clone()), rev(e_ver[0].clone())]),
        t.add_twire(vec![rev(b_ed[2].clone()), e_ver[2].clone(), t_ed[2].clone(), rev(e_ver[3].clone())]),
        t.add_twire(vec![rev(b_ed[3].clone()), e_ver[0].clone(), t_ed[3].clone(), rev(e_ver[3].clone())]),
        t.add_twire(vec![b_ed[1].clone(), e_ver[2].clone(), rev(t_ed[1].clone()), rev(e_ver[1].clone())]),
    ];

    // Surfaces.  The two caps are planes; each lateral face is the linear
    // extrusion of its base edge along +Z (OCCT Geom_SurfaceOfLinearExtrusion).
    // Profile directions are chosen so that tangent x direction gives the
    // outward-pointing surface normal, matching the face orientation below.
    let pln = |pt: DVec3, n: DVec3, u: DVec3| Surface3::Plane(Plane {
        origin: pt,
        normal: n,
        u_dir: u,
        v_dir: n.cross(u).normalize_or_zero(),
    });
    let lext = |a: DVec3, b: DVec3| Surface3::LinearExtrusion(LinearExtrusionSurface {
        profile: Box::new(Curve3::Line(Line3::new(a, b - a))),
        direction: DVec3::Z,
    });
    // y=0 face: profile -X -> normal +Y; reversed -> outward -Y.
    // x=w face: profile +Y -> normal +X; forward  -> outward +X.
    // y=d face: profile -X -> normal +Y; forward  -> outward +Y.
    // x=0 face: profile +Y -> normal +X; reversed -> outward -X.
    let faces = [
        rev(t.add_tface(Some(pln(DVec3::ZERO, DVec3::Z, DVec3::X)), wires[0].clone(), vec![], None, None, vec![], true)),
        t.add_tface(Some(pln(DVec3::new(0.0, 0.0, h), DVec3::Z, DVec3::X)), wires[1].clone(), vec![], None, None, vec![], true),
        rev(t.add_tface(Some(lext(DVec3::new(w, 0.0, 0.0), DVec3::new(0.0, 0.0, 0.0))), wires[2].clone(), vec![], None, None, vec![], true)),
        t.add_tface(Some(lext(DVec3::new(w, d, 0.0), DVec3::new(0.0, d, 0.0))), wires[3].clone(), vec![], None, None, vec![], true),
        rev(t.add_tface(Some(lext(DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, d, 0.0))), wires[4].clone(), vec![], None, None, vec![], true)),
        t.add_tface(Some(lext(DVec3::new(w, 0.0, 0.0), DVec3::new(w, d, 0.0))), wires[5].clone(), vec![], None, None, vec![], true),
    ];
    let shell = t.add_tshell(faces.to_vec());
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn prism_brep(width: f64, depth: f64, height: f64) -> Result<BRep, crate::BuildError> {
    make_prism_brep(width, depth, height)
}
