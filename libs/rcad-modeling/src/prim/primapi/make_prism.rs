// OCCT BRepPrimAPI_MakePrism(face, gp_Vec(0,0,h)) 1:1 translation.
//
// The OCCT reference prism (pavefiller_stage_dump.cpp "PRISM") is built as a
// rectangular face in the XOY plane extruded by (0,0,h).  The 4 lateral faces
// of the sweep are Geom_Plane faces (BRepAdaptor_Surface reports Geom_Plane /
// DynamicType Geom_Plane — the sweep of a line edge is a plane).  The FF
// classifies them as planes and intersects them analytically (Plane x quadric).
// This builder reproduces that structure: 2 planar caps + 4 planar lateral faces.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Plane, Surface3};
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
    // box (make_box.rs) lateral face topology.  The even-index wires (0,2,4)
    // belong to REVERSED faces — OCCT BRepPrim_GWedge::Face() reverses them
    // (i%2==0, ReverseFace) and TopoDS_Builder::Add then reverses the added
    // wire into the face, so the wires are stored REVERSED (see make_box.rs).
    let wires = [
        rev(t.add_twire(vec![b_ed[0].clone(), b_ed[1].clone(), b_ed[2].clone(), b_ed[3].clone()])),
        t.add_twire(vec![t_ed[0].clone(), rev(t_ed[3].clone()), rev(t_ed[2].clone()), rev(t_ed[1].clone())]),
        rev(t.add_twire(vec![b_ed[0].clone(), e_ver[1].clone(), rev(t_ed[0].clone()), rev(e_ver[0].clone())])),
        t.add_twire(vec![rev(b_ed[2].clone()), e_ver[2].clone(), t_ed[2].clone(), rev(e_ver[3].clone())]),
        rev(t.add_twire(vec![rev(b_ed[3].clone()), e_ver[0].clone(), t_ed[3].clone(), rev(e_ver[3].clone())])),
        t.add_twire(vec![b_ed[1].clone(), e_ver[2].clone(), rev(t_ed[1].clone()), rev(e_ver[1].clone())]),
    ];

    // Surfaces.  The two caps are planes; each lateral face is the planar image
    // of its base edge swept along +Z.  OCCT BRepPrimAPI_MakePrism represents
    // these as Geom_Plane faces (the sweep of a line edge is a plane), so rcad
    // builds them as Plane surfaces too — the FF pipeline then intersects them
    // analytically (Plane x quadric) instead of walking an extrusion surface.
    let pln = |pt: DVec3, n: DVec3, u: DVec3| Surface3::Plane(Plane {
        origin: pt,
        normal: n,
        u_dir: u,
        v_dir: n.cross(u).normalize_or_zero(),
    });
    // y=0 face: profile -X -> normal +Y; reversed -> outward -Y.
    // x=w face: profile +Y -> normal +X; forward  -> outward +X.
    // y=d face: profile -X -> normal +Y; forward  -> outward +Y.
    // x=0 face: profile +Y -> normal +X; reversed -> outward -X.
    // UV domains: the caps are [0,w]x[0,d]; each lateral face is the profile
    // parameter (the edge length) x the extrusion height.  OCCT prism faces
    // carry finite UV domains (BRepPrim_Build::MakeFace).
    let cap_uv = [0.0, w, 0.0, d];
    // OCCT BRepPrimAPI_MakePrism lateral faces carry the v-parameter running
    // from -1 (top) to 0 (bottom): the base-z=0 edge lies at v=0 (the uv bounds
    // are (0,1,-1,0)).  To reproduce that, the lateral v-dir is reversed here
    // (v=-z), so the z=0 boundary is at vmax=0.  This matches the reference
    // prism's FF domain classification (ClassifyLin2d rejects lines coincident
    // with the v=vmax boundary).
    let faces = [
        rev(t.add_tface(Some(pln(DVec3::ZERO, DVec3::Z, DVec3::X)), wires[0].clone(), vec![], None, Some(cap_uv), vec![], true)),
        t.add_tface(Some(pln(DVec3::new(0.0, 0.0, h), DVec3::Z, DVec3::X)), wires[1].clone(), vec![], None, Some(cap_uv), vec![], true),
        rev(t.add_tface(Some(pln(DVec3::new(w, 0.0, 0.0), DVec3::Y, DVec3::X)), wires[2].clone(), vec![], None, Some([0.0, w, -h, 0.0]), vec![], true)),
        t.add_tface(Some(pln(DVec3::new(w, d, 0.0), DVec3::Y, DVec3::X)), wires[3].clone(), vec![], None, Some([0.0, w, -h, 0.0]), vec![], true),
        rev(t.add_tface(Some(pln(DVec3::ZERO, DVec3::X, -DVec3::Y)), wires[4].clone(), vec![], None, Some([0.0, d, -h, 0.0]), vec![], true)),
        t.add_tface(Some(pln(DVec3::new(w, 0.0, 0.0), DVec3::X, -DVec3::Y)), wires[5].clone(), vec![], None, Some([0.0, d, -h, 0.0]), vec![], true),
    ];
    let shell = t.add_tshell(faces.to_vec());
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn prism_brep(width: f64, depth: f64, height: f64) -> Result<BRep, crate::BuildError> {
    make_prism_brep(width, depth, height)
}
