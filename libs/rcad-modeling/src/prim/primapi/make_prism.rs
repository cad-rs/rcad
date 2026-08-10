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
    // OCCT BRepLib_MakeEdge forces V1=FWD, V2=REV (BRepLib_MakeEdge.cxx L772-774);
    // BRepBuilderAPI_MakePrism edges are built via MakeEdge, so the last vertex
    // of every prism edge carries the REVERSED orientation.
    let rev_v = |v: &Shape| Shape { orientation: rcad_kernel::topods::Orientation::Reversed, ..v.clone() };
    // Bottom edges (z=0).
    let b_ed = [
        t.add_tedge(Some(ln(DVec3::new(0.0, 0.0, 0.0), DVec3::new(w, 0.0, 0.0))), v[0].clone(), rev_v(&v[1]), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(w, 0.0, 0.0), DVec3::new(w, d, 0.0))), v[1].clone(), rev_v(&v[2]), [0.0, d]),
        t.add_tedge(Some(ln(DVec3::new(w, d, 0.0), DVec3::new(0.0, d, 0.0))), v[2].clone(), rev_v(&v[3]), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(0.0, d, 0.0), DVec3::new(0.0, 0.0, 0.0))), v[3].clone(), rev_v(&v[0]), [0.0, d]),
    ];
    // Vertical edges.
    let e_ver = [
        t.add_tedge(Some(ln(DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.0, 0.0, h))), v[0].clone(), rev_v(&v[4]), [0.0, h]),
        t.add_tedge(Some(ln(DVec3::new(w, 0.0, 0.0), DVec3::new(w, 0.0, h))), v[1].clone(), rev_v(&v[5]), [0.0, h]),
        t.add_tedge(Some(ln(DVec3::new(w, d, 0.0), DVec3::new(w, d, h))), v[2].clone(), rev_v(&v[6]), [0.0, h]),
        t.add_tedge(Some(ln(DVec3::new(0.0, d, 0.0), DVec3::new(0.0, d, h))), v[3].clone(), rev_v(&v[7]), [0.0, h]),
    ];
    // Top edges (z=h).
    let t_ed = [
        t.add_tedge(Some(ln(DVec3::new(0.0, 0.0, h), DVec3::new(w, 0.0, h))), v[4].clone(), rev_v(&v[5]), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(w, 0.0, h), DVec3::new(w, d, h))), v[5].clone(), rev_v(&v[6]), [0.0, d]),
        t.add_tedge(Some(ln(DVec3::new(w, d, h), DVec3::new(0.0, d, h))), v[6].clone(), rev_v(&v[7]), [0.0, w]),
        t.add_tedge(Some(ln(DVec3::new(0.0, d, h), DVec3::new(0.0, 0.0, h))), v[7].clone(), rev_v(&v[4]), [0.0, d]),
    ];

    // Wires.  Each lateral wire stores [vertical edge, next vertical edge,
    // base edge, top edge]; the caps store their 4 boundary edges.  Edge flags
    // reproduce OCCT BRepPrimAPI_MakePrism exactly (verified against the
    // OCCT DS face dumps, P029 prism 1x1x1 x cylinder):
    //   y=0 face (Rev): [e_ver0, rev(e_ver1), rev(b_ed0), t_ed0] -> R,F,F,R
    //   x=w face (Rev): [e_ver1, rev(e_ver2), rev(b_ed1), t_ed1] -> R,F,F,R
    //   y=d face (Rev): [e_ver2, rev(e_ver3), rev(b_ed2), t_ed2] -> R,F,F,R
    //   x=0 face (Rev): [e_ver3, rev(e_ver0), rev(b_ed3), t_ed3] -> R,F,F,R
    //   z=0 face (Rev): [b_ed0, b_ed1, b_ed2, b_ed3]             -> R,R,R,R
    //   z=h face (Fwd): [t_ed0, t_ed1, t_ed2, t_ed3]             -> F,F,F,F
    // The composed orientation of every shared edge is opposite on the two
    // adjacent faces, which is what the BuilderSolid ShellSplitter's
    // GetEdgeOff relies on.
    let wires = [
        t.add_twire(vec![e_ver[0].clone(), rev(e_ver[1].clone()), rev(b_ed[0].clone()), t_ed[0].clone()]),
        t.add_twire(vec![e_ver[1].clone(), rev(e_ver[2].clone()), rev(b_ed[1].clone()), t_ed[1].clone()]),
        t.add_twire(vec![e_ver[2].clone(), rev(e_ver[3].clone()), rev(b_ed[2].clone()), t_ed[2].clone()]),
        t.add_twire(vec![e_ver[3].clone(), rev(e_ver[0].clone()), rev(b_ed[3].clone()), t_ed[3].clone()]),
        t.add_twire(vec![b_ed[0].clone(), b_ed[1].clone(), b_ed[2].clone(), b_ed[3].clone()]),
        t.add_twire(vec![t_ed[0].clone(), t_ed[1].clone(), t_ed[2].clone(), t_ed[3].clone()]),
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
        // y=0 face (REVERSED, outward -Y). Surface per OCCT MakePrism:
        // N(+Y), u=X, v=-Z (BRepSweep_Prism profile sweep).
        rev(t.add_tface(Some(pln(DVec3::new(w, 0.0, 0.0), DVec3::Y, DVec3::X)), wires[0].clone(), vec![], None, Some([0.0, w, -h, 0.0]), vec![], true)),
        // x=w face (REVERSED, outward +X). Surface per OCCT MakePrism:
        // N(-X), u=Y, v=-Z — the sweep of the +Y profile edge along +Z gives
        // normal -X (X = profile dir × sweep dir), face REVERSED to flip out.
        rev(t.add_tface(Some(pln(DVec3::new(w, 0.0, 0.0), -DVec3::X, DVec3::Y)), wires[1].clone(), vec![], None, Some([0.0, d, -h, 0.0]), vec![], true)),
        // y=d face (REVERSED, outward +Y). Surface per OCCT MakePrism:
        // N(-Y), u=-X, v=-Z.
        rev(t.add_tface(Some(pln(DVec3::new(w, d, 0.0), -DVec3::Y, -DVec3::X)), wires[2].clone(), vec![], None, Some([0.0, w, -h, 0.0]), vec![], true)),
        // x=0 face (REVERSED, outward -X). Surface per OCCT MakePrism:
        // N(+X), u=-Y, v=-Z.
        rev(t.add_tface(Some(pln(DVec3::ZERO, DVec3::X, -DVec3::Y)), wires[3].clone(), vec![], None, Some([0.0, d, -h, 0.0]), vec![], true)),
        // z=0 face (REVERSED, outward -Z). Surface per OCCT MakePrism:
        // N(+Z), u=X, v=Y.
        rev(t.add_tface(Some(pln(DVec3::ZERO, DVec3::Z, DVec3::X)), wires[4].clone(), vec![], None, Some(cap_uv), vec![], true)),
        // z=h face (FORWARD, outward +Z). Surface per OCCT MakePrism:
        // N(+Z), u=X, v=Y.
        t.add_tface(Some(pln(DVec3::new(0.0, 0.0, h), DVec3::Z, DVec3::X)), wires[5].clone(), vec![], None, Some(cap_uv), vec![], true),
    ];
    let shell = t.add_tshell(faces.to_vec());
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn prism_brep(width: f64, depth: f64, height: f64) -> Result<BRep, crate::BuildError> {
    make_prism_brep(width, depth, height)
}
