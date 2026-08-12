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
    // x=w face: profile +Y -> normal -X (profile dir × sweep dir); reversed
    //           -> outward +X.
    // y=d face: profile -X -> normal -Y; reversed -> outward +Y.
    // x=0 face: profile +Y -> normal +X; reversed -> outward -X.
    // OCCT BRepPrimAPI_MakePrism lateral faces are Geom_Plane with infinite
    // UV bounds (Geom_Plane::Bounds = RealLast); uv_domain is None to match
    // (make_box.rs uses None too).  The v-parameter runs with the sweep: the
    // base-z=0 edge lies at v=0 (v=-z, so the top edge is at v=-h).
    let faces = [
        // y=0 face (REVERSED, outward -Y). Surface per OCCT MakePrism:
        // N(+Y), u=X, v=-Z; origin (0,0,0); infinite UV bounds (Geom_Plane).
        rev(t.add_tface(Some(pln(DVec3::new(0.0, 0.0, 0.0), DVec3::Y, DVec3::X)), wires[0].clone(), vec![], None, None, vec![], true)),
        // x=w face (REVERSED, outward +X). Surface per OCCT MakePrism:
        // N(-X), u=Y, v=-Z — the sweep of the +Y profile edge along +Z gives
        // normal -X (X = profile dir × sweep dir), face REVERSED to flip out.
        rev(t.add_tface(Some(pln(DVec3::new(w, 0.0, 0.0), -DVec3::X, DVec3::Y)), wires[1].clone(), vec![], None, None, vec![], true)),
        // y=d face (REVERSED, outward +Y). Surface per OCCT MakePrism:
        // N(-Y), u=-X, v=-Z; origin (w,d,0).
        rev(t.add_tface(Some(pln(DVec3::new(w, d, 0.0), -DVec3::Y, -DVec3::X)), wires[2].clone(), vec![], None, None, vec![], true)),
        // x=0 face (REVERSED, outward -X). Surface per OCCT MakePrism:
        // N(+X), u=-Y, v=-Z; origin (0,d,0).
        rev(t.add_tface(Some(pln(DVec3::new(0.0, d, 0.0), DVec3::X, -DVec3::Y)), wires[3].clone(), vec![], None, None, vec![], true)),
        // z=0 face (REVERSED, outward -Z). Surface per OCCT MakePrism:
        // N(+Z), u=X, v=Y; origin (w/2,d/2,0) — the cap centre.
        rev(t.add_tface(Some(pln(DVec3::new(0.5 * w, 0.5 * d, 0.0), DVec3::Z, DVec3::X)), wires[4].clone(), vec![], None, None, vec![], true)),
        // z=h face (FORWARD, outward +Z). Surface per OCCT MakePrism:
        // N(+Z), u=X, v=Y; origin (w/2,d/2,h) — the cap centre.
        t.add_tface(Some(pln(DVec3::new(0.5 * w, 0.5 * d, h), DVec3::Z, DVec3::X)), wires[5].clone(), vec![], None, None, vec![], true),
    ];
    let shell = t.add_tshell(faces.to_vec());
    t.add_tsolid(vec![shell]);
    Ok(t)
}

pub fn prism_brep(width: f64, depth: f64, height: f64) -> Result<BRep, crate::BuildError> {
    make_prism_brep(width, depth, height)
}

/// OCCT BRepPrimAPI_MakePrism(face, vec) 1:1 translation for a rectangular
/// profile face. The profile face lies in the plane spanned by `x_dir`/`y_dir`
/// at `origin` with extents `width`×`height`; it is extruded by `extr`.
///
/// The extruded (top) cap reuses the source cap's TShape with a Location
/// translation — exactly like OCCT BRepSweep_Trsf::Process (TShape + Location,
/// no copy): the top vertices/edges are `Shape { location, ..source }` sharing
/// the source TShape. The lateral faces are new planar faces whose wires are
/// [sweep_i, rev(sweep_{i+1}), rev(base_i), top_i] with the same orientations
/// as OCCT BRepPrimAPI_MakePrism (verified against DRAW `prism` f1: all faces
/// FORWARD, lateral normal = extr_dir × edge_dir).
pub fn make_prism_from_face_brep(
    origin: DVec3,
    x_dir: DVec3,
    y_dir: DVec3,
    width: f64,
    height: f64,
    extr: DVec3,
) -> Result<BRep, crate::BuildError> {
    let xd = x_dir.normalize_or_zero();
    let yd = y_dir.normalize_or_zero();
    let ext_len = extr.length();
    if xd == DVec3::ZERO || yd == DVec3::ZERO || ext_len < 1e-12 {
        return Err(crate::BuildError::DegenerateGeometry("prism zero axis"));
    }
    let mut t = BRep::new();
    // Location of the extruded cap (OCCT BRepSweep_Prism::Location).
    let loc = t.add_location(glam::DAffine3::from_translation(extr));
    let rev = |sr: Shape| Shape { orientation: Orientation::Reversed, ..sr };
    let rev_v = |v: &Shape| Shape { orientation: Orientation::Reversed, ..v.clone() };
    let vpt = |v: &Shape| v.as_vertex().expect("prism profile vertex").point;

    // Profile face vertices (OCCT f1: (4,2,0),(4,3,0),(4,3,1),(4,2,1)).
    let v = [
        t.add_tvertex(origin),
        t.add_tvertex(origin + xd * width),
        t.add_tvertex(origin + xd * width + yd * height),
        t.add_tvertex(origin + yd * height),
    ];
    // Extruded vertices: same TShape + Location (OCCT BRepSweep_Trsf::Process).
    let ve = [
        Shape { location: loc, ..v[0].clone() },
        Shape { location: loc, ..v[1].clone() },
        Shape { location: loc, ..v[2].clone() },
        Shape { location: loc, ..v[3].clone() },
    ];
    let ln = |a: DVec3, b: DVec3| Curve3::Line(Line3::new(a, b - a));
    // Base (profile) edges: v_i -> v_{i+1}.
    let b_ed = [
        t.add_tedge(Some(ln(vpt(&v[0]), vpt(&v[1]))), v[0].clone(), rev_v(&v[1]), [0.0, width]),
        t.add_tedge(Some(ln(vpt(&v[1]), vpt(&v[2]))), v[1].clone(), rev_v(&v[2]), [0.0, height]),
        t.add_tedge(Some(ln(vpt(&v[2]), vpt(&v[3]))), v[2].clone(), rev_v(&v[3]), [0.0, width]),
        t.add_tedge(Some(ln(vpt(&v[3]), vpt(&v[0]))), v[3].clone(), rev_v(&v[0]), [0.0, height]),
    ];
    // Top edges: same TShape + Location.
    let t_ed = [
        Shape { location: loc, ..b_ed[0].clone() },
        Shape { location: loc, ..b_ed[1].clone() },
        Shape { location: loc, ..b_ed[2].clone() },
        Shape { location: loc, ..b_ed[3].clone() },
    ];
    // Sweep edges: v_i -> v_i + extr (new TShapes).
    let e_ver = [
        t.add_tedge(Some(ln(vpt(&v[0]), vpt(&v[0]) + extr)), v[0].clone(), rev_v(&ve[0]), [0.0, ext_len]),
        t.add_tedge(Some(ln(vpt(&v[1]), vpt(&v[1]) + extr)), v[1].clone(), rev_v(&ve[1]), [0.0, ext_len]),
        t.add_tedge(Some(ln(vpt(&v[2]), vpt(&v[2]) + extr)), v[2].clone(), rev_v(&ve[2]), [0.0, ext_len]),
        t.add_tedge(Some(ln(vpt(&v[3]), vpt(&v[3]) + extr)), v[3].clone(), rev_v(&ve[3]), [0.0, ext_len]),
    ];
    // Lateral wires (OCCT BRepPrimAPI_MakePrism order):
    //   [e_ver_i, rev(e_ver_{i+1}), rev(b_ed_i), t_ed_i] — all FORWARD faces.
    let wires = [
        t.add_twire(vec![e_ver[0].clone(), rev(e_ver[1].clone()), rev(b_ed[0].clone()), t_ed[0].clone()]),
        t.add_twire(vec![e_ver[1].clone(), rev(e_ver[2].clone()), rev(b_ed[1].clone()), t_ed[1].clone()]),
        t.add_twire(vec![e_ver[2].clone(), rev(e_ver[3].clone()), rev(b_ed[2].clone()), t_ed[2].clone()]),
        t.add_twire(vec![e_ver[3].clone(), rev(e_ver[0].clone()), rev(b_ed[3].clone()), t_ed[3].clone()]),
        t.add_twire(vec![b_ed[0].clone(), b_ed[1].clone(), b_ed[2].clone(), b_ed[3].clone()]),
        t.add_twire(vec![t_ed[0].clone(), t_ed[1].clone(), t_ed[2].clone(), t_ed[3].clone()]),
    ];
    // Surfaces. Profile cap: normal = x_dir × y_dir. Lateral face i: the
    // image of base edge i swept along extr — normal = extr_dir × edge_dir
    // (OCCT BRepPrimAPI_MakePrism, verified on DRAW `prism f1`).
    let pln = |pt: DVec3, n: DVec3, u: DVec3| Surface3::Plane(Plane {
        origin: pt,
        normal: n,
        u_dir: u,
        v_dir: n.cross(u).normalize_or_zero(),
    });
    let cap_n = xd.cross(yd).normalize_or_zero();
    let extr_dir = extr.normalize_or_zero();
    let e_dir = [
        xd,
        yd,
        -xd,
        -yd,
    ];
    let e_pt = [
        vpt(&v[0]),
        vpt(&v[0]),
        vpt(&v[0]),
        vpt(&v[0]),
    ];
    let mut faces: Vec<Shape> = Vec::new();
    for i in 0..4 {
        let n = extr_dir.cross(e_dir[i]).normalize_or_zero();
        let surf = pln(e_pt[i], n, e_dir[i]);
        faces.push(t.add_tface(Some(surf), wires[i].clone(), vec![], None, None, vec![], false));
    }
    // Source cap (FORWARD, outward +cap_n) and extruded cap (FORWARD; the
    // face itself is a new TShape whose wire reuses the located top edges).
    faces.push(t.add_tface(Some(pln(origin, cap_n, xd)), wires[4].clone(), vec![], None, None, vec![], false));
    faces.push(t.add_tface(Some(pln(origin + extr, cap_n, xd)), wires[5].clone(), vec![], None, None, vec![], false));
    let shell = t.add_tshell(faces);
    t.add_tsolid(vec![shell]);
    Ok(t)
}


#[cfg(test)]
mod shell_test {
    use super::*;

    #[test]
    fn prism_shell_closure() {
        let b = make_prism_from_face_brep(
            DVec3::new(4.0, 2.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            1.0, 1.0,
            DVec3::new(-1.0, 0.0, 0.0),
        ).unwrap();
        let mut edge_usage: std::collections::HashMap<usize, Vec<(usize, bool)>> = std::collections::HashMap::new();
        for (fi, s) in b.solids().iter().enumerate() {
            for sh in &s.shells {
                for f in &sh.faces {
                    for we in &f.outer_wire.edges {
                        edge_usage.entry(we.idx).or_default().push((fi, we.forward));
                    }
                }
            }
        }
        let mut bad = 0;
        for (e, uses) in &edge_usage {
            if uses.len() != 2 { println!("edge {} used {} times", e, uses.len()); bad += 1; }
            else if uses[0].1 == uses[1].1 { println!("edge {} same dir on both faces", e); bad += 1; }
        }
        println!("bad edges: {}", bad);
        assert_eq!(bad, 0, "shell must be closed with opposite edge orientations");
    }
}


#[cfg(test)]
mod face_dir_test {
    use super::*;
    use rcad_kernel::topods::TShape;

    #[test]
    fn prism_face_directions() {
        let b = make_prism_from_face_brep(
            DVec3::new(4.0, 2.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
            DVec3::new(0.0, 0.0, 1.0),
            1.0, 1.0,
            DVec3::new(-1.0, 0.0, 0.0),
        ).unwrap();
        for ts in &b.tshapes {
            if let TShape::Face(fd) = ts.as_ref() {
                if let Some(rcad_kernel::geom::Surface3::Plane(p)) = &fd.surface {
                    let mut es = String::new();
                    if let TShape::Wire(wd) = &*fd.outer_wire.data {
                        for e in &wd.edges {
                            es.push_str(&format!(" {}:{:?}", e.ptr_id(), e.orientation));
                        }
                    }
                    println!("face normal=({:.2},{:.2},{:.2}) edges:{}", p.normal.x, p.normal.y, p.normal.z, es);
                }
            }
        }
    }
}
