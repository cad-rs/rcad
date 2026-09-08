//! OCCT BRepAlgo_FaceRestrictor — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepAlgo/
//!         BRepAlgo_FaceRestrictor.cxx (L39-454) +
//!         BRepAlgo_FaceRestrictor.hxx (L36-92).
//!
//! Architecture differences:
//! 1. NCollection_DataMap<TopoDS_Shape, List, ShapeMapHasher> ->
//!    HashMap<ShapeKey, Vec<Shape>> (the key shape is never read back).
//! 2. BRepTopAdaptor_FClass2d -> topalgo::brep_top_adaptor::fclass2d::
//!    FClass2d over FaceShapeSource (loc_ope_wires_on_shape_b.rs arch.
//!    diff. #8 precedent).
//! 3. TopOpeBRepBuild_WireToFace is a GAP-pending re-host (see the struct
//!    comment).
//! 4. BRep_Builder edits are in-place Arc::make_mut mutations (tool.rs);
//!    OCCT TShape sharing makes them visible through every handle, rcad
//!    mutates the owning copy.

use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_first_curve_on_surface, brep_tool_surface,
    builder_add_face_wire, builder_update_edge_pcurve, empty_copied,
    explorer, reversed, shape_is_closed, shape_key, top_exp_vertices_wire, ShapeKey,
};
use rcad_kernel::geom::{Curve2dEval, Curve3, Surface3, TrimmedCurve3};
use rcad_kernel::precision::{CONFUSION, INFINITE_VALUE, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, State};
use std::collections::HashMap;

/// OCCT BRepAlgo_FaceRestrictor (BRepAlgo_FaceRestrictor.hxx L36-92) —
/// builds all the faces limited with a set of non jointing and planars
/// wires.
pub struct BRepAlgoFaceRestrictor {
    /// OCCT: myDone.
    my_done: bool,
    /// OCCT: modeProj.
    mode_proj: bool,
    /// OCCT: myFace.
    my_face: Shape,
    /// OCCT: wires.
    wires: Vec<Shape>,
    /// OCCT: faces.
    faces: Vec<Shape>,
    /// OCCT: myCorrection.
    my_correction: bool,
    /// OCCT: keyIsIn.
    key_is_in: HashMap<ShapeKey, Vec<Shape>>,
    /// OCCT: keyContains.
    key_contains: HashMap<ShapeKey, Vec<Shape>>,
}

impl Default for BRepAlgoFaceRestrictor {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAlgoFaceRestrictor {
    /// OCCT BRepAlgo_FaceRestrictor::BRepAlgo_FaceRestrictor() (cxx L39).
    pub fn new() -> Self {
        BRepAlgoFaceRestrictor {
            my_done: false,
            mode_proj: false,
            my_face: Shape::null(),
            wires: Vec::new(),
            faces: Vec::new(),
            my_correction: false,
            key_is_in: HashMap::new(),
            key_contains: HashMap::new(),
        }
    }

    /// OCCT BRepAlgo_FaceRestrictor::Init(F, Proj, CorrectionOrientation)
    /// (cxx L43-50) — the surface of F will be the surface of each new face
    /// built; Proj is used to update pcurves on edges if necessary.
    pub fn init(&mut self, f: &Shape, proj: bool, correction_orientation: bool) {
        self.my_face = f.clone();
        self.mode_proj = proj;
        self.my_correction = correction_orientation;
    }

    /// OCCT BRepAlgo_FaceRestrictor::Add(W) (cxx L54-57) — adds the wire W
    /// to the set of wires.  The wires must be closed.
    pub fn add(&mut self, w: &Shape) {
        self.wires.push(w.clone());
    }

    /// OCCT BRepAlgo_FaceRestrictor::Clear() (cxx L61-65) — removes all the
    /// wires.
    pub fn clear(&mut self) {
        self.wires.clear();
        self.faces.clear();
    }

    /// OCCT BRepAlgo_FaceRestrictor::Perform() (cxx L109-179) — evaluates
    /// all the faces limited by the set of wires.
    pub fn perform(&mut self) {
        if self.my_correction {
            self.perform_with_correction();
            return;
        }

        self.my_done = false;

        //--------------------------------------------------------------------
        // return geometry of the reference face (OCCT L121-125).
        //--------------------------------------------------------------------
        let s = brep_tool_surface(&self.my_face);
        let Some(s) = s else {
            // OCCT would carry the null surface handle into the pcurve
            // queries below (all of which then fail); rcad returns (marked).
            return;
        };

        //-----------------------------------------------------------------------
        // test if edges are on S. otherwise add S to the first pcurve.
        // or projection of the edge on F (OCCT L127-135).
        //----------------------------------------------------------------------
        let mut wtf = TopOpeBRepBuildWireToFace::new();

        for it in &self.wires {
            // update the surface on edges.
            let w = it; // TopoDS::Wire(it.Value())

            for e in explorer(w, ShapeType::Edge, ShapeType::Shape) {
                let mut e = e;
                let c2 = brep_tool_curve_on_surface(&e, &self.my_face);

                if c2.is_none() {
                    // no pcurve on the reference surface (OCCT L148-171).
                    if self.mode_proj {
                        // Projection of the 3D curve on surface.
                        if !proj_curve_3d(&mut e, &s, &self.my_face) {
                            return;
                        }
                    } else {
                        // return the first pcurve glued on <S>.
                        let ya_pcurve = change_pcurve(&mut e, &s, &self.my_face);
                        if !ya_pcurve {
                            if !proj_curve_3d(&mut e, &s, &self.my_face) {
                                return;
                            }
                        }
                    }
                }
            }
            wtf.add_wire(w);
        }

        wtf.make_faces(&self.my_face, &mut self.faces);

        self.my_done = true;
    }

    /// OCCT BRepAlgo_FaceRestrictor::IsDone() const (cxx L183-186).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepAlgo_FaceRestrictor::More() const (cxx L190-193).
    pub fn more(&self) -> bool {
        !self.faces.is_empty()
    }

    /// OCCT BRepAlgo_FaceRestrictor::Next() (cxx L197-200).
    pub fn next(&mut self) {
        if !self.faces.is_empty() {
            self.faces.remove(0);
        }
    }

    /// OCCT BRepAlgo_FaceRestrictor::Current() const (cxx L204-207).
    pub fn current(&self) -> Shape {
        self.faces.first().expect("faces.First()").clone()
    }

    /// OCCT BRepAlgo_FaceRestrictor::PerformWithCorrection() (cxx L370-454)
    /// — evaluates all the faces limited by the set of wires.
    fn perform_with_correction(&mut self) {
        self.my_done = false;
        //---------------------------------------------------------
        // Reorientation of all closed wires to the left (OCCT L376-396).
        //---------------------------------------------------------
        for i in 0..self.wires.len() {
            let w = &mut self.wires[i]; // TopoDS::Wire(it.ChangeValue())
            let mut nf = empty_copied(&self.my_face);
            nf.orientation = Orientation::Forward;
            builder_add_face_wire(&mut nf, w);

            if is_closed(w) {
                let Some(surf) = brep_tool_surface(&nf) else {
                    continue;
                };
                let src = crate::topalgo::shape_source::FaceShapeSource::new(
                    &nf,
                    surf,
                    &[glam::DAffine3::IDENTITY],
                );
                let fclass2d = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
                    &src,
                    0,
                    PCONFUSION,
                );
                if fclass2d.perform_infinite_point(&src) != State::Out {
                    *w = reversed(w);
                }
            }
        }
        //---------------------------------------------------------
        // Classification of wires ones compared to the others (OCCT
        // L397-424).
        //---------------------------------------------------------
        for i1 in 0..self.wires.len() {
            let w1 = self.wires[i1].clone(); // TopoDS::Wire(it.Value())

            if is_closed(&w1) {
                let Some(surf) = brep_tool_surface(&self.my_face) else {
                    continue;
                };
                let mut nf = empty_copied(&self.my_face);
                nf.orientation = Orientation::Forward;
                builder_add_face_wire(&mut nf, &w1);

                let src = crate::topalgo::shape_source::FaceShapeSource::new(
                    &nf,
                    surf,
                    &[glam::DAffine3::IDENTITY],
                );
                let fclass2d = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
                    &src,
                    0,
                    PCONFUSION,
                );
                for i2 in 0..self.wires.len() {
                    let w2 = &self.wires[i2];
                    if !w1.is_same(w2) && is_inside(w2, &nf, &src, &fclass2d) {
                        store(w2, &w1, &mut self.key_is_in, &mut self.key_contains);
                    }
                }
            }
        }
        let mut wire_ext: Vec<Shape> = Vec::new();

        for it in &self.wires {
            let w = it; // TopoDS::Wire(it.Value())
            if !self.key_is_in.contains_key(&shape_key(w))
                || self.key_is_in.get(&shape_key(w)).map_or(true, |l| l.is_empty())
            {
                wire_ext.push(w.clone());
            }
        }

        for it in &wire_ext {
            let w = it; // TopoDS::Wire(it.Value())
            if !self.key_is_in.contains_key(&shape_key(w))
                || self.key_is_in.get(&shape_key(w)).map_or(true, |l| l.is_empty())
            {
                let mut new_face = empty_copied(&self.my_face);
                new_face.orientation = Orientation::Forward;
                builder_add_face_wire(&mut new_face, w);
                self.faces.push(new_face.clone());
                //--------------------------------------------
                // Construction of a face by exterior wire (OCCT L447-451).
                //--------------------------------------------
                build_face_in(
                    &mut new_face,
                    w,
                    &mut self.key_contains,
                    &mut self.key_is_in,
                    Orientation::Forward,
                    &mut self.faces,
                );
            }
        }
        self.my_done = true;
    }
}

/// OCCT TopOpeBRepBuild_WireToFace (TopOpeBRepBuild_WireToFace.cxx L32-64) —
/// the AddWire/MakeFaces vehicle of the non-correction Perform path.
///
/// GAP (architecture difference #3): MakeFaces drives the TKBool
/// TopOpeBRepBuild chain (WireEdgeSet / FaceBuilder with ForceClass=true /
/// TopOpeBRepBuild_Builder::MakeFaces) which is not translated yet; the
/// output list stays empty so the callers take the empty-faces path.  Closes
/// with the TKBool/TopOpeBRepBuild batch.
pub struct TopOpeBRepBuildWireToFace {
    /// OCCT: myLW.
    my_lw: Vec<Shape>,
}

impl Default for TopOpeBRepBuildWireToFace {
    fn default() -> Self {
        Self::new()
    }
}

impl TopOpeBRepBuildWireToFace {
    /// OCCT TopOpeBRepBuild_WireToFace::TopOpeBRepBuild_WireToFace()
    /// (cxx L34-35).
    pub fn new() -> Self {
        TopOpeBRepBuildWireToFace { my_lw: Vec::new() }
    }

    /// OCCT TopOpeBRepBuild_WireToFace::Init() (cxx L38-41).
    pub fn init(&mut self) {
        self.my_lw.clear();
    }

    /// OCCT TopOpeBRepBuild_WireToFace::AddWire(W) (cxx L44-47).
    pub fn add_wire(&mut self, w: &Shape) {
        self.my_lw.push(w.clone());
    }

    /// OCCT TopOpeBRepBuild_WireToFace::MakeFaces(F, LF) (cxx L50-64).
    pub fn make_faces(&mut self, _f: &Shape, lf: &mut Vec<Shape>) {
        lf.clear();

        // OCCT L55-62:
        //   TopOpeBRepBuild_WireEdgeSet wes(F);
        //   for (it(myLW)) wes.AddShape(it.Value());
        //   bool ForceClass = true;
        //   TopOpeBRepBuild_FaceBuilder FB;
        //   FB.InitFaceBuilder(wes, F, ForceClass);
        //   TopOpeBRepDS_BuildTool BT(TopOpeBRepTool_APPROX);
        //   TopOpeBRepBuild_Builder B(BT);
        //   B.MakeFaces(F, FB, LF);
        // (pending — see the struct comment.)
    }
}

/// OCCT static ChangePCurve (cxx L69-83).
///
/// OCCT L77: BRep_Tool::CurveOnSurface(E, C2, SE, LE, f, l, 1) — the FIRST
/// pcurve whatever its surface; OCCT L80 binds it on (S, L) through
/// BRep_Builder::UpdateEdge.  The rcad pcurve map is keyed by the face
/// shape, so the binding goes through the reference face theFace (the
/// surface carrier of this pipeline — architecture difference, the
/// Geom_Surface handle identity maps to the face shape key).
fn change_pcurve(e: &mut Shape, s: &Surface3, the_face: &Shape) -> bool {
    let c2 = brep_tool_first_curve_on_surface(e);
    let _ = s;
    match &c2 {
        Some((c2d, _, _)) => {
            builder_update_edge_pcurve(e, c2d, the_face, CONFUSION);
        }
        None => {}
    }
    c2.is_none()
}

/// OCCT static ProjCurve3d (cxx L87-105).
///
/// OCCT L103 binds the projected pcurve on (S, L) through
/// BRep_Builder::UpdateEdge; the rcad pcurve map is keyed by the face shape,
/// so the binding goes through the reference face theFace (architecture
/// difference, same as change_pcurve).
fn proj_curve_3d(e: &mut Shape, s: &Surface3, the_face: &Shape) -> bool {
    // OCCT L92: C = BRep_Tool::Curve(E, LE, f, l).
    let Some((c, f, l)) = crate::brep_algo::tool::brep_tool_curve(e) else {
        return false;
    };
    // OCCT L97: CT = new Geom_TrimmedCurve(C, f, l).
    let _ct = Curve3::Trimmed(TrimmedCurve3::new(c.clone(), f, l));
    // OCCT L99-100: LL = L.Inverted().Multiplied(LE); CT->Transform(...) —
    // locations are identity in this pipeline (arch. diff. #1), no-op.
    // OCCT L102: C2 = GeomProjLib::Curve2d(CT, S).
    let Some(c2) = rcad_kernel::base::geom_proj_lib::curve2d_simple(&c, f, l, s) else {
        // OCCT binds the null C2 through UpdateEdge (the null handle is
        // ignored on use); rcad reports the null binding the same way.
        return true;
    };
    builder_update_edge_pcurve(e, &c2, the_face, CONFUSION);
    true
}

/// OCCT static IsClosed (cxx L211-221).
fn is_closed(w: &Shape) -> bool {
    if shape_is_closed(w) {
        return true;
    }
    let (v1, v2) = top_exp_vertices_wire(w);
    match (&v1, &v2) {
        (Some(a), Some(b)) => a.is_same(b),
        // OCCT: V1.IsSame(V2) with null vertices is true; both-null ->
        // closed is not reachable for wires with edges (the OCCT null
        // handles compare equal).
        (None, None) => true,
        _ => false,
    }
}

/// OCCT static IsInside (cxx L225-264).
fn is_inside(
    wir: &Shape,
    f: &Shape,
    src: &crate::topalgo::shape_source::FaceShapeSource,
    _fclass2d: &crate::topalgo::brep_top_adaptor::fclass2d::FClass2d,
) -> bool {
    let mut exp = explorer(wir, ShapeType::Edge, ShapeType::Shape);
    if !exp.is_empty() {
        let edg = exp.remove(0); // TopoDS::Edge(exp.Current())
        let Some((c2d, f, l)) = brep_tool_curve_on_surface(&edg, f) else {
            // OCCT dereferences the null C2d (cxx L258); rcad returns false
            // (marked).
            return false;
        };
        let prm: f64;

        if !is_negative_infinite(f) && !is_positive_infinite(l) {
            prm = (f + l) / 2.0;
        } else if is_negative_infinite(f) && is_positive_infinite(l) {
            prm = 0.0;
        } else if is_negative_infinite(f) {
            prm = l - 1.0;
        } else {
            prm = f + 1.0;
        }

        let pt2d = c2d.point_at(prm);
        // OCCT L259: a local BRepTopAdaptor_FClass2d(F, PConfusion()) is
        // built here (the passed-in one is ignored).
        let f_class2d = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(src, 0, PCONFUSION);
        let st2 = f_class2d.perform(src, pt2d, false);
        return st2 == State::In;
    }
    false
}

/// OCCT Precision::IsNegativeInfinite.
fn is_negative_infinite(v: f64) -> bool {
    v <= -INFINITE_VALUE
}

/// OCCT Precision::IsPositiveInfinite.
fn is_positive_infinite(v: f64) -> bool {
    v >= INFINITE_VALUE
}

/// OCCT static Store (cxx L268-288).
fn store(
    w2: &Shape,
    w1: &Shape,
    key_is_in: &mut HashMap<ShapeKey, Vec<Shape>>,
    key_contains: &mut HashMap<ShapeKey, Vec<Shape>>,
) {
    let k2 = shape_key(w2);
    if !key_is_in.contains_key(&k2) {
        key_is_in.insert(k2, Vec::new());
    }
    key_is_in.get_mut(&k2).expect("keyIsIn(W2)").push(w1.clone());
    let k1 = shape_key(w1);
    if !key_contains.contains_key(&k1) {
        key_contains.insert(k1, Vec::new());
    }
    key_contains.get_mut(&k1).expect("keyContains(W1)").push(w2.clone());
}

/// OCCT static BuildFaceIn (cxx L292-366).
fn build_face_in(
    f: &mut Shape,
    w: &Shape,
    key_contains: &mut HashMap<ShapeKey, Vec<Shape>>,
    key_is_in: &mut HashMap<ShapeKey, Vec<Shape>>,
    orientation: Orientation,
    faces: &mut Vec<Shape>,
) {
    if !key_contains.contains_key(&shape_key(w))
        || key_contains.get(&shape_key(w)).map_or(true, |l| l.is_empty())
    {
        return;
    }

    // Removal of W in KeyIsIn (OCCT L309-326).
    let contains_list = key_contains.get(&shape_key(w)).cloned().unwrap_or_default();
    for it in &contains_list {
        let wi = it; // TopoDS::Wire(it.Value())
        let kwi = shape_key(wi);
        if let Some(l2) = key_is_in.get_mut(&kwi) {
            if let Some(pos) = l2.iter().position(|x| x.is_same(w)) {
                l2.remove(pos);
            }
        }
    }

    let mut wire_ext: Vec<Shape> = Vec::new();

    for it in &contains_list {
        let wi = it; // TopoDS::Wire(it.Value())
        let l2_empty = key_is_in
            .get(&shape_key(wi))
            .map_or(true, |l| l.is_empty());
        if l2_empty {
            wire_ext.push(wi.clone());
        }
    }

    for it in &wire_ext {
        let wi = it; // TopoDS::Wire(it.Value())
        let l2_empty = key_is_in
            .get(&shape_key(wi))
            .map_or(true, |l| l.is_empty());
        if l2_empty {
            if orientation == Orientation::Forward {
                let nwi = reversed(wi);
                // OCCT L349-351: NWI.Reverse() — TopAbs::Reverse of WI.
                builder_add_face_wire(f, &nwi);
                build_face_in(f, wi, key_contains, key_is_in, Orientation::Reversed, faces);
            } else {
                // OCCT L357-362: NF = TopoDS::Face(Faces.First().EmptyCopied()).
                let first = faces.first().expect("Faces.First()").clone();
                let mut nf = empty_copied(&first);
                builder_add_face_wire(&mut nf, wi);
                faces.push(nf.clone());
                build_face_in(&mut nf, wi, key_contains, key_is_in, Orientation::Forward, faces);
            }
        }
    }
}
