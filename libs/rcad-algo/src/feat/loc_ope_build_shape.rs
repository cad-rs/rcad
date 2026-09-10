// OCCT LocOpe_BuildShape.hxx L26-45 + LocOpe_BuildShape.cxx L40-465 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildShape.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildShape.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_BuildShape.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. BRep_Builder (TopoDS_Builder) — rcad's BRepBuilder
//    (rcad_kernel::topods) carries the OCCT MakeXxx/Add surface; the
//    container is an owned rcad topods::BRep pool local to Perform and the
//    resulting Shape references it (OCCT: TopoDS compound + incremental
//    Add). Same vehicle as the rcad ResultBuilder pattern (AGENTS.md A1).
// 2. BRep_Tool::IsClosed(shell) (BRep_Tool.cxx) is re-hosted below
//    (brep_tool_is_closed): the CLOSED flag short-circuit plus the
//    every-edge-shared-twice walk.
// 3. NCollection_IndexedMap / NCollection_IndexedDataMap of shapes map to
//    insertion-ordered IndexMaps keyed by (TShape ptr, Location) — the
//    TopTools_ShapeMapHasher identity, the same model as bop/algo/section.rs.
//    NCollection_Map<int> maps to a HashSet<i32>; the OCCT bucket iteration
//    order is not reproduced (Propagate iterates mapIf — the order
//    sensitivity across multiple candidate faces follows the HashSet order).
// 4. TopExp_Explorer / TopExp::MapShapesAndAncestors are the re-hosts of
//    feat::brep_feat_builder::explorer and loc_ope_glued_shape's
//    map_shapes_and_ancestors (same crate).
// 5. BRepClass3d_SolidClassifier is rcad's
//    topalgo::brep_class3d::SolidClassifier (state mapping: 0 = TopAbs_IN,
//    1 = TopAbs_OUT, 2 = TopAbs_ON).
// 6. The OCCT dead local `mapSh` (cxx L88, declared but never used) is not
//    translated.
//
// first consumer: BRepFeat_Form family (3b) — LocOpe_Generator /
// LocOpe_Gluer::Perform (cxx L227) build the resulting shape through
// LocOpe_BuildShape.

use crate::feat::brep_feat_builder::explorer;
use indexmap::IndexMap;
use rcad_kernel::topo::topods::{ tshape_flags, BRep, BRepBuilder, TShape };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, State };
use std::collections::HashSet;
use std::sync::Arc;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopAbs::Reverse (TopAbs.hxx) — FORWARD<->REVERSED, INTERNAL/
/// EXTERNAL unchanged.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// OCCT TopoDS_Shape::Reversed().
fn shape_reversed(s: &Shape) -> Shape {
    let mut c = s.clone();
    c.orientation = top_abs_reverse(c.orientation);
    c
}

/// OCCT BRep_Tool::Pnt(vtx).
fn brep_tool_pnt(vtx: &Shape) -> glam::DVec3 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(vtx).
fn brep_tool_tolerance(vtx: &Shape) -> f64 {
    match vtx.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::IsClosed(theShape) re-host (architecture difference #2).
fn brep_tool_is_closed(the_shape: &Shape) -> bool {
    // OCCT: the CLOSED flag short-circuits.
    let closed_flag = match the_shape.data.as_ref() {
        TShape::Wire(wd) => wd.flags & tshape_flags::CLOSED != 0,
        TShape::Shell(sd) => sd.flags & tshape_flags::CLOSED != 0,
        TShape::Solid(sd) => sd.flags & tshape_flags::CLOSED != 0,
        _ => false,
    };
    if closed_flag {
        return true;
    }
    // OCCT: otherwise every edge must be shared exactly twice.
    let mut counts: std::collections::HashMap<(u64, u32), i32> = std::collections::HashMap::new();
    for edg in explorer(the_shape, ShapeType::Edge, ShapeType::Shape) {
        *counts.entry(shape_key(&edg)).or_insert(0) += 1;
    }
    counts.values().all(|&c| c == 2)
}

/// OCCT TopoDS_Shape::Closed(theFlag) — set/clear the CLOSED bit on the
/// TShape (TopoDS_Shape.cxx).
fn shape_set_closed(pool: &mut BRep, s: &Shape, the_flag: bool) {
    let idx = s.index;
    if idx >= pool.tshapes.len() {
        return;
    }
    let ts = Arc::make_mut(&mut pool.tshapes[idx]);
    let flags = match ts {
        TShape::Wire(wd) => &mut wd.flags,
        TShape::Shell(sd) => &mut sd.flags,
        TShape::Solid(sd) => &mut sd.flags,
        _ => return,
    };
    if the_flag {
        *flags |= tshape_flags::CLOSED;
    } else {
        *flags &= !tshape_flags::CLOSED;
    }
}

/// OCCT TopoDS_Shape::Oriented(theOr) — a copy carrying the orientation.
fn with_orientation(s: &Shape, the_or: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = the_or;
    c
}

/// OCCT static Add(ind, mapI, mapF, mapEF) (cxx L324-362).
fn add(
    the_ind: i32,
    the_map_i: &mut HashSet<i32>,
    the_map_f: &mut IndexMap<(u64, u32), Shape>,
    the_map_ef: &IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
) {
    // OCCT cxx L332-335: if (!mapI.Add(ind)) throw.
    if !the_map_i.insert(the_ind) {
        panic!("Standard_ConstructionError");
    }

    // OCCT cxx L337-361.
    let (_, faces) = the_map_ef
        .get_index((the_ind - 1) as usize)
        .expect("mapEF entry")
        .1
        .clone();
    for f in faces.iter() {
        if !the_map_f.contains_key(&shape_key(f)) {
            the_map_f.insert(shape_key(f), f.clone());
            // OCCT cxx L343-359: edges of the face.
            for edg in explorer(f, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L350-354: FindIndex(edg) == 0 -> throw.
                let Some(ind0) = the_map_ef.get_index_of(&shape_key(&edg)) else {
                    panic!("Standard_ConstructionError");
                };
                // OCCT cxx L355-358: recursion on unvisited edges.
                if !the_map_i.contains(&((ind0 + 1) as i32)) {
                    add((ind0 + 1) as i32, the_map_i, the_map_f, the_map_ef);
                }
            }
        }
    }
}

/// OCCT static Propagate(F, Sh, mapF, mapIf) (cxx L366-436). The OCCT shell
/// argument is a TopoDS_Shape built incrementally through the builder; rcad
/// collects the faces into the_vec (architecture difference #1) and the
/// caller Adds them in order.
fn propagate(
    the_f: &Shape,
    the_vec: &mut Vec<Shape>,
    the_map_f: &IndexMap<(u64, u32), Shape>,
    the_map_if: &mut HashSet<i32>,
) {
    // OCCT cxx L372-377: indf = mapF.FindIndex(F); if (!mapIf.Contains)
    // return; mapIf.Remove(indf).
    let Some(indf) = the_map_f.get_index_of(&shape_key(the_f)) else {
        return;
    };
    if !the_map_if.contains(&((indf + 1) as i32)) {
        return;
    }
    the_map_if.remove(&((indf + 1) as i32));
    // OCCT cxx L378-381: if (mapIf.Extent() == 0) return.
    if the_map_if.is_empty() {
        return;
    }

    // OCCT cxx L383-435.
    for edg in explorer(the_f, ShapeType::Edge, ShapeType::Shape) {
        // OCCT cxx L389.
        let ored1 = edg.orientation;

        // OCCT cxx L391-394.
        if ored1 == Orientation::Internal || ored1 == Orientation::External {
            continue;
        }
        // OCCT cxx L395-414: find the mapIf face sharing the edge (the OCCT
        // NCollection_Map<int> iterator order is not reproduced — arch.
        // diff. #3).
        let mut ored2 = Orientation::Forward;
        let mut found_key: Option<i32> = None;
        for &k in the_map_if.iter() {
            let new_f = &the_map_f.get_index((k - 1) as usize).expect("mapF entry").1;
            // OCCT cxx L400-408: explorer edges of newF.
            let mut hit: Option<Orientation> = None;
            for e2 in explorer(new_f, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L404: exp2.Current().IsSame(edg).
                if shape_key(&e2) == shape_key(&edg) {
                    hit = Some(e2.orientation);
                    break;
                }
            }
            // OCCT cxx L409-413: if (exp2.More()) { ored2 = ...; break; }
            if let Some(o2) = hit {
                ored2 = o2;
                found_key = Some(k);
                break;
            }
        }
        // OCCT cxx L415-434: if (itm.More()).
        if let Some(k) = found_key {
            let f_to_add = the_map_f
                .get_index((k - 1) as usize)
                .expect("mapF entry")
                .1
                .clone();
            let mut added = false;
            if ored2 == ored1 {
                // OCCT cxx L420-423: FtoAdd.Reverse(); B.Add(Sh, FtoAdd).
                let f_add = shape_reversed(&f_to_add);
                the_vec.push(f_add);
                added = true;
            } else if ored2 == top_abs_reverse(ored1) {
                // OCCT cxx L425-428: B.Add(Sh, FtoAdd).
                the_vec.push(f_to_add.clone());
                added = true;
            }
            // OCCT cxx L430-433: if (added) Propagate(FtoAdd, Sh, ...).
            if added {
                propagate(&f_to_add, the_vec, the_map_f, the_map_if);
            }
        }
    }
}

/// OCCT static IsInside(S1, S2) (cxx L440-465).
fn is_inside(the_s1: &Shape, the_s2: &Shape) -> bool {
    // OCCT cxx L442: BRepClass3d_SolidClassifier Class(S2) (arch. diff. #5).
    let mut class_ = crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(the_s2);
    // OCCT cxx L443-459.
    for vtx in explorer(the_s1, ShapeType::Vertex, ShapeType::Shape) {
        let pttest = brep_tool_pnt(&vtx);
        let tol = brep_tool_tolerance(&vtx);
        class_.perform(pttest, tol);
        let state = match class_.state() {
            0 => State::In,
            1 => State::Out,
            2 => State::On,
            _ => State::Out,
        };
        if state == State::In {
            return true;
        } else if state == State::Out {
            return false;
        }
    }
    // OCCT cxx L460-464: "Classification impossible sur vertex" -> true.
    true
}

/// OCCT LocOpe_BuildShape (LocOpe_BuildShape.hxx L26-45).
pub struct LocOpeBuildShape {
    my_res: Option<Shape>, // OCCT: myRes (TopoDS_Shape; None = null)
}

impl Default for LocOpeBuildShape {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeBuildShape {
    /// OCCT LocOpe_BuildShape::LocOpe_BuildShape() (hxx L31).
    pub fn new() -> Self {
        LocOpeBuildShape { my_res: None }
    }

    /// OCCT LocOpe_BuildShape::LocOpe_BuildShape(L) (hxx L35) — builds
    /// shape(s) from the list <the_l>. Uses only the faces of <the_l>.
    pub fn with_faces(the_l: &[Shape]) -> Self {
        let mut res = LocOpeBuildShape::new();
        res.perform(the_l);
        res
    }

    /// OCCT LocOpe_BuildShape::Perform(L) (cxx L56-320).
    pub fn perform(&mut self, the_l: &[Shape]) {
        // OCCT cxx L58-61: i/j/k declared; myRes.Nullify().
        self.my_res = None;

        // OCCT cxx L63-65: compound C + builder B (architecture difference
        // #1 — the pool carries the container; C feeds MapShapesAndAncestors
        // and B.Add appends into it).
        let mut pool = BRep::new();
        let b = BRepBuilder::new();
        let mut b = b;
        let c = b.make_compound(&mut pool, Vec::new());

        // OCCT cxx L67-77: mapF + the face fill.
        let mut map_f: IndexMap<(u64, u32), Shape> = IndexMap::new();
        for s in the_l {
            if s.shape_type() == ShapeType::Face && !map_f.contains_key(&shape_key(s)) {
                map_f.insert(shape_key(s), s.clone());
                // OCCT cxx L75: B.Add(C, itl.Value()).
                b.add_to_compound(&mut pool, c.clone(), s.clone());
            }
        }

        // OCCT cxx L79-82: if (mapF.Extent() == 0) return. — no face.
        if map_f.is_empty() {
            return;
        }

        // OCCT cxx L84-86: theMapEF.
        let mut the_map_ef: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &c,
            ShapeType::Edge,
            ShapeType::Face,
            &mut the_map_ef,
        );

        // OCCT cxx L88-90: mapI, mapIf, Nbedges (the dead mapSh local of the
        // OCCT source is not translated — arch. diff. #6).
        let mut map_i: HashSet<i32> = HashSet::new();
        let mut map_if: HashSet<i32> = HashSet::new();
        let nbedges = the_map_ef.len() as i32;

        // OCCT cxx L92-93.
        let mut lshell: Vec<Shape> = Vec::new();
        let mut lresult: Vec<Shape> = Vec::new();

        // OCCT cxx L95-211: do { ... } while (mapI.Extent() < Nbedges).
        loop {
            // OCCT cxx L97-104: first unprocessed edge.
            let mut i = 1;
            while i <= nbedges {
                if !map_i.contains(&i) {
                    break;
                }
                i += 1;
            }
            // OCCT cxx L105: if (i <= Nbedges).
            if i <= nbedges {
                // OCCT cxx L107-109.
                map_f.clear();
                map_if.clear();
                add(i, &mut map_i, &mut map_f, &the_map_ef);
                let mut manifold = true;
                let mut face_ref: Option<Shape> = None;

                // OCCT cxx L114-126.
                for j in 1..=(map_f.len() as i32) {
                    let (_, f) = map_f.get_index((j - 1) as usize).expect("mapF entry");
                    let orient = f.orientation;
                    if orient == Orientation::Internal || orient == Orientation::External {
                        manifold = false;
                    } else if face_ref.is_none() {
                        face_ref = Some(f.clone());
                    }
                    map_if.insert(j);
                }

                // OCCT cxx L128-129: newSh.
                let new_sh = b.make_shell(&mut pool);

                // OCCT cxx L130-156.
                if !manifold && face_ref.is_none() {
                    // "on a un paquet de faces. pas d'orientation possible ?"
                    for j in 1..=(map_f.len() as i32) {
                        let (_, f) = map_f.get_index((j - 1) as usize).expect("mapF entry");
                        // OCCT cxx L135: B.Add(newSh, mapF(j)).
                        b.add_to_shell(&mut pool, new_sh.clone(), f.clone());
                    }
                } else {
                    // "orienter ce qu'on peut"
                    // OCCT cxx L141-152.
                    if !manifold {
                        for j in 1..=(map_f.len() as i32) {
                            let (_, f) = map_f.get_index((j - 1) as usize).expect("mapF entry");
                            if f.orientation == Orientation::Internal
                                || f.orientation == Orientation::External
                            {
                                // OCCT cxx L148-149.
                                b.add_to_shell(&mut pool, new_sh.clone(), f.clone());
                                map_if.remove(&j);
                            }
                        }
                    }

                    // OCCT cxx L154-155: B.Add(newSh, FaceRef); Propagate(...).
                    let f_ref = face_ref.expect("FaceRef non-null on this path");
                    b.add_to_shell(&mut pool, new_sh.clone(), f_ref.clone());
                    let mut propagated: Vec<Shape> = Vec::new();
                    propagate(&f_ref, &mut propagated, &map_f, &mut map_if);
                    for f in &propagated {
                        b.add_to_shell(&mut pool, new_sh.clone(), f.clone());
                    }
                }

                // OCCT cxx L157: newSh.Closed(BRep_Tool::IsClosed(newSh)).
                shape_set_closed(&mut pool, &new_sh, brep_tool_is_closed(&new_sh));

                // OCCT cxx L158-209.
                if !manifold {
                    // OCCT cxx L160.
                    lshell.push(with_orientation(&new_sh, Orientation::Internal));
                } else {
                    // OCCT cxx L164-168: theMapEFbis.
                    let mut the_map_ef_bis: IndexMap<(u64, u32), (Shape, Vec<Shape>)> =
                        IndexMap::new();
                    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
                        &new_sh,
                        ShapeType::Edge,
                        ShapeType::Face,
                        &mut the_map_ef_bis,
                    );
                    // OCCT cxx L169-188.
                    let mut k = 1usize;
                    while k <= the_map_ef_bis.len() {
                        let (ed, list) = the_map_ef_bis.get_index(k - 1).expect("mapEFbis entry");
                        let ori_ed = list.0.orientation;
                        if ori_ed != Orientation::Internal && ori_ed != Orientation::External {
                            let nb_fac = list.1.len();
                            if nb_fac > 2 {
                                // OCCT cxx L176-179: "peu probable" -> break.
                                break;
                            } else if nb_fac == 1 {
                                // OCCT cxx L180-186.
                                if !brep_tool_degenerated(&list.0) {
                                    break;
                                }
                            }
                        }
                        k += 1;
                    }
                    // OCCT cxx L189-208.
                    if k > the_map_ef_bis.len() {
                        // OCCT cxx L191-193: newSo solid with the FORWARD shell.
                        let new_so = b.make_solid(&mut pool, vec![new_sh.clone()]);
                        // OCCT cxx L194-195.
                        let mut class_ =
                            crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(&new_so);
                        class_.perform_infinite_point(rcad_kernel::precision::CONFUSION);
                        // OCCT cxx L196-203: state mapping (arch. diff. #5).
                        if class_.state() == 0 {
                            // TopAbs_IN
                            lresult.push(with_orientation(&new_sh, Orientation::Reversed));
                        } else {
                            lresult.push(new_sh.clone());
                        }
                    } else {
                        // OCCT cxx L207.
                        lshell.push(with_orientation(&new_sh, Orientation::Internal));
                    }
                }
            }
            // OCCT cxx L211: while (mapI.Extent() < Nbedges).
            if !((map_i.len() as i32) < nbedges) {
                break;
            }
        }

        // on a une list de shells dans lresult. on suppose qu'ils ne
        // s'intersectent pas. il faut classifier les shells orientes pour en
        // faire des solides... on n'accepte qu'1 niveau d'imbrication

        // OCCT cxx L217-245: imbSh + LIntern.
        let mut imb_sh: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        let mut l_intern: Vec<Shape> = Vec::new();

        for sh in &lresult {
            // OCCT cxx L223-225: tempo solid around sh.
            let tempo = b.make_solid(&mut pool, vec![sh.clone()]);

            imb_sh.insert(shape_key(sh), (sh.clone(), Vec::new()));
            for sh2 in &lresult {
                if shape_key(sh2) != shape_key(sh) {
                    // OCCT cxx L238: IsInside(sh2, tempo).
                    if is_inside(sh2, &tempo) {
                        l_intern.push(sh2.clone());
                        imb_sh
                            .get_mut(&shape_key(sh))
                            .expect("imb entry")
                            .1
                            .push(sh2.clone());
                    }
                }
            }
        }

        // OCCT cxx L247-256: "on vire les shells imbriques comme etant aussi
        // des solides a part entiere."
        for sh in &l_intern {
            imb_sh.shift_remove(&shape_key(sh));
        }

        // OCCT cxx L258-295: lsolid assembly.
        let mut lsolid: Vec<Shape> = Vec::new();
        loop {
            // OCCT cxx L263-271: first entry with a non-empty list.
            let mut hit: Option<(u64, u32)> = None;
            for (key, (_, list)) in imb_sh.iter() {
                if !list.is_empty() {
                    hit = Some(*key);
                    break;
                }
            }
            // OCCT cxx L272-283.
            if let Some(key) = hit {
                let (outer, inners) = imb_sh.get(&key).expect("imb entry").clone();
                let mut shells = vec![outer.clone()];
                for inner in &inners {
                    // OCCT cxx L279: B.Add(newSo, itl.Value().Reversed()).
                    shells.push(shape_reversed(inner));
                }
                let new_so = b.make_solid(&mut pool, shells);
                lsolid.push(new_so);
                imb_sh.shift_remove(&key);
            } else {
                // OCCT cxx L284-294.
                let keys: Vec<(u64, u32)> = imb_sh.keys().copied().collect();
                for key in keys {
                    let (outer, _) = imb_sh.get(&key).expect("imb entry").clone();
                    let new_so = b.make_solid(&mut pool, vec![outer]);
                    lsolid.push(new_so);
                }
                imb_sh.clear();
            }
            // OCCT cxx L295: while (imbSh.Extent() != 0).
            if imb_sh.is_empty() {
                break;
            }
        }

        // OCCT cxx L297-319.
        let nbsol = lsolid.len() as i32;
        let nbshl = lshell.len() as i32;

        if nbsol == 1 && nbshl == 0 {
            // OCCT cxx L302.
            self.my_res = Some(lsolid.first().expect("lsolid").clone());
        } else if nbsol == 0 && nbshl == 1 {
            // OCCT cxx L306.
            self.my_res = Some(lshell.first().expect("lshell").clone());
        } else {
            // OCCT cxx L310: B.MakeCompound(TopoDS::Compound(myRes)).
            let res = b.make_compound(&mut pool, Vec::new());
            for s in &lsolid {
                // OCCT cxx L313.
                b.add_to_compound(&mut pool, res.clone(), s.clone());
            }
            for s in &lshell {
                // OCCT cxx L317.
                b.add_to_compound(&mut pool, res.clone(), s.clone());
            }
            self.my_res = Some(res);
        }
    }

    /// OCCT LocOpe_BuildShape::Shape() (lxx) — the resulting shape; None
    /// carries the OCCT null shape.
    pub fn shape(&self) -> Option<&Shape> {
        self.my_res.as_ref()
    }
}
