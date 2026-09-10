//! OCCT BRepFill_Generator::Perform and accessors — split from
//! generator.rs (file-size rule): the inherent impl of the class continues
//! here; the struct itself lives in generator.rs.

use glam::{DVec2, DVec3};
use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{BezierCurve3, Curve2d, Curve3, Line2d, Surface3, TrimmedCurve3};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

use super::generator::{
    create_k_part, create_new_edge, curve_first_last, curve_reversed, detect_k_part, is_edge_closed,
    is_same_parameter, is_same_range, make_edge, surface_bounds, surface_uiiso, shape_key, shape_reversed,
    shape_oriented, top_exp_vertices, top_exp_wire_vertices, bezier2, bind_pcurve, bind_seam_pcurves,
    BRepFillGenerator, BRepFillThruSectionErrorStatus, GeomFillGenerator, ShapeKey,
};
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;

/// OCCT Precision::Confusion().
const TOL_CONFUSION: f64 = CONFUSION;
/// OCCT Precision::PConfusion().
const TOL_PCONFUSION: f64 = rcad_kernel::core::precision::PCONFUSION;

impl BRepFillGenerator {
    /// OCCT BRepFill_Generator::BRepFill_Generator (L598-602).
    pub fn new() -> Self {
        BRepFillGenerator {
            my_wires: Vec::new(),
            my_shell: None,
            my_map: HashMap::new(),
            my_old_new_shapes: HashMap::new(),
            my_reshaper: ShapeBuildReShape::new(),
            my_mutable_input: true,
            my_status: BRepFillThruSectionErrorStatus::NotDone,
        }
    }

    /// OCCT BRepFill_Generator::AddWire (L606-609).
    pub fn add_wire(&mut self, wire: Shape) {
        self.my_wires.push(wire);
    }

    /// OCCT BRepFill_Generator::Shell (hxx L50).
    pub fn shell(&self) -> Option<Shape> {
        self.my_shell.clone()
    }

    /// OCCT BRepFill_Generator::GeneratedShapes (L1190-1202).
    pub fn generated_shapes(&self, ssection: &Shape) -> Vec<Shape> {
        self.my_map
            .get(&shape_key(ssection))
            .cloned()
            .unwrap_or_default()
    }

    /// OCCT BRepFill_Generator::Generated (L1206-1210).
    pub fn generated(&self) -> &HashMap<ShapeKey, Vec<Shape>> {
        &self.my_map
    }

    /// OCCT BRepFill_Generator::ResultShape (L1214-1225).
    pub fn result_shape(&mut self, brep: &mut BRep, the_shape: &Shape) -> Shape {
        let mut a_new_shape = match self.my_old_new_shapes.get(&shape_key(the_shape)) {
            Some(s) => s.clone(),
            None => the_shape.clone(),
        };
        loop {
            let a_prev_shape = a_new_shape.clone();
            a_new_shape = self.my_reshaper.value(brep, &a_new_shape);
            if a_new_shape.ptr_id() == a_prev_shape.ptr_id() {
                break;
            }
        }
        a_new_shape
    }

    /// OCCT BRepFill_Generator::SetMutableInput (L1229-1232).
    pub fn set_mutable_input(&mut self, the_is_mutable_input: bool) {
        self.my_mutable_input = the_is_mutable_input;
    }

    /// OCCT BRepFill_Generator::IsMutableInput (L1236-1239).
    pub fn is_mutable_input(&self) -> bool {
        self.my_mutable_input
    }

    /// OCCT BRepFill_Generator::GetStatus (hxx L76).
    pub fn get_status(&self) -> BRepFillThruSectionErrorStatus {
        self.my_status
    }

    /// OCCT BRepFill_Generator::Perform (L613-1186).
    pub fn perform(&mut self, brep: &mut BRep) {
        self.my_status = BRepFillThruSectionErrorStatus::Done;

        let nb = self.my_wires.len();
        // NCollection_IndexedMap aModifWires — insertion-ordered set.
        let mut a_modif_wires: Vec<Shape> = Vec::new();

        let mut a_first_wire = true;
        let mut w_point1;
        let mut w_point2;
        let mut u_closed = false;
        let mut degen_first = false;
        let mut degen_last = false;

        // Faces accumulated for myShell (B.Add(myShell, Face) per face; the
        // shell TShape is built once at the end — architecture note).
        let mut faces: Vec<Shape> = Vec::new();

        // OCCT BRepTools_WireExplorer — TWireData::edges is the wire-ordered
        // edge list, each entry carrying its cumulated orientation.
        let wire_edges = |b: &BRep, w: &Shape| -> Vec<Shape> {
            match w.data.as_ref() {
                TShape::Wire(wd) => wd.edges.clone(),
                _ => Vec::new(),
            }
        };

        for i in 0..nb.saturating_sub(1) {
            let wire1 = self.my_wires[i].clone();
            let wire2 = self.my_wires[i + 1].clone();

            w_point1 = false;
            if i == 0 {
                w_point1 = true;
                for e in wire_edges(brep, &wire1) {
                    w_point1 = w_point1 && brep.edge(e).degenerated;
                }
                degen_first = w_point1;

                let (v1, v2) = top_exp_wire_vertices(brep, &wire1);
                u_closed = v1.is_same(&v2);
            }

            w_point2 = false;
            if i == nb - 2 {
                w_point2 = true;
                for e in wire_edges(brep, &wire2) {
                    w_point2 = w_point2 && brep.edge(e).degenerated;
                }
                degen_last = w_point2;
            }

            let ex1 = wire_edges(brep, &wire1);
            let ex2 = wire_edges(brep, &wire2);
            let mut ix1 = 0usize;
            let mut ix2 = 0usize;

            let mut map: HashMap<ShapeKey, Shape> = HashMap::new();

            let mut tantque = ix1 < ex1.len() && ix2 < ex2.len();

            while tantque {
                let an_or_edge1 = ex1[ix1].clone();
                let an_or_edge2 = ex2[ix2].clone();

                let degen1 = brep.edge(an_or_edge1.clone()).degenerated;
                let degen2 = brep.edge(an_or_edge2.clone()).degenerated;

                let mut edge1: Shape;
                let mut edge2: Shape;

                if degen1 {
                    // TopoDS_Shape aLocalShape = anOrEdge1.EmptyCopied();
                    edge1 = brep.empty_copied(&an_or_edge1);
                    let (v1f, v1l) = top_exp_vertices(brep, &an_or_edge1);
                    // B.Add(Edge1, V1f FORWARD); B.Add(Edge1, V1l REVERSED);
                    // B.Range(Edge1, 0, 1)
                    let e1 = brep.edge_mut_inplace(edge1.clone());
                    e1.first = shape_oriented(&v1f, Orientation::Forward);
                    e1.last = shape_oriented(&v1l, Orientation::Reversed);
                    e1.my_shapes = vec![e1.first.clone(), e1.last.clone()];
                    e1.range = [0.0, 1.0];
                    self.my_old_new_shapes
                        .insert(shape_key(&an_or_edge1), edge1.clone());
                } else {
                    edge1 = an_or_edge1.clone();
                }

                if degen2 {
                    edge2 = brep.empty_copied(&an_or_edge2);
                    let (v2f, v2l) = top_exp_vertices(brep, &an_or_edge2);
                    let e2 = brep.edge_mut_inplace(edge2.clone());
                    e2.first = shape_oriented(&v2f, Orientation::Forward);
                    e2.last = shape_oriented(&v2l, Orientation::Reversed);
                    e2.my_shapes = vec![e2.first.clone(), e2.last.clone()];
                    e2.range = [0.0, 1.0];
                    self.my_old_new_shapes
                        .insert(shape_key(&an_or_edge2), edge2.clone());
                } else {
                    edge2 = an_or_edge2.clone();
                }

                // bool Periodic = (IsClosed(Edge1) || degen1) && (IsClosed(Edge2) || degen2);
                let periodic = (is_edge_closed(brep, &edge1) || degen1)
                    && (is_edge_closed(brep, &edge2) || degen2);

                // ATTENTION : a non-punctual wire should not
                //             contain a punctual edge
                if !w_point1 {
                    ix1 += 1;
                }
                if !w_point2 {
                    ix2 += 1;
                }

                // initialization of vertices
                let mut f1 = 0.0f64;
                let mut l1 = 1.0f64;
                let mut f2 = 0.0f64;
                let mut l2 = 1.0f64;
                let (v1f, v1l) = if edge1.orientation == Orientation::Reversed {
                    // TopExp::Vertices(Edge1, V1l, V1f)
                    let ed = brep.edge(edge1.clone());
                    (ed.last.clone(), ed.first.clone())
                } else {
                    let ed = brep.edge(edge1.clone());
                    (ed.first.clone(), ed.last.clone())
                };
                let (v2f, v2l) = if edge2.orientation == Orientation::Reversed {
                    let ed = brep.edge(edge2.clone());
                    (ed.last.clone(), ed.first.clone())
                } else {
                    let ed = brep.edge(edge2.clone());
                    (ed.first.clone(), ed.last.clone())
                };
                let vf_tomap = if degen1 { v2f.clone() } else { v1f.clone() };
                let vl_tomap = if degen1 { v2l.clone() } else { v1l.clone() };

                // processing of KPart
                let itype = detect_k_part(brep, &edge1, &edge2);
                if itype == -1 {
                    self.my_status = BRepFillThruSectionErrorStatus::Null3DCurve;
                    return;
                }

                let mut surf: Option<Surface3> = None;
                if itype == 0 {
                    // no part cases
                    let mut c1: Curve3;
                    if degen1 {
                        let e1 = brep.vertex(v1f.clone()).point;
                        let e2 = brep.vertex(v1l.clone()).point;
                        c1 = Curve3::Bezier(BezierCurve3 {
                            control_points: vec![e1, e2],
                            weights: vec![1.0, 1.0],
                        });
                        f1 = 0.0;
                        l1 = 1.0;
                    } else {
                        let ed = brep.edge(edge1.clone());
                        match &ed.curve {
                            Some(c) => {
                                c1 = c.clone();
                                f1 = ed.range[0];
                                l1 = ed.range[1];
                            }
                            None => {
                                self.my_status = BRepFillThruSectionErrorStatus::Null3DCurve;
                                return;
                            }
                        }
                    }
                    let mut c2: Curve3;
                    if degen2 {
                        let e1 = brep.vertex(v2l.clone()).point;
                        let e2 = brep.vertex(v2f.clone()).point;
                        c2 = Curve3::Bezier(BezierCurve3 {
                            control_points: vec![e1, e2],
                            weights: vec![1.0, 1.0],
                        });
                        f2 = 0.0;
                        l2 = 1.0;
                    } else {
                        let ed = brep.edge(edge2.clone());
                        match &ed.curve {
                            Some(c) => {
                                c2 = c.clone();
                                f2 = ed.range[0];
                                l2 = ed.range[1];
                            }
                            None => {
                                self.my_status = BRepFillThruSectionErrorStatus::Null3DCurve;
                                return;
                            }
                        }
                    }

                    // compute the location
                    let same_loc = false;

                    // transform and trim the curves
                    {
                        let (fp, lp) = curve_first_last(&c1);
                        if (f1 - fp).abs() > TOL_PCONFUSION || (l1 - lp).abs() > TOL_PCONFUSION {
                            c1 = Curve3::Trimmed(TrimmedCurve3::new(c1, f1, l1));
                        }
                        // else: C1->Copy() — rcad Curve3 values are copies.
                    }
                    if !same_loc && edge1.orientation == Orientation::Reversed {
                        // OCCT: C1->Transform(L1.Transformation()) (identity)
                        // then C1->Reverse().
                        c1 = curve_reversed(&c1);
                    }

                    {
                        let (fp, lp) = curve_first_last(&c2);
                        if (f2 - fp).abs() > TOL_PCONFUSION || (l2 - lp).abs() > TOL_PCONFUSION {
                            c2 = Curve3::Trimmed(TrimmedCurve3::new(c2, f2, l2));
                        }
                    }
                    if !same_loc && edge2.orientation == Orientation::Reversed {
                        c2 = curve_reversed(&c2);
                    }

                    let mut generator = GeomFillGenerator::new();
                    generator.add_curve(&c1);
                    generator.add_curve(&c2);
                    generator.perform(TOL_PCONFUSION);

                    match generator.surface() {
                        Some(s) => surf = Some(Surface3::BSpline(s)),
                        None => {
                            self.my_status = BRepFillThruSectionErrorStatus::Null3DCurve;
                            return;
                        }
                    }
                } else {
                    // particular case
                    match create_k_part(brep, &edge1, &edge2, itype) {
                        Some(s) => surf = Some(s),
                        None => {
                            self.my_status = BRepFillThruSectionErrorStatus::Null3DCurve;
                            return;
                        }
                    }
                }

                // make the missing edges
                let surf_ref = surf.clone().unwrap();
                let (bf1, bl1, bf2, bl2) = surface_bounds(&surf_ref);
                f1 = bf1;
                l1 = bl1;
                f2 = bf2;
                l2 = bl2;
                let (first, last) = if itype == 0 { (f2, l2) } else { (0.0, 1.0) };

                let mut edge3: Shape;
                if let Some(prev) = map.get(&shape_key(&vf_tomap)) {
                    // Edge3 = Map(Vf_toMap).Reversed()
                    edge3 = shape_reversed(prev);
                } else {
                    if v1f.is_same(&v2f) {
                        // B.MakeEdge(Edge3); B.Degenerated(Edge3, true);
                        edge3 = make_edge(
                            brep,
                            None,
                            shape_oriented(&v1f, Orientation::Forward),
                            shape_oriented(&v2f, Orientation::Reversed),
                            [first, last],
                        );
                        brep.edge_mut_inplace(edge3.clone()).degenerated = true;
                    } else {
                        let cc: Curve3 = if itype == 0 {
                            // general case : Edge3 corresponds to iso U=f1
                            surface_uiiso(&surf_ref, f1)
                        } else {
                            // particular case : it is required to calculate the curve 3d
                            let p1 = brep.vertex(v1f.clone()).point;
                            let p2 = brep.vertex(v2f.clone()).point;
                            Curve3::Bezier(BezierCurve3 {
                                control_points: vec![p1, p2],
                                weights: vec![1.0, 1.0],
                            })
                        };
                        edge3 = make_edge(
                            brep,
                            Some(cc),
                            shape_oriented(&v1f, Orientation::Forward),
                            shape_oriented(&v2f, Orientation::Reversed),
                            [first, last],
                        );
                    }
                    // Edge3.Reverse(); Map.Bind(Vf_toMap, Edge3);
                    map.insert(shape_key(&vf_tomap), shape_reversed(&edge3));
                }

                let mut common_edge = false;
                if let Some(prev) = map.get(&shape_key(&vl_tomap)) {
                    let common_e = shape_reversed(prev);
                    let (vv1, vv2) = top_exp_vertices(brep, &common_e);
                    common_edge = vv1.is_same(&v1l) && vv2.is_same(&v2l);
                }
                let mut edge4: Shape;
                if common_edge {
                    // Edge4 = Map(Vl_toMap).Reversed()
                    let prev = map.get(&shape_key(&vl_tomap)).unwrap();
                    edge4 = shape_reversed(prev);
                } else {
                    if v1l.is_same(&v2l) {
                        edge4 = make_edge(
                            brep,
                            None,
                            shape_oriented(&v1l, Orientation::Forward),
                            shape_oriented(&v2l, Orientation::Reversed),
                            [first, last],
                        );
                        brep.edge_mut_inplace(edge4.clone()).degenerated = true;
                    } else {
                        let cc: Curve3 = if itype == 0 {
                            // general case : Edge4 corresponds to iso U=l1
                            surface_uiiso(&surf_ref, l1)
                        } else {
                            let p1 = brep.vertex(v1l.clone()).point;
                            let p2 = brep.vertex(v2l.clone()).point;
                            Curve3::Bezier(BezierCurve3 {
                                control_points: vec![p1, p2],
                                weights: vec![1.0, 1.0],
                            })
                        };
                        edge4 = make_edge(
                            brep,
                            Some(cc),
                            shape_oriented(&v1l, Orientation::Forward),
                            shape_oriented(&v2l, Orientation::Reversed),
                            [first, last],
                        );
                    }
                    // Map.Bind(Vl_toMap, Edge4);
                    map.insert(shape_key(&vl_tomap), edge4.clone());
                }

                if !self.my_mutable_input {
                    if !degen1 {
                        // if true=>already empty-copied
                        if let Some(new_ed) = self.my_old_new_shapes.get(&shape_key(&edge1)) {
                            edge1 = new_ed.clone();
                        } else if a_first_wire
                            && (itype != 4
                                || is_same_parameter(brep, &edge1)
                                || is_same_range(brep, &edge1))
                        {
                            // if such expression is true and mutableInput is false => pre-copy the edge to prevent
                            // a following modifying (see code below)
                            edge1 = create_new_edge(
                                brep,
                                &edge1,
                                &mut self.my_old_new_shapes,
                                &wire1,
                                &mut a_modif_wires,
                            );
                        }
                    }
                    if !degen2 {
                        if let Some(new_ed) = self.my_old_new_shapes.get(&shape_key(&edge2)) {
                            edge2 = new_ed.clone();
                        } else if itype != 4
                            || is_same_parameter(brep, &edge2)
                            || is_same_range(brep, &edge2)
                        {
                            edge2 = create_new_edge(
                                brep,
                                &edge2,
                                &mut self.my_old_new_shapes,
                                &wire2,
                                &mut a_modif_wires,
                            );
                        }
                    }
                }

                // set the pcurves
                let _t = TOL_CONFUSION;

                // make the wire (OCCT B.MakeWire(aWire) + the four B.Add)
                let mut a_wire_edges: Vec<Shape> = Vec::new();
                if !degen1 || itype != 4 {
                    a_wire_edges.push(edge1.clone());
                }
                a_wire_edges.push(edge4.clone());
                if !degen2 || itype != 4 {
                    a_wire_edges.push(shape_reversed(&edge2));
                }
                a_wire_edges.push(edge3.clone());
                let a_wire = brep.add_twire(a_wire_edges);

                // OCCT B.MakeFace(Face, Surf, Precision::Confusion()) +
                // B.Add(Face, aWire) — architecture note: rcad builds the
                // face from the wire in one step (see the file header).
                let face = brep.add_tface_tol(
                    Some(surf_ref.clone()),
                    a_wire,
                    vec![],
                    None,
                    None,
                    vec![],
                    false,
                    TOL_CONFUSION,
                );

                let fkey = (face.ptr_id(), face.location);
                if itype != 4 {
                    // not plane
                    if edge1.orientation == Orientation::Reversed {
                        let pc = Curve2d::Line(Line2d::new(
                            DVec2::new(0.0, f2),
                            DVec2::new(-1.0, 0.0),
                        ));
                        bind_pcurve(brep, &edge1, fkey, pc, -l1, -f1);
                    } else {
                        let pc = Curve2d::Line(Line2d::new(
                            DVec2::new(0.0, f2),
                            DVec2::new(1.0, 0.0),
                        ));
                        bind_pcurve(brep, &edge1, fkey, pc, f1, l1);
                    }

                    if edge2.orientation == Orientation::Reversed {
                        let pc = Curve2d::Line(Line2d::new(
                            DVec2::new(0.0, l2),
                            DVec2::new(-1.0, 0.0),
                        ));
                        bind_pcurve(brep, &edge2, fkey, pc, -l1, -f1);
                    } else {
                        let pc = Curve2d::Line(Line2d::new(
                            DVec2::new(0.0, l2),
                            DVec2::new(1.0, 0.0),
                        ));
                        bind_pcurve(brep, &edge2, fkey, pc, f1, l1);
                    }
                }

                if itype == 0 {
                    if periodic {
                        let pc1 = Curve2d::Line(Line2d::new(
                            DVec2::new(l1, 0.0),
                            DVec2::new(0.0, 1.0),
                        ));
                        let pc2 = Curve2d::Line(Line2d::new(
                            DVec2::new(f1, 0.0),
                            DVec2::new(0.0, 1.0),
                        ));
                        bind_seam_pcurves(brep, &edge3, fkey, pc1, pc2, 0.0, 1.0);
                    } else {
                        let pc3 = Curve2d::Line(Line2d::new(
                            DVec2::new(f1, 0.0),
                            DVec2::new(0.0, 1.0),
                        ));
                        let pc4 = Curve2d::Line(Line2d::new(
                            DVec2::new(l1, 0.0),
                            DVec2::new(0.0, 1.0),
                        ));
                        bind_pcurve(brep, &edge3, fkey, pc3, 0.0, 1.0);
                        bind_pcurve(brep, &edge4, fkey, pc4, 0.0, 1.0);
                    }
                } else {
                    // KPart
                    if periodic {
                        let e1 = bezier2([l1, f2], [l1, l2]);
                        let e2 = bezier2([f1, f2], [f1, l2]);
                        bind_seam_pcurves(brep, &edge3, fkey, e1, e2, 0.0, 1.0);
                    } else if itype != 4 {
                        // not plane
                        let e2 = bezier2([f1, f2], [f1, l2]);
                        let e1 = bezier2([l1, f2], [l1, l2]);
                        bind_pcurve(brep, &edge3, fkey, e2, 0.0, 1.0);
                        bind_pcurve(brep, &edge4, fkey, e1, 0.0, 1.0);
                    }
                }

                // Set the non parameter flag;
                for e in [&edge1, &edge2, &edge3, &edge4] {
                    let ed = brep.edge_mut_inplace(e.clone());
                    ed.same_parameter = false;
                    ed.same_range = false;
                }

                // B.Add(Face, aWire); B.Add(myShell, Face) — the wire is
                // carried by the face and the face is collected for the shell.
                faces.push(face.clone());

                // complete myMap for edge1
                if !degen1 || itype != 4 {
                    let a_red = if degen1 { edge1.clone() } else { an_or_edge1.clone() };
                    self.my_map
                        .entry(shape_key(&a_red))
                        .or_default()
                        .push(face.clone());
                }

                tantque = ix1 < ex1.len() && ix2 < ex2.len();
                if w_point1 {
                    tantque = ix2 < ex2.len();
                }
                if w_point2 {
                    tantque = ix1 < ex1.len();
                }
            }
            a_first_wire = false;
        }

        // B.MakeShell(myShell)
        let shell = brep.add_tshell(std::mem::take(&mut faces));
        self.my_shell = Some(shell.clone());

        // all vertices from myShell are the part of orig. section wires
        if self.my_mutable_input {
            // BRepLib::SameParameter(myShell)
            brep.same_parameter();
        } else {
            // for each (K, V) in myOldNewShapes: myReshaper.Replace(K, V)
            let pairs: Vec<(Shape, Shape)> = self
                .my_old_new_shapes
                .iter()
                .filter_map(|(k, v)| find_by_key(brep, *k).map(|s| (s, v.clone())))
                .collect();
            for (a_k, a_val) in pairs {
                self.my_reshaper.replace(brep, &a_k, &a_val);
            }
            // BRepLib::SameParameter(myShell, myReshaper) — GAP: the
            // reshaper-aware overload is not ported; the plain walk is used.
            brep.same_parameter();
            let new_shell = self
                .my_reshaper
                .apply(brep, &shell, ShapeType::Shell);
            self.my_shell = Some(new_shell);
        }

        if u_closed && degen_first && degen_last {
            // myShell.Closed(true)
            if let Some(sh) = &self.my_shell {
                let idx = sh.index;
                if let TShape::Shell(sd) = Arc::make_mut(&mut brep.tshapes[idx]) {
                    sd.flags |= tshape_flags::CLOSED;
                }
            }
        }

        // update wire's history
        let modif = std::mem::take(&mut a_modif_wires);
        for a_cur_wire in &modif {
            let wd_edges = match a_cur_wire.data.as_ref() {
                TShape::Wire(wd) => wd.edges.clone(),
                _ => Vec::new(),
            };
            let mut new_edges: Vec<Shape> = Vec::with_capacity(wd_edges.len());
            for a_cur_edge in &wd_edges {
                // edges only
                let a_new_edge = self.result_shape(brep, a_cur_edge);
                new_edges.push(a_new_edge);
            }
            let a_new_wire = brep.add_twire(new_edges);
            // aNewWire.Free(aCurWire.Free()); ... .Convex(aCurWire.Convex())
            if let TShape::Wire(cwd) = a_cur_wire.data.as_ref() {
                let idx = a_new_wire.index;
                if let TShape::Wire(nwd) = Arc::make_mut(&mut brep.tshapes[idx]) {
                    nwd.flags = cwd.flags;
                }
            }
            self.my_old_new_shapes.insert(shape_key(a_cur_wire), a_new_wire);
        }
    }
}

fn find_by_key(brep: &BRep, k: ShapeKey) -> Option<Shape> {
    for (i, ts) in brep.tshapes.iter().enumerate() {
        if Arc::as_ptr(ts) as u64 == k.0 {
            return Some(brep.shape_at(i));
        }
    }
    None
}
