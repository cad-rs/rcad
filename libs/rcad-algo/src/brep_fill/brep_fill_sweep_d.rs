//! OCCT BRepFill_Sweep.cxx, part D (members, cxx L3213-3530 + L3587-3906 +
//! L4047-4151) — Build, PerformCorner and RebuildTopOrBottomEdge, the
//! companion of [`super::brep_fill_sweep`] (part A), [`super::
//! brep_fill_sweep_b`] (part-B statics) and [`super::brep_fill_sweep_c`]
//! (CorrectApproxParameters / BuildShell).
//!
//! Architecture differences:
//! - OCCT `BRepFill_Right` maps to the rcad `RightCorner` variant (the
//!   part-A enum-value mapping of BRepFill_TransitionStyle).
//! - `TopExp::MapShapesAndAncestors(shell, EDGE, FACE, EFmap)` maps to the
//!   local face/edge explorer walk (same ancestor lists, insertion order).
//! - `TopoDS_Shape::Free(bool)` maps to the FREE TShape flag.

use std::collections::HashSet;

use glam::DVec3;

use rcad_kernel::core::precision::{p_confusion, CONFUSION};
use rcad_kernel::geom::{CurveEval, Plane, Surface3};
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::topo::topods::{tshape_flags, BRep, BRepBuilder, GeomAbsShape, Orientation, Shape, TShape};

use crate::brep_fill::brep_fill_pipe_shell_b::{BRepFillTransitionStyle, ShapeHArray2, ShapeToArray2Map};
use crate::brep_fill::brep_fill_sweep::{translate, BRepFillSweep};
use crate::brep_fill::generator::ShapeKey;

use super::brep_fill_sweep_b::GpAx2;
use super::brep_fill_sweep_b::gp_vec_angle;

/// The surface value of a pool face resolved by its identity (the OCCT
/// BRep_GCurve::Surface read through the rcad representation key).
fn brep_face_surface(brep: &BRep, face_ptr: u64) -> Option<Surface3> {
    let idx = brep.index_by_ptr(face_ptr)?;
    match brep.tshapes[idx].as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT TopExp::MapShapesAndAncestors(shell, EDGE, FACE, EFmap) — the
/// (edge, ancestor faces) pairs in the shell traversal order.
fn map_shapes_and_ancestors_edge_face(shell: &Shape) -> Vec<(Shape, Vec<Shape>)> {
    let mut efmap: Vec<(Shape, Vec<Shape>)> = Vec::new();
    let faces = crate::brep_algo::tool::explorer(
        shell,
        rcad_kernel::topo::topods::ShapeType::Face,
        rcad_kernel::topo::topods::ShapeType::Compound,
    );
    for face in faces {
        let edges = crate::brep_algo::tool::explorer(
            &face,
            rcad_kernel::topo::topods::ShapeType::Edge,
            rcad_kernel::topo::topods::ShapeType::Compound,
        );
        for edge in edges {
            if let Some(entry) = efmap.iter_mut().find(|(e, _)| e.is_same(&edge)) {
                entry.1.push(face.clone());
            } else {
                efmap.push((edge, vec![face.clone()]));
            }
        }
    }
    efmap
}

/// OCCT GeomAbs_Shape — the rcad dual-enum bridge (the kernel carries the
/// approximation-side `math::GeomAbsShape` C0..CN set and the topology-side
/// `topods::GeomAbsShape` C0..CN + G1/G2 set; the Build Continuity feeds the
/// approximation stack, so the G-variants fold onto their C levels).
fn to_math_geom_abs(c: GeomAbsShape) -> rcad_kernel::math::GeomAbsShape {
    match c {
        GeomAbsShape::C0 => rcad_kernel::math::GeomAbsShape::C0,
        GeomAbsShape::G1 | GeomAbsShape::C1 => rcad_kernel::math::GeomAbsShape::C1,
        GeomAbsShape::G2 | GeomAbsShape::C2 => rcad_kernel::math::GeomAbsShape::C2,
        GeomAbsShape::C3 => rcad_kernel::math::GeomAbsShape::C3,
        GeomAbsShape::CN => rcad_kernel::math::GeomAbsShape::CN,
    }
}

impl BRepFillSweep {
    /// OCCT Build (cxx L3213-3530) — construct the result of sweeping.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        &mut self,
        brep: &mut BRep,
        reversed_edges: &mut HashSet<ShapeKey>,
        tapes: &mut ShapeToArray2Map,
        rails: &mut ShapeToArray2Map,
        transition: BRepFillTransitionStyle,
        continuity: GeomAbsShape,
        approx: crate::geomalgo::geomfill::sweep::GeomFillApproxStyle,
        degmax: i32,
        segmax: i32,
    ) {
        // myContinuity = Continuity; myApproxStyle = Approx;
        // myDegmax = Degmax; mySegmax = Segmax;
        self.my_continuity = to_math_geom_abs(continuity);
        self.my_approx_style = approx;
        self.my_degmax = degmax;
        self.my_segmax = segmax;

        // CorrectApproxParameters();
        self.correct_approx_parameters(brep);

        // Wire
        if self.my_sec.borrow().is_vertex() {
            // isDone = BuildWire(Transition);
            self.is_done = self.build_wire(brep, transition);
        } else {
            // Shell
            let nb_trous = self.my_loc.borrow_mut().base_mut().nb_holes(brep, self.my_tol3d);
            let nb_path = self.my_loc.borrow().base().nb_law();
            let nb_law = self.my_sec.borrow().base().nb_law();
            let mut nb_part = 1i32;
            let mut b = BRepBuilder::new();
            // myUEdges = new (1, NbLaw + 1, 1, NbPath);
            self.my_u_edges =
                Some(vec![vec![Shape::null(); nb_path as usize]; (nb_law + 1) as usize]);
            // myVEdges = new (1, NbLaw, 1, NbPath + 1);
            self.my_v_edges =
                Some(vec![vec![Shape::null(); (nb_path + 1) as usize]; nb_law as usize]);
            // myFaces = new (1, NbLaw, 1, NbPath);
            self.my_faces = Some(vec![vec![Shape::null(); nb_path as usize]; nb_law as usize]);
            // myTapes = new (1, NbLaw);
            self.my_tapes = Some(vec![Shape::null(); nb_law as usize]);
            let mut bb = BRepBuilder::new();
            for i in 1..=nb_law {
                // TopoDS_Shell aShell; BB.MakeShell(aShell);
                let a_shell = bb.make_shell(brep);
                // myTapes->ChangeValue(i) = aShell;
                self.my_tapes.as_mut().expect("myTapes")[i as usize - 1] = a_shell;
            }
            // occ::handle<NCollection_HArray2<TopoDS_Shape>> Bounds =
            //   new (1, NbLaw, 1, 2);
            let mut bounds: ShapeHArray2 = vec![vec![Shape::null(); 2]; nb_law as usize];

            // occ::handle<NCollection_HArray1<int>> Trous;
            let mut trous: Vec<i32> = Vec::new();

            if nb_trous > 0 {
                // How many sub-parts ?
                trous = vec![0i32; nb_trous as usize];
                // myLoc->Holes(Trous->ChangeArray1());
                self.my_loc.borrow().base().holes(&mut trous);
                nb_part += nb_trous;
                if trous[nb_trous as usize - 1] == nb_path + 1 {
                    nb_part -= 1;
                }
            }
            if nb_part == 1 {
                // This is done at once
                let mut extend = 0.0f64;
                if nb_trous == 1 {
                    // Extend = EvalExtrapol(1, Transition);
                    extend = self.eval_extrapol(brep, 1, transition);
                }
                // isDone = BuildShell(Transition, 1, NbPath + 1, ...);
                self.is_done = self.build_shell(
                    brep,
                    transition,
                    1,
                    nb_path + 1,
                    reversed_edges,
                    tapes,
                    rails,
                    extend,
                    extend,
                );
            } else {
                // This is done piece by piece
                let mut ifirst = 1i32;
                let mut ilast;
                let mut ii = 1i32;
                self.is_done = true;
                while ii <= nb_part && self.is_done {
                    if ii > nb_trous {
                        ilast = nb_path + 1;
                    } else {
                        ilast = trous[ii as usize - 1];
                    }
                    // BuildShell(Transition, IFirst, ILast, ...,
                    //            EvalExtrapol(IFirst, Transition),
                    //            EvalExtrapol(ILast, Transition));
                    let extrap_first = self.eval_extrapol(brep, ifirst, transition);
                    let extrap_last = self.eval_extrapol(brep, ilast, transition);
                    self.is_done = self.build_shell(
                        brep,
                        transition,
                        ifirst,
                        ilast,
                        reversed_edges,
                        tapes,
                        rails,
                        extrap_first,
                        extrap_last,
                    );
                    if ifirst > 1 {
                        // Translate(myVEdges, IFirst, Bounds, 2);
                        let my_v_edges = self.my_v_edges.as_ref().expect("myVEdges").clone();
                        translate(&my_v_edges, ifirst, &mut bounds, 2);
                        // if (!PerformCorner(IFirst, Transition, Bounds))
                        if !self.perform_corner(brep, ifirst, transition, &bounds) {
                            self.is_done = false;
                            return;
                        }
                    }
                    ifirst = ilast;
                    // Translate(myVEdges, IFirst, Bounds, 1);
                    let my_v_edges = self.my_v_edges.as_ref().expect("myVEdges").clone();
                    translate(&my_v_edges, ifirst, &mut bounds, 1);
                    ii += 1;
                }
            }
            // Management of looping ends
            if (nb_trous > 0)
                && self.my_loc.borrow().base().is_closed(brep)
                && trous[nb_trous as usize - 1] == nb_path + 1
            {
                // Translate(myVEdges, NbPath + 1, Bounds, 1);
                let my_v_edges = self.my_v_edges.as_ref().expect("myVEdges").clone();
                translate(&my_v_edges, nb_path + 1, &mut bounds, 1);
                // Translate(myVEdges, 1, Bounds, 2);
                let my_v_edges = self.my_v_edges.as_ref().expect("myVEdges").clone();
                translate(&my_v_edges, 1, &mut bounds, 2);
                // if (!PerformCorner(1, Transition, Bounds))
                if !self.perform_corner(brep, 1, transition, &bounds) {
                    self.is_done = false;
                    return;
                }
                // Translate(myVEdges, 1, myVEdges, NbPath + 1);
                let my_v_edges = self.my_v_edges.as_ref().expect("myVEdges").clone();
                translate(
                    &my_v_edges,
                    1,
                    self.my_v_edges.as_mut().expect("myVEdges"),
                    nb_path + 1,
                );
            }

            // Construction of the shell
            // TopoDS_Shell shell; B.MakeShell(shell);
            let mut shell = b.make_shell(brep);
            let mut a_nb_faces = 0i32;
            for ipath in 1..=nb_path {
                for isec in 1..=nb_law {
                    let face = self.my_faces.as_ref().expect("myFaces")[isec as usize - 1]
                        [(ipath - 1) as usize]
                        .clone();
                    if !face.is_null() && matches!(face.data.as_ref(), TShape::Face(_)) {
                        // B.Add(shell, face);
                        b.add_to_shell(brep, shell.clone(), face);
                        a_nb_faces += 1;
                    }
                }
            }

            // if ((NbTrous > 0) ? (aNbFaces < NbLaw) : (aNbFaces == 0))
            let faces_check = if nb_trous > 0 {
                a_nb_faces < nb_law
            } else {
                a_nb_faces == 0
            };
            if faces_check {
                self.is_done = false;
                return;
            }

            // NCollection_List<TopoDS_Shape>::Iterator It(myAuxShape);
            for face in self.my_aux_shape.clone() {
                if !face.is_null() && matches!(face.data.as_ref(), TShape::Face(_)) {
                    // B.Add(shell, face);
                    b.add_to_shell(brep, shell.clone(), face);
                }
            }
            // Set common Uedges to faces
            // BRepTools_Substitution aSubstitute;
            let mut a_substitute =
                crate::topalgo::brep_tools_substitution::BRepToolsSubstitution::new();
            // NCollection_DataMap<...>::Iterator mapit(myVEdgesModified);
            // (the OCCT map iteration order is unspecified; the rcad HashMap
            // order stands in)
            let v_edges_modified = self.my_v_edges_modified.clone();
            for (old_key, new_edge) in v_edges_modified.iter() {
                // const TopoDS_Edge& OldEdge = TopoDS::Edge(mapit.Key());
                let old_edge = match brep.index_by_ptr(old_key.0) {
                    Some(idx) => brep.shape_at(idx),
                    None => Shape::null(),
                };
                if old_edge.is_null() {
                    continue;
                }
                // Substitute(aSubstitute, OldEdge, NewEdge);
                super::brep_fill_sweep_b::substitute(brep, &mut a_substitute, &old_edge, new_edge);
            }
            // aSubstitute.Build(shell);
            a_substitute.build(brep, &shell);
            if a_substitute.is_copied(&shell) {
                // const NCollection_List<TopoDS_Shape>& listSh =
                //   aSubstitute.Copy(shell);
                // shell = TopoDS::Shell(listSh.First());
                let list_sh = a_substitute.copy(&shell);
                shell = list_sh.first().expect("listSh.First()").clone();
            }

            // myFaces replacement loop
            for ii in 0..self.my_faces.as_ref().expect("myFaces").len() {
                for jj in 0..self.my_faces.as_ref().expect("myFaces")[ii].len() {
                    let a_local_shape = self.my_faces.as_ref().expect("myFaces")[ii][jj].clone();
                    if !a_local_shape.is_null() && a_substitute.is_copied(&a_local_shape) {
                        let a_list = a_substitute.copy(&a_local_shape);
                        if !a_list.is_empty() {
                            // myFaces->ChangeValue(ii, jj) = aList.First();
                            self.my_faces.as_mut().expect("myFaces")[ii][jj] = a_list[0].clone();
                        }
                    }
                }
            }

            // myVEdges replacement loop
            for ii in 0..self.my_v_edges.as_ref().expect("myVEdges").len() {
                for jj in 0..self.my_v_edges.as_ref().expect("myVEdges")[ii].len() {
                    let a_local_shape = self.my_v_edges.as_ref().expect("myVEdges")[ii][jj].clone();
                    if !a_local_shape.is_null() && a_substitute.is_copied(&a_local_shape) {
                        let a_list = a_substitute.copy(&a_local_shape);
                        if !a_list.is_empty() {
                            self.my_v_edges.as_mut().expect("myVEdges")[ii][jj] = a_list[0].clone();
                        }
                    }
                }
            }

            // myUEdges replacement loop
            for ii in 0..self.my_u_edges.as_ref().expect("myUEdges").len() {
                for jj in 0..self.my_u_edges.as_ref().expect("myUEdges")[ii].len() {
                    let a_local_shape = self.my_u_edges.as_ref().expect("myUEdges")[ii][jj].clone();
                    if !a_local_shape.is_null() && a_substitute.is_copied(&a_local_shape) {
                        let a_list = a_substitute.copy(&a_local_shape);
                        if !a_list.is_empty() {
                            self.my_u_edges.as_mut().expect("myUEdges")[ii][jj] = a_list[0].clone();
                        }
                    }
                }
            }

            // Ensure Same Parameter on U-edges
            let u_edges = self.my_u_edges.as_ref().expect("myUEdges").clone();
            for ii in 0..u_edges.len() {
                // if (mySec->IsUClosed() && ii == myUEdges->UpperRow()) continue;
                if self.my_sec.borrow().base().is_uclosed() && ii == u_edges.len() - 1 {
                    continue;
                }
                for jj in 0..u_edges[ii].len() {
                    let an_edge = u_edges[ii][jj].clone();
                    if an_edge.is_null() || brep.edge(an_edge.clone()).degenerated {
                        continue;
                    }
                    // TopoDS_Face Face1, Face2; int i1 = ii - 1, i2 = ii;
                    let mut face1 = Shape::null();
                    let mut face2 = Shape::null();
                    // (the OCCT 1-based rows ii/ii+1 map to the rcad 0-based
                    // rows of the myFaces table)
                    let mut i1 = ii; // OCCT ii - 1
                    let mut i2 = ii + 1; // OCCT ii
                    let upper_row = self.my_faces.as_ref().expect("myFaces").len();
                    if i1 == 0 && self.my_sec.borrow().base().is_uclosed() {
                        // i1 = myFaces->UpperRow();
                        i1 = upper_row;
                    }
                    if i2 > upper_row {
                        // i2 = 0;
                        i2 = 0;
                    }
                    if i1 != 0 {
                        let a_shape1 = self.my_faces.as_ref().expect("myFaces")[i1 - 1][jj]
                            .clone();
                        if matches!(a_shape1.data.as_ref(), TShape::Face(_)) {
                            face1 = a_shape1;
                        }
                    }
                    if i2 != 0 {
                        let a_shape2 = self.my_faces.as_ref().expect("myFaces")[i2 - 1][jj]
                            .clone();
                        if matches!(a_shape2.data.as_ref(), TShape::Face(_)) {
                            face2 = a_shape2;
                        }
                    }
                    if !face1.is_null() && !face2.is_null() {
                        // CorrectSameParameter(anEdge, Face1, Face2);
                        super::brep_fill_sweep_b::correct_same_parameter(
                            brep, &an_edge, &face1, &face2,
                        );
                    }
                }
            }

            for ii in 1..=nb_law {
                for jj in 1..=nb_path {
                    let a_face = self.my_faces.as_ref().expect("myFaces")[ii as usize - 1]
                        [(jj - 1) as usize]
                        .clone();
                    if !a_face.is_null() && matches!(a_face.data.as_ref(), TShape::Face(_)) {
                        // BB.Add(myTapes->ChangeValue(ii), aFace);
                        let tape =
                            self.my_tapes.as_ref().expect("myTapes")[ii as usize - 1].clone();
                        bb.add_to_shell(brep, tape, a_face);
                    }
                }
            }

            // Is it Closed ?
            if self.my_loc.borrow().base().is_closed(brep)
                && self.my_sec.borrow().base().is_uclosed()
            {
                // Check
                let mut closed = true;
                // TopExp::MapShapesAndAncestors(shell, EDGE, FACE, EFmap);
                let efmap = map_shapes_and_ancestors_edge_face(&shell);
                let mut iedge = 1usize;
                while iedge <= efmap.len() && closed {
                    let (the_edge, faces) = &efmap[iedge - 1];
                    if brep.edge(the_edge.clone()).degenerated {
                        iedge += 1;
                        continue;
                    }
                    // closed = (EFmap(iedge).Extent() > 1);
                    closed = faces.len() > 1;
                    iedge += 1;
                }
                // shell.Closed(closed);
                {
                    let sd = brep.shell_mut(shell.clone());
                    if closed {
                        sd.flags |= tshape_flags::CLOSED;
                    } else {
                        sd.flags &= !tshape_flags::CLOSED;
                    }
                }
            }
            // myShape = shell;
            self.my_shape = shell;
        }
    }

    /// OCCT PerformCorner (cxx L3587-3906).
    pub fn perform_corner(
        &mut self,
        brep: &mut BRep,
        index: i32,
        transition: BRepFillTransitionStyle,
        bounds: &ShapeHArray2,
    ) -> bool {
        // if (Transition == BRepFill_Modified) return true; // Do nothing.
        if transition == BRepFillTransitionStyle::Modified {
            return true;
        }

        // const double anAngularTol = 0.025;
        let an_angular_tol = 0.025;

        // OCCT BRepFill_Right maps to the rcad RightCorner variant.
        let mut the_transition = transition;
        let mut is_tangent = false;
        let mut f = 0.0f64;
        let mut l = 0.0f64;
        let i1;
        let i2;
        let mut p1;
        let p2;
        let mut t1 = DVec3::ZERO;
        let mut t2 = DVec3::ZERO;
        let mut tang;
        let mut sortant = DVec3::ZERO;

        if index > 1 {
            i1 = index - 1;
            i2 = index;
        } else {
            i1 = self.my_loc.borrow().base().nb_law();
            i2 = 1;
        }

        // Construct an axis supported by the bissectrice
        // myLoc->Law(I1)->GetDomain(F, L);
        // myLoc->Law(I1)->GetCurve()->D1(L, P1, T1); T1.Normalize();
        {
            let law = self.my_loc.borrow().base().law(i1);
            law.borrow().get_domain(&mut f, &mut l);
            let curve = law.borrow().get_curve().expect("GetCurve");
            p1 = curve.point_at(l);
            t1 = curve.derivative_at(l).normalize_or_zero();
        }
        {
            let law = self.my_loc.borrow().base().law(i2);
            law.borrow().get_domain(&mut f, &mut l);
            let curve = law.borrow().get_curve().expect("GetCurve");
            p2 = curve.point_at(f);
            t2 = curve.derivative_at(f).normalize_or_zero();
        }

        if gp_vec_angle(t1, t2) < self.my_ang_min {
            is_tangent = true;
            // gp_Vec t1, t2, V; gp_Mat M;
            let mut m1 = crate::geomalgo::geomfill::gp_mat::GpMat::identity();
            let mut m2 = crate::geomalgo::geomfill::gp_mat::GpMat::identity();
            let mut v = DVec3::ZERO;
            {
                let law = self.my_loc.borrow().base().law(i1);
                law.borrow().get_domain(&mut f, &mut l);
                law.borrow().d0(l, &mut m1, &mut v);
            }
            // t1 = M.Column(3);
            let mt1 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m1, 3);
            {
                let law = self.my_loc.borrow().base().law(i2);
                law.borrow().get_domain(&mut f, &mut l);
                law.borrow().d0(l, &mut m2, &mut v);
            }
            let mt2 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m2, 3);

            if gp_vec_angle(mt1, mt2) < self.my_ang_min {
                // "BRepFill_Sweep::PerformCorner : This is not a corner !"
                return true;
            }
            // Sortant = t2 - t1;
            sortant = mt2 - mt1;
        }

        if gp_vec_angle(t1, t2) >= std::f64::consts::PI - an_angular_tol {
            return false;
        }
        // if ((TheTransition == BRepFill_Right) && (T1.Angle(T2) > myAngMax))
        if the_transition == BRepFillTransitionStyle::RightCorner
            && gp_vec_angle(t1, t2) > self.my_ang_max
        {
            // TheTransition = BRepFill_Round;
            the_transition = BRepFillTransitionStyle::Round;
        }

        // Tang = T1 + T2; // Average direction
        tang = t1 + t2;
        // gp_Dir NormalOfBisPlane = Tang;
        let normal_of_bis_plane = tang.normalize_or_zero();
        // gp_Vec anIntersectPointCrossDirection = T1.Crossed(T2);
        let an_intersect_point_cross_direction = t1.cross(t2);
        if is_tangent {
            // Sortant -= Tang.Dot(Tang) * Tang;
            sortant -= tang.dot(tang) * tang;
        } else {
            // Sortant = T2 - T1; // Direction input
            sortant = t2 - t1;
            // Sortant *= -1; // " " output
            sortant *= -1.0;
            // Tang -= (Tang.Dot(T2)) * T2;
            tang -= tang.dot(t2) * t2;
        }

        // P1.BaryCenter(0.5, P2, 0.5);
        p1 = (1.0 - 0.5) * p1 + 0.5 * p2;
        // gp_Dir N(Sortant); gp_Dir Dx(Tang);
        let n = sortant.normalize_or_zero();
        let dx = tang.normalize_or_zero();

        // gp_Ax2 Axe(P1, N, Dx); gp_Ax2 AxeOfBisPlane(P1, NormalOfBisPlane);
        let axe = GpAx2::new(p1, n, dx);
        let axe_of_bis_plane = GpAx2::from_normal(p1, normal_of_bis_plane);

        // Construct 2 intersecting Shells
        // UEdges = new (1, mySec->NbLaw() + 1, 1, myLoc->NbLaw());
        // UEdges->ChangeArray2() = myUEdges->Array2();
        let u_edges: ShapeHArray2 = self.my_u_edges.as_ref().expect("myUEdges").clone();

        // aFaces = new (myFaces->LowerRow(), myFaces->UpperRow(), 1, 2);
        let my_faces = self.my_faces.as_ref().expect("myFaces").clone();
        let mut a_faces: ShapeHArray2 = vec![vec![Shape::null(); 2]; my_faces.len()];
        // Translate(myFaces, I1, aFaces, 1);
        translate(&my_faces, i1, &mut a_faces, 1);
        // Translate(myFaces, I2, aFaces, 2);
        translate(&my_faces, i2, &mut a_faces, 2);

        // aUEdges = new (myUEdges->LowerRow(), myUEdges->UpperRow(), 1, 2);
        let mut a_u_edges: ShapeHArray2 = vec![vec![Shape::null(); 2]; u_edges.len()];
        // Translate(myUEdges, I1, aUEdges, 1);
        translate(&u_edges, i1, &mut a_u_edges, 1);
        // Translate(myUEdges, I2, aUEdges, 2);
        translate(&u_edges, i2, &mut a_u_edges, 2);

        // gp_Vec aNormal = T2 + T1; TopoDS_Face aPlaneF;
        let a_normal = t2 + t1;
        let mut a_plane_f = Shape::null();

        if a_normal.length() > GP_RESOLUTION {
            // gp_Pln pl(P1, gp_Dir(aNormal)); BRepLib_MakeFace aFMaker(pl);
            let a_dir = a_normal.normalize_or_zero();
            // (the OCCT gp_Pln default X direction is solver-picked)
            let u_dir = if a_dir.x.abs() < 0.9 {
                DVec3::X.cross(a_dir)
            } else {
                DVec3::Y.cross(a_dir)
            }
            .normalize_or_zero();
            let pl = Surface3::Plane(Plane {
                origin: p1,
                normal: a_dir,
                u_dir,
                v_dir: a_dir.cross(u_dir).normalize_or_zero(),
            });
            let mut a_fmaker = BRepBuilder::new();
            let a_fmaker_face = a_fmaker.make_face(brep, Some(pl), Shape::null());
            // if (aFMaker.Error() == BRepLib_FaceDone) — the rcad make_face
            // is infallible at this layer.
            {
                a_plane_f = a_fmaker_face;
                // BRep_Builder aBB; aBB.UpdateFace(aPlaneF, PConfusion*10);
                BRepBuilder::new().update_face_tolerance(
                    brep,
                    a_plane_f.clone(),
                    p_confusion() * 10.0,
                );
            }
        }

        // BRepFill_TrimShellCorner aTrim(aFaces, Transition, AxeOfBisPlane,
        //                                anIntersectPointCrossDirection);
        let mut a_trim =
            crate::brep_fill::brep_fill_trim_shell_corner::BRepFillTrimShellCorner::new(
                &a_faces,
                transition,
                (
                    axe_of_bis_plane.location,
                    axe_of_bis_plane.direction,
                    axe_of_bis_plane.x_direction,
                ),
                an_intersect_point_cross_direction,
            );
        // aTrim.AddBounds(Bounds);
        a_trim.add_bounds(bounds);
        // aTrim.AddUEdges(aUEdges);
        a_trim.add_u_edges(&a_u_edges);
        // aTrim.AddVEdges(myVEdges, Index);
        let my_v_edges = self.my_v_edges.as_ref().expect("myVEdges").clone();
        a_trim.add_v_edges(&my_v_edges, index);
        // aTrim.Perform();
        a_trim.perform();

        if a_trim.is_done() {
            let nb_law = self.my_sec.borrow().base().nb_law();
            for ii in 1..=nb_law {
                // aTrim.Modified(myVEdges->Value(ii, Index), listmodif);
                let modified = a_trim.modified(
                    &self.my_v_edges.as_ref().expect("myVEdges")[ii as usize - 1]
                        [(index - 1) as usize],
                );
                if modified.is_empty() {
                    // TopoDS_Edge NullEdge; myVEdges->SetValue(ii, Index, NullEdge);
                    self.my_v_edges.as_mut().expect("myVEdges")[ii as usize - 1]
                        [(index - 1) as usize] = Shape::null();
                } else {
                    // myVEdges->SetValue(ii, Index, listmodif.First());
                    self.my_v_edges.as_mut().expect("myVEdges")[ii as usize - 1]
                        [(index - 1) as usize] = modified[0].clone();
                }
            }

            let mut iit = 0i32;
            while iit < 2 {
                // int II = (iit == 0) ? I1 : I2;
                let ii2 = if iit == 0 { i1 } else { i2 };

                for ii in 1..=nb_law {
                    // aTrim.Modified(myFaces->Value(ii, II), listmodif);
                    let modified = a_trim.modified(
                        &self.my_faces.as_ref().expect("myFaces")[ii as usize - 1]
                            [(ii2 - 1) as usize],
                    );
                    if !modified.is_empty() {
                        // myFaces->SetValue(ii, II, listmodif.First());
                        self.my_faces.as_mut().expect("myFaces")[ii as usize - 1]
                            [(ii2 - 1) as usize] = modified[0].clone();
                    }
                }

                let u_edges = self.my_u_edges.as_ref().expect("myUEdges").clone();
                for ii in 0..u_edges.len() {
                    // aTrim.Modified(myUEdges->Value(ii, II), listmodif);
                    let modified = a_trim.modified(&u_edges[ii][(ii2 - 1) as usize]);
                    if !modified.is_empty() {
                        // myUEdges->SetValue(ii, II, listmodif.First());
                        self.my_u_edges.as_mut().expect("myUEdges")[ii][(ii2 - 1) as usize] =
                            modified[0].clone();
                    }
                }
                iit += 1;
            }
        } else if (the_transition == BRepFillTransitionStyle::RightCorner)
            || a_trim.has_section()
        {
            // "Fail of TrimCorner"
            return true; // Nothing is touched
        }

        if self.my_sec.borrow().base().is_uclosed() {
            let nb_law = self.my_sec.borrow().base().nb_law();
            // myUEdges->SetValue(1, I1, myUEdges->Value(mySec->NbLaw() + 1, I1));
            let tail1 = self.my_u_edges.as_ref().expect("myUEdges")[nb_law as usize]
                [(i1 - 1) as usize]
                .clone();
            self.my_u_edges.as_mut().expect("myUEdges")[0][(i1 - 1) as usize] = tail1;
            // myUEdges->SetValue(1, I2, myUEdges->Value(mySec->NbLaw() + 1, I2));
            let tail2 = self.my_u_edges.as_ref().expect("myUEdges")[nb_law as usize]
                [(i2 - 1) as usize]
                .clone();
            self.my_u_edges.as_mut().expect("myUEdges")[0][(i2 - 1) as usize] = tail2;
        }

        // if (TheTransition == BRepFill_Round)
        if the_transition == BRepFillTransitionStyle::Round {
            // Filling
            let mut list1: Vec<Shape> = Vec::new();
            let mut list2: Vec<Shape> = Vec::new();
            let mut bord1 = Shape::null();
            let mut bord2 = Shape::null();
            let mut bord_first = Shape::null();
            let mut has_filling = false;
            let mut ff = Shape::null();
            let nb_law = self.my_sec.borrow().base().nb_law();
            for ii in 1..=nb_law {
                // KeepEdge(myFaces->Value(ii, I1), Bounds->Value(ii, 1), list1);
                let f1 = self.my_faces.as_ref().expect("myFaces")[ii as usize - 1]
                    [(i1 - 1) as usize]
                    .clone();
                let b1 = bounds[ii as usize - 1][0].clone();
                list1 = super::brep_fill_sweep_b::keep_edge(brep, &f1, &b1);
                // KeepEdge(myFaces->Value(ii, I2), Bounds->Value(ii, 2), list2);
                let f2 = self.my_faces.as_ref().expect("myFaces")[ii as usize - 1]
                    [(i2 - 1) as usize]
                    .clone();
                let b2 = bounds[ii as usize - 1][1].clone();
                list2 = super::brep_fill_sweep_b::keep_edge(brep, &f2, &b2);
                // if (list1.Extent() == list2.Extent())
                if list1.len() == list2.len() {
                    // (the OCCT parallel It1/It2 list walk maps to the index
                    // walk over the two rcad lists)
                    for it in 0..list1.len() {
                        // TopoDS_Edge E = TopoDS::Edge(It1.Value());
                        let e = list1[it].clone();
                        if has_filling {
                            // Transversal choice of constraints
                            let (vf, _vl) =
                                crate::brep_fill::brep_fill_pipe::top_exp_vertices(&e);
                            let mut e1 = Shape::null();
                            let mut e2 = Shape::null();
                            // if (!Bord1.IsNull() && TopExp::CommonVertex(E, Bord1, VC))
                            if !bord1.is_null() {
                                if let Some(vc) =
                                    crate::fillet::chfi3d_builder_0::topexp_common_vertex(&e, &bord1)
                                {
                                    if vc.is_same(&vf) {
                                        e1 = bord1.clone();
                                    } else {
                                        e2 = bord1.clone();
                                    }
                                }
                            }
                            if !bord2.is_null() {
                                if let Some(vc) =
                                    crate::fillet::chfi3d_builder_0::topexp_common_vertex(&e, &bord2)
                                {
                                    if vc.is_same(&vf) {
                                        e1 = bord2.clone();
                                    } else {
                                        e2 = bord2.clone();
                                    }
                                }
                            }
                            if !bord_first.is_null() {
                                if let Some(vc) = crate::fillet::chfi3d_builder_0::topexp_common_vertex(
                                    &e, &bord_first,
                                ) {
                                    if vc.is_same(&vf) {
                                        e1 = bord_first.clone();
                                    } else {
                                        e2 = bord_first.clone();
                                    }
                                }
                            }
                            // Bord1 = E1; Bord2 = E2;
                            bord1 = e1;
                            bord2 = e2;
                        }

                        // Filling
                        let f1 = self.my_faces.as_ref().expect("myFaces")[ii as usize - 1]
                            [(i1 - 1) as usize]
                            .clone();
                        let f2 = self.my_faces.as_ref().expect("myFaces")[ii as usize - 1]
                            [(i2 - 1) as usize]
                            .clone();
                        let mut bord1_mut = bord1.clone();
                        let mut bord2_mut = bord2.clone();
                        let b_ok = super::brep_fill_sweep_b2::filling(
                            brep,
                            &list1[it],
                            &f1,
                            &list2[it],
                            &f2,
                            &mut self.my_v_edges_modified,
                            self.my_tol3d,
                            &axe,
                            t1,
                            &mut bord1_mut,
                            &mut bord2_mut,
                            &mut ff,
                        );
                        bord1 = bord1_mut;
                        bord2 = bord2_mut;

                        if b_ok {
                            // myAuxShape.Append(FF);
                            self.my_aux_shape.push(ff.clone());
                            let mut bb = BRepBuilder::new();
                            // TopoDS_Shape aVshape = myVEdges->Value(ii, I2);
                            let a_vshape = self.my_v_edges.as_ref().expect("myVEdges")
                                [ii as usize - 1][(i2 - 1) as usize]
                                .clone();
                            // TopoDS_Compound aCompound; BB.MakeCompound(aCompound);
                            let a_compound = bb.make_compound(brep, vec![]);
                            if !a_vshape.is_null() {
                                // BB.Add(aCompound, aVshape);
                                bb.add_to_compound(brep, a_compound.clone(), a_vshape);
                            }
                            // BB.Add(aCompound, FF);
                            bb.add_to_compound(brep, a_compound.clone(), ff.clone());
                            // myVEdges->ChangeValue(ii, I2) = aCompound;
                            self.my_v_edges.as_mut().expect("myVEdges")[ii as usize - 1]
                                [(i2 - 1) as usize] = a_compound;

                            // BB.Add(myTapes->ChangeValue(ii), FF);
                            let tape = self.my_tapes.as_ref().expect("myTapes")
                                [ii as usize - 1]
                                .clone();
                            bb.add_to_shell(brep, tape, ff.clone());
                            // HasFilling = true;
                            has_filling = true;
                        }
                        if ii == 1 {
                            // BordFirst = Bord1;
                            bord_first = bord1.clone();
                        }
                    }
                }
                // else "PerformCorner : Unsymmetry of free border"
            }
        }
        true
    }

    /// OCCT RebuildTopOrBottomEdge (cxx L4047-4151) — rebuild v-iso edge of
    /// top or bottom section inserting new 3d and 2d curves taken from
    /// swept surfaces.
    pub fn rebuild_top_or_bottom_edge(
        &self,
        brep: &mut BRep,
        a_new_edge: &Shape,
        an_edge: &mut Shape,
        reversed_edges: &mut HashSet<ShapeKey>,
    ) {
        let mut bb = BRepBuilder::new();

        // occ::handle<Geom_Curve> aNewCurve = BRep_Tool::Curve(aNewEdge, fpar, lpar);
        let new_ed = brep.edge(a_new_edge.clone());
        let a_new_curve = new_ed
            .curve
            .clone()
            .expect("RebuildTopOrBottomEdge: new edge without 3d curve");
        let (mut fpar, mut lpar) = (new_ed.range[0], new_ed.range[1]);

        // bool ToReverse = false;
        let mut to_reverse = false;
        // bool IsDegen = BRep_Tool::Degenerated(aNewEdge);
        let is_degen = new_ed.degenerated;
        if is_degen {
            // BRep_Tool::Range(aNewEdge, fpar, lpar);
            // (the rcad range read above carries the edge range)
        } else {
            // TopExp::Vertices(anEdge, V1, V2);
            let (v1, v2) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(&an_edge.clone());
            if !v1.is_same(&v2) {
                // TopExp::Vertices(aNewEdge, NewV1, NewV2);
                let (new_v1, new_v2) =
                    crate::brep_fill::brep_fill_pipe::top_exp_vertices(a_new_edge);
                // V1.Location(Identity);
                let mut v1 = v1;
                v1.location = 0;
                if !v1.is_same(&new_v1) {
                    if v1.is_same(&new_v2) {
                        // ToReverse = true;
                        to_reverse = true;
                    } else {
                        // gp_Pnt Pnt1 = BRep_Tool::Pnt(V1);
                        let pnt1 = brep.vertex(v1.clone()).point;
                        // gp_Pnt NewPnt1 = BRep_Tool::Pnt(NewV1);
                        let new_pnt1 = brep.vertex(new_v1.clone()).point;
                        // double TolSum = Tolerance(V1) + Tolerance(NewV1);
                        let tol_sum = brep.vertex(v1.clone()).tolerance
                            + brep.vertex(new_v1.clone()).tolerance;
                        // if (!Pnt1.IsEqual(NewPnt1, TolSum)) ToReverse = true;
                        if pnt1.distance(new_pnt1) > tol_sum {
                            to_reverse = true;
                        }
                    }
                }
            } else {
                // double OldFirst, OldLast; OldCurve = BRep_Tool::Curve(anEdge, ...);
                let old_ed = brep.edge(an_edge.clone());
                let old_curve = old_ed
                    .curve
                    .clone()
                    .expect("RebuildTopOrBottomEdge: edge without 3d curve");
                let (old_first, old_last) = (old_ed.range[0], old_ed.range[1]);
                // OldCurve->D1(0.5 * (OldFirst + OldLast), MidPnt, OldD1);
                let old_d1 = old_curve.derivative_at(0.5 * (old_first + old_last));
                // aNewCurve->D1(0.5 * (fpar + lpar), MidPnt, NewD1);
                let new_d1 = a_new_curve.derivative_at(0.5 * (fpar + lpar));
                // if (OldD1 * NewD1 < 0.) ToReverse = true;
                if old_d1.dot(new_d1) < 0.0 {
                    to_reverse = true;
                }
            }
        }

        // anEdge.Location(Identity);
        an_edge.location = 0;
        // const BRep_TEdge& TEdge = anEdge.TShape();
        // TEdge->Tolerance(BRep_Tool::Tolerance(aNewEdge));
        super::brep_fill_sweep_b::set_edge_tolerance(
            brep,
            an_edge,
            brep.edge(a_new_edge.clone()).tolerance,
        );
        // BB.Range(anEdge, fpar, lpar);
        bb.set_edge_range(brep, an_edge.clone(), fpar, lpar);
        // BB.UpdateEdge(anEdge, aNewCurve, Precision::Confusion());
        {
            let ed = brep.edge_mut_inplace(an_edge.clone());
            ed.curve = Some(a_new_curve.clone());
            ed.tolerance = ed.tolerance.max(CONFUSION);
        }
        // const BRep_TEdge& TE = aNewEdge.TShape();
        // const NCollection_List<...>& lcr = TE->Curves();
        {
            let new_ed = brep.edge(a_new_edge.clone());
            let lcr = new_ed.representations.clone();
            for curve_rep in lcr {
                // if (CurveRep->IsCurveOnSurface())
                match curve_rep {
                    rcad_kernel::topo::topods::CurveRepresentation::CurveOnSurface {
                        face,
                        pcurve,
                        ..
                    } => {
                        // occ::handle<Geom_Surface> aSurf = GC->Surface();
                        // TopLoc_Location aLoc = aNewEdge.Location() * GC->Location();
                        let a_loc = face.1;
                        // (the surface value is resolved before the mutable
                        // update — the rcad pool borrow split)
                        let a_surf = brep_face_surface(brep, face.0);
                        if let Some(a_surf) = a_surf {
                            // BB.UpdateEdge(anEdge, aPCurve, aSurf, aLoc,
                            //               Precision::Confusion());
                            super::brep_fill_sweep_b::update_edge_pcurve_on_surf(
                                brep, an_edge, &pcurve, &a_surf, a_loc, CONFUSION,
                            );
                        }
                    }
                    rcad_kernel::topo::topods::CurveRepresentation::CurveOnClosedSurface {
                        face,
                        pcurve1,
                        pcurve2,
                        ..
                    } => {
                        let a_loc = face.1;
                        let a_surf = brep_face_surface(brep, face.0);
                        if let Some(a_surf) = a_surf {
                            // the closed-surface rep carries both
                            // pcurves (BRep_GCurve::PCurve/PCurve2)
                            super::brep_fill_sweep_b::update_edge_pcurve2_on_surf(
                                brep, an_edge, &pcurve1, &pcurve2, &a_surf, a_loc, CONFUSION,
                            );
                        }
                    }
                    _ => {}
                }
            }
        }

        // anEdge.Free(true);
        {
            let ed = brep.edge_mut_inplace(an_edge.clone());
            ed.flags |= tshape_flags::FREE;
        }

        // TopExp::Vertices(anEdge, V1, V2);
        let (v1, v2) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(&an_edge.clone());

        // TopoDS_Shape anEdgeFORWARD = anEdge.Oriented(TopAbs_FORWARD);
        let mut an_edge_forward = an_edge.clone();
        an_edge_forward.orientation = Orientation::Forward;

        // BB.Remove(anEdgeFORWARD, V1); BB.Remove(anEdgeFORWARD, V2);
        bb.remove_from_edge(brep, an_edge_forward.clone(), v1.clone());
        bb.remove_from_edge(brep, an_edge_forward.clone(), v2.clone());

        // V1.Location(Identity); V2.Location(Identity);
        let mut v1 = v1;
        v1.location = 0;
        let mut v2 = v2;
        v2.location = 0;
        if to_reverse {
            // V2.Orientation(TopAbs_FORWARD); V1.Orientation(TopAbs_REVERSED);
            v2.orientation = Orientation::Forward;
            v1.orientation = Orientation::Reversed;
        }
        // BB.Add(anEdgeFORWARD, V1); BB.Add(anEdgeFORWARD, V2);
        bb.add_to_edge(brep, an_edge_forward.clone(), v1);
        bb.add_to_edge(brep, an_edge_forward.clone(), v2);

        if to_reverse {
            // anEdge.Reverse();
            an_edge.orientation = match an_edge.orientation {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
                o => o,
            };
            // ReversedEdges.Add(anEdge);
            reversed_edges.insert(ShapeKey(an_edge.ptr_id()));
        }

        // BB.Degenerated(anEdge, IsDegen);
        bb.set_edge_degenerated(brep, an_edge.clone(), is_degen);
    }
}
