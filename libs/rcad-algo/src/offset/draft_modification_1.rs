//! OCCT Draft_Modification_1.cxx — the private member functions of
//! Draft_Modification: InternalAdd (L116-420), Propagate (L424-704),
//! Perform (L708-1776), and the private NewSurface (L1780-1959) /
//! NewCurve (L1963-2043) geometry mappers.
//!
//! Perform (the single OCCT member) is carried as the entry plus the three
//! OCCT-commented sections — Calculate eventual faces (L723-824), Calculate
//! new edges (L826-1443, in offset/draft_modification_1_c.rs), Calculate new
//! vertices (L1445-1636, in draft_modification_1_c.rs) and the small loop of
//! validation/protection (L1639-1775, in draft_modification_1_c.rs) — the
//! split follows the single-file <2000-line rule only.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/
//!         Draft_Modification_1.cxx
//!
//! Architecture differences (continuing draft_modification.rs):
//! 5. Geom_Surface::Transformed(Loc.Transformation()) /
//!    Geom_Curve::Transformed(Loc) — the rcad locations are identity in the
//!    draft flows (the loc_ope_split_drafts.rs #12 precedent).
//! 6. The `S != NewS` comparison (the Geom_Surface handle identity) — the
//!    rcad surface_same stand-in (topods::surface_same).
//! 7. BRep_Builder::MakeFace(TheNewFace, S2, Tol) — the marker-wire
//!    add_tface_tol form (the loc_ope_split_drafts.rs #6 precedent); the
//!    BRep_Builder local is a pool-local BRep (arch. difference #4 of
//!    brep_offset_make_simple_offset.rs).
//! 8. BRep_Tool::Parameter / Range / Curve / Surface — the brep_algo::tool
//!    re-hosts.
//! 9. GeomAdaptor_Curve / GeomAdaptor_Surface — the
//!    bop::int_tools::bean_face_intersector::{BRepAdaptorCurve,
//!    BRepAdaptorSurface} vehicles; IntCurveSurface_HInter is the
//!    bean_face_intersector::IntCurveSurfaceHInter.

use glam::DVec3;
use rcad_kernel::geom::{
    Circle3, ConicalSurface, CylindricalSurface, Curve3, CurveEval, Line3, LinearExtrusionSurface,
    Plane, Surface3, SurfaceEval, TrimmedSurface,
};
use rcad_kernel::math::gp::Ax3;
use rcad_kernel::precision::{
    is_negative_infinite_value, is_positive_infinite_value, ANGULAR, CONFUSION,
};
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{surface_same, GeomAbsShape, Orientation, ShapeType};

use crate::brep_algo::tool::{builder_make_wire, brep_tool_curve, brep_tool_parameter, brep_tool_surface, explorer};

use super::draft_edge_info::DraftEdgeInfo;
use super::draft_error_status::DraftErrorStatus;
use super::draft_face_info::DraftFaceInfo;
use super::draft_modification::{
    brep_tool_continuity, brep_tools_is_really_closed, brep_tools_uv_bounds, DraftModification,
};
use super::draft_modification_1_b::{find_rotation, orientation, IntAnaIntConicQuad};
use super::draft_vertex_info::DraftVertexInfo;

/// OCCT gp_Dir::IsParallel(Other, AngularTolerance) (gp_Dir.cxx) — the
/// parallel test consumed by InternalAdd and Propagate.
pub(crate) fn dir_is_parallel(the_d1: DVec3, the_d2: DVec3, the_tol: f64) -> bool {
    let ang = the_d1.dot(the_d2).clamp(-1.0, 1.0).acos();
    ang <= the_tol || std::f64::consts::PI - ang <= the_tol
}

impl DraftModification {
    /// OCCT Draft_Modification::InternalAdd(F, Direction, Angle, NeutralPlane,
    /// Flag) (Draft_Modification_1.cxx L116-420).
    pub(crate) fn internal_add(
        &mut self,
        the_f: &Shape,
        the_direction: DVec3,
        the_angle: f64,
        the_neutral_plane: &Plane,
        the_flag: bool,
    ) -> bool {
        // OCCT L123-126.
        if self.fmap().contains(the_f) {
            return self.bad_shape().is_null();
        }

        // OCCT L128: oris = Orientation(myShape, F);
        let oris = orientation(self.my_shape(), the_f);
        // OCCT L129: TopLoc_Location Lo; — identity.
        // OCCT L132-133: S = BRep_Tool::Surface(F, Lo);
        //                S = down_cast<Geom_Surface>(S->Transformed(Lo.Transformation()));
        // (architecture difference #5: the location re-application is identity).
        let mut s = brep_tool_surface(the_f).expect("BRep_Tool::Surface");
        // OCCT L134-137.
        if let Surface3::Trimmed(t) = &s {
            s = t.basis.as_ref().clone();
        }
        // OCCT L138-139: NewS; theCircle.
        let mut new_s: Option<Surface3> = None;
        let mut the_circle: Option<Curve3> = None;

        // OCCT L141-198.
        let mut postponed = !the_flag;
        if postponed {
            let is_cyl_or_le =
                matches!(&s, Surface3::Cylinder(_) | Surface3::LinearExtrusion(_));
            if is_cyl_or_le {
                let mut cir: Circle3;
                if let Surface3::Cylinder(cyl) = &s {
                    // OCCT L151-155.
                    let axcyl_dir = cyl.axis;
                    // OCCT: Cir = ElSLib::CylinderVIso(cyl.Position(), cyl.Radius(), 0.)
                    // — the V=0 iso circle on the cylinder frame.
                    cir = Circle3 {
                        center: cyl.origin,
                        normal: cyl.axis,
                        x_dir: cyl.ref_dir,
                        y_dir: cyl.y_axis(),
                        radius: cyl.radius,
                    };
                    // OCCT L154-155: VV = Vec(cyl.Location(), NeutralPlane.Location());
                    // Cir.Translate(VV.Dot(axcyl.Direction()) * axcyl.Direction());
                    let vv = the_neutral_plane.origin - cyl.origin;
                    let d = vv.dot(axcyl_dir);
                    cir.center += axcyl_dir * d;
                } else {
                    // OCCT L157-190: the linear extrusion branch.
                    let sle = match &s {
                        Surface3::LinearExtrusion(x) => x,
                        _ => unreachable!(),
                    };
                    let mut cbas = sle.profile.as_ref().clone();
                    let the_dirextr = sle.direction;

                    // OCCT L163-166.
                    if let Curve3::Trimmed(t) = cbas {
                        cbas = t.curve.as_ref().clone();
                    }
                    // OCCT L167-183.
                    match &cbas {
                        Curve3::Circle(c) => {
                            cir = *c;
                            let dircir = cir.normal;
                            if !dir_is_parallel(the_direction, dircir, ANGULAR) {
                                self.bad_shape_mut().clone_from(the_f);
                                self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                                return false;
                            }
                        }
                        _ => {
                            self.bad_shape_mut().clone_from(the_f);
                            self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                            return false;
                        }
                    }

                    // OCCT L185-189.
                    let axis = the_neutral_plane; // gp_Ax3 Axis = NeutralPlane.Position();
                    let l = (cir.center - axis.origin).dot(axis.normal);
                    let cos_ = the_dirextr.dot(axis.normal);
                    let vv = the_dirextr * (l / cos_);
                    cir.center += vv;
                }
                // OCCT L192: theCircle = new Geom_Circle(Cir);
                the_circle = Some(Curve3::Circle(cir));
            } else {
                // OCCT L196: postponed = false;
                postponed = false;
            }
        }

        // OCCT L200-262.
        if !postponed {
            // OCCT L202: NewS = NewSurface(S, oris, Direction, Angle, NeutralPlane);
            new_s = self.new_surface_geom(&s, oris, the_direction, the_angle, the_neutral_plane);
            // OCCT L203-208.
            if new_s.is_none() {
                self.bad_shape_mut().clone_from(the_f);
                self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                return false;
            }
            // OCCT L210-261: the infinite-restriction guard on the
            // cylinder/cone results.
            let typs_cyl = matches!(new_s, Some(Surface3::Cylinder(_)));
            let typs_con = matches!(new_s, Some(Surface3::Cone(_)));
            if typs_cyl || typs_con {
                // OCCT L214-215: BRepTools::UVBounds(F, umin, umax, vmin, vmax);
                let bounds = brep_tools_uv_bounds(the_f);
                let (umin, umax, mut vmin, mut vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                let _ = (umin, umax); // the OCCT umin/umax are unused at this site
                // OCCT L216.
                if !is_negative_infinite_value(vmin) && !is_positive_infinite_value(vmax) {
                    let deltav = 10.0 * (vmax - vmin);
                    if typs_cyl {
                        // OCCT L219-223.
                        vmin -= deltav;
                        vmax += deltav;
                    } else {
                        // OCCT L224-258: the cone apex guard.
                        let co: ConicalSurface = match new_s {
                            Some(Surface3::Cone(c)) => c,
                            _ => unreachable!(),
                        };
                        let vapex = -co.radius / co.half_angle_rad.sin();
                        if vmin < vapex {
                            // vmax should not exceed Vapex
                            if vmax + deltav > vapex {
                                vmax = vapex;
                                vmin = vmin - 10.0 * (vmax - vmin);
                                // JAG debug to avoid apex
                                vmax -= CONFUSION;
                            } else {
                                vmin -= deltav;
                                vmax += deltav;
                            }
                        } else {
                            // Vapex <= vmin < vmax
                            if vmin - deltav < vapex {
                                vmin = vapex;
                                vmax += 10.0 * (vmax - vmin);
                                // JAG debug to avoid apex
                                vmin += CONFUSION;
                            } else {
                                vmin -= deltav;
                                vmax += deltav;
                            }
                        }
                    }
                    // OCCT L259.
                    new_s = Some(Surface3::Trimmed(TrimmedSurface::new(
                        new_s.expect("NewS"),
                        0.0,
                        2.0 * std::f64::consts::PI,
                        vmin,
                        vmax,
                    )));
                }
            }
        }

        // OCCT L264-273: if (postponed || S != NewS).
        let same_geometry = match &new_s {
            Some(ns) => surface_same(&s, ns),
            None => false,
        };
        if postponed || !same_geometry {
            // OCCT L266: Draft_FaceInfo FI(NewS, true); — the postponed case
            // carries the null geometry.
            let mut fi = DraftFaceInfo::new_with_geometry(new_s.as_ref(), true);
            // OCCT L267: FI.RootFace(curFace);
            fi.set_root_face(self.cur_face());
            // OCCT L268: myFMap.Add(F, FI);
            self.fmap_mut().add(the_f, fi);
            if postponed {
                // OCCT L271: myFMap.ChangeFromKey(F).ChangeCurve() = theCircle;
                *self.fmap_mut().change_from_key(the_f).change_curve() = the_circle;
            }
        }

        // OCCT L275-418: the edge loop.
        // OCCT L276: NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> MapOfE;
        let mut map_of_e: Vec<Shape> = Vec::new();
        // OCCT L275: TopExp_Explorer aExp(F, TopAbs_EDGE);
        let a_exp = explorer(the_f, ShapeType::Edge, ShapeType::Shape);
        let mut a_i = 0usize;
        // OCCT L277: while (aExp.More() && badShape.IsNull()).
        while a_i < a_exp.len() && self.bad_shape().is_null() {
            // OCCT L279: edg = TopoDS::Edge(aExp.Current());
            let edg = a_exp[a_i].clone();
            if !self.emap().contains(&edg) {
                // OCCT L282-284.
                let mut addedg = false;
                let mut addface = false;
                let mut other_f = Shape::null();
                // OCCT L286-290.
                if brep_tools_is_really_closed(&edg, the_f) {
                    addedg = true;
                    addface = false;
                } else {
                    // OCCT L292-308: find the other face containing the edge.
                    let mut nbother = 0;
                    for a_v in self.efmap_find(&edg) {
                        if !a_v.is_same(the_f) {
                            if other_f.is_null() {
                                other_f = a_v.clone(); // TopoDS::Face(it.Value())
                            }
                            nbother += 1;
                        }
                    }
                    // OCCT L309-313.
                    if nbother >= 2 {
                        self.bad_shape_mut().clone_from(&edg);
                        self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                    } else if !other_f.is_null()
                        && brep_tool_continuity(&edg, the_f, &other_f) >= GeomAbsShape::G1
                    {
                        // OCCT L314-318.
                        addface = true;
                        addedg = true;
                    } else if nbother == 0 {
                        // OCCT L319-322: (the empty branch).
                    }
                }
                // OCCT L324-415.
                if addedg {
                    if postponed {
                        // OCCT L326-329: myFMap.ChangeFromKey(F).Add(OtherF);
                        self.fmap_mut().change_from_key(the_f).add(&other_f);
                    }
                    // OCCT L330-333.
                    let c0 = brep_tool_curve(&edg).expect("BRep_Tool::Curve(edg)");
                    let mut c = c0.0; // C->Transformed(L) — identity (#5)
                    // OCCT L334-337.
                    if let Curve3::Trimmed(t) = &c {
                        c = t.curve.as_ref().clone();
                    }
                    // OCCT L338-339.
                    let mut new_c: Option<Curve3> = None;
                    let mut e_inf = DraftEdgeInfo::new_with_geometry(true);
                    if postponed {
                        // OCCT L342-343.
                        e_inf.add(the_f);
                        e_inf.add(&other_f);

                        // OCCT L345-366: find the fixed point.
                        match &c {
                            Curve3::Line(l) => {
                                // OCCT L355.
                                let ilipl = IntAnaIntConicQuad::new(l, the_neutral_plane, ANGULAR);
                                // OCCT L356-364.
                                if ilipl.is_done() && ilipl.nb_points() != 0 {
                                    e_inf.tangent(ilipl.point(1));
                                } else {
                                    self.bad_shape_mut().clone_from(&edg);
                                    self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                                }
                            }
                            _ => {
                                // OCCT L346-351: aLocalGeom null (not a line).
                                self.bad_shape_mut().clone_from(&edg);
                                self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                            }
                        }
                    } else {
                        // OCCT L369.
                        new_c = self.new_curve_geom(
                            &c,
                            &s,
                            oris,
                            the_direction,
                            the_angle,
                            the_neutral_plane,
                            the_flag,
                        );
                        // OCCT L370-374.
                        if new_c.is_none() {
                            self.bad_shape_mut().clone_from(&edg);
                            self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                        }
                    }

                    // OCCT L377-381.
                    if let Some(Curve3::Trimmed(t)) = &new_c {
                        new_c = Some(t.curve.as_ref().clone());
                    }
                    // OCCT L382: EInf.ChangeGeometry() = NewC;
                    *e_inf.change_geometry() = new_c;
                    // OCCT L384: EInf.RootFace(curFace);
                    e_inf.set_root_face(self.cur_face());
                    // OCCT L385: myEMap.Add(edg, EInf);
                    self.emap_mut().add(&edg, e_inf);
                    // OCCT L386: MapOfE.Add(edg);
                    if !map_of_e.iter().any(|e| e.is_same(&edg)) {
                        map_of_e.push(edg.clone());
                    }
                    // OCCT L387-414.
                    if addface {
                        // OCCT L389: bool Fl = Flag;
                        let fl = the_flag;
                        // OCCT L390-395.
                        let mut alocal_surface = brep_tool_surface(&other_f);
                        if let Some(Surface3::Trimmed(t)) = &alocal_surface {
                            alocal_surface = Some(t.basis.as_ref().clone());
                        }
                        // OCCT L396-398.
                        let typs_cyl_or_le = matches!(
                            &alocal_surface,
                            Some(Surface3::Cylinder(_)) | Some(Surface3::LinearExtrusion(_))
                        );
                        if typs_cyl_or_le {
                            // OCCT L400-411.
                            if self.fmap().contains(the_f) && !the_flag && !postponed {
                                // OCCT L404: myFMap.RemoveKey(F);
                                self.fmap_mut().remove_key(the_f);
                                // OCCT L405-409.
                                for itm in &map_of_e {
                                    self.emap_mut().remove_key(itm); // TopoDS::Edge(itm.Key())
                                }
                            }
                        }
                        // OCCT L413.
                        self.internal_add(&other_f, the_direction, the_angle, the_neutral_plane, fl);
                    }
                }
            }
            // OCCT L417: aExp.Next();
            a_i += 1;
        }
        // OCCT L419: return (badShape.IsNull());
        self.bad_shape().is_null()
    }

    /// OCCT Draft_Modification::Propagate() (Draft_Modification_1.cxx
    /// L424-704).
    pub(crate) fn propagate(&mut self) -> bool {
        // OCCT L427-430.
        if !self.bad_shape().is_null() {
            return false;
        }

        // OCCT L432-473: set all edges and vertices of modified faces.
        for i in 1..=self.fmap().extent() {
            // OCCT L441: const TopoDS_Face& Fc = myFMap.FindKey(i);
            let fc = self.fmap().find_key(i).clone();

            // OCCT L443-444: exploration of the edges of the face.
            let f_edges = explorer(&fc, ShapeType::Edge, ShapeType::Shape);
            for e in f_edges {
                // OCCT L447: E = TopoDS::Edge(editer.Current());
                // OCCT L449-453.
                if !self.emap().contains(&e) {
                    let e_inf = DraftEdgeInfo::new_with_geometry(true);
                    self.emap_mut().add(&e, e_inf);
                }
                // OCCT L454: myEMap.ChangeFromKey(E).Add(Fc);
                self.emap_mut().change_from_key(&e).add(&fc);

                // OCCT L456-458: exploration of the vertices of the edge.
                let e_vertices = explorer(&e, ShapeType::Vertex, ShapeType::Shape);
                for v in e_vertices {
                    // OCCT L460: V = TopoDS::Vertex(vtiter.Current());
                    // OCCT L461-465.
                    if !self.vmap().contains(&v) {
                        let v_inf = DraftVertexInfo::new();
                        self.vmap_mut().add(&v, v_inf);
                    }
                    // OCCT L467-468.
                    self.vmap_mut().change_from_key(&v).add(&e);
                    let prm = brep_tool_parameter(&v, &e);
                    *self.vmap_mut().change_from_key(&v).change_parameter(&e) = prm;
                }
            }
        }

        // OCCT L475-513: set edges containing modified vertices.
        for i in 1..=self.vmap().extent() {
            // OCCT L482: const TopoDS_Vertex& Vt = myVMap.FindKey(i);
            let vt = self.vmap().find_key(i).clone();

            // OCCT L484-485: exploration of the ancestors of the vertex —
            // all edges of the shape.
            let anc = explorer(self.my_shape(), ShapeType::Edge, ShapeType::Shape);
            for e in anc {
                // OCCT L490-500: found = Vt among the vertices of E.
                let mut found = false;
                let vtiter = explorer(&e, ShapeType::Vertex, ShapeType::Shape);
                for v_shape in vtiter {
                    if vt.is_same(&v_shape) {
                        found = true;
                        break;
                    }
                }
                // OCCT L501-509.
                if found {
                    if !self.emap().contains(&e) {
                        let e_inf = DraftEdgeInfo::new_with_geometry(false);
                        self.emap_mut().add(&e, e_inf);
                    }
                    self.vmap_mut().change_from_key(&vt).add(&e);
                    let prm = brep_tool_parameter(&vt, &e);
                    *self.vmap_mut().change_from_key(&vt).change_parameter(&e) = prm;
                }
            }
        }

        // OCCT L515-550: set faces containing modified edges.
        for i in 1..=self.emap().extent() {
            // OCCT L518: const TopoDS_Edge& Ed = myEMap.FindKey(i);
            let ed = self.emap().find_key(i).clone();
            // OCCT L519-520: the ancestor faces of Ed.
            for f in self.efmap_find(&ed) {
                // TopoDS::Face(it.Value())
                // OCCT L523-547.
                if !self.fmap().contains(&f) {
                    // OCCT L525-528: S = BRep_Tool::Surface(F, L);
                    // NewS = S->Transformed(L.Transformation()); — identity.
                    let s = brep_tool_surface(&f).expect("BRep_Tool::Surface");
                    let mut new_s = s.clone();
                    // OCCT L530-543.
                    let typs_cyl = matches!(s, Surface3::Cylinder(_));
                    let typs_con = matches!(s, Surface3::Cone(_));
                    if typs_cyl || typs_con {
                        let bounds = brep_tools_uv_bounds(&f);
                        let (umin, umax, mut vmin, mut vmax) =
                            (bounds[0], bounds[1], bounds[2], bounds[3]);
                        let _ = (umin, umax); // the OCCT umin/umax are unused here
                        if !is_negative_infinite_value(vmin) && !is_positive_infinite_value(vmax)
                        {
                            let deltav = 10.0 * (vmax - vmin);
                            vmin -= deltav;
                            vmax += deltav;
                            new_s = Surface3::Trimmed(TrimmedSurface::new(
                                new_s,
                                0.0,
                                2.0 * std::f64::consts::PI,
                                vmin,
                                vmax,
                            ));
                        }
                    }
                    // OCCT L545-546.
                    let f_inf = DraftFaceInfo::new_with_geometry(Some(&new_s), false);
                    self.fmap_mut().add(&f, f_inf);
                }
                // OCCT L548: myEMap.ChangeFromKey(Ed).Add(F);
                self.emap_mut().change_from_key(&ed).add(&f);
            }
        }

        // OCCT L552-702: try to add faces for free borders (JAG 09.11.95).
        for i in 1..=self.emap().extent() {
            // OCCT L556: Draft_EdgeInfo& Einf = myEMap.ChangeFromIndex(i);
            let einf_needs_face = {
                let einf = self.emap_mut().change_from_index(i);
                einf.new_geometry() && einf.geometry().is_none() && einf.second_face().is_null()
            };
            if einf_needs_face {
                // OCCT L560-563: S1 = BRep_Tool::Surface(Einf.FirstFace(), Loc);
                // S1 = S1->Transformed(Loc.Transformation()); — identity.
                let first_face = self.emap().find_from_index(i).first_face().clone();
                let s1 = brep_tool_surface(&first_face).expect("BRep_Tool::Surface");
                let mut s2: Option<Surface3> = None;

                // OCCT L565-572: C = BRep_Tool::Curve(EK, Loc, f, l);
                let ek = self.emap().find_key(i).clone(); // OCCT: FindKey(i)
                let c0 = brep_tool_curve(&ek).expect("BRep_Tool::Curve(EK)");
                let mut c = c0.0; // identity location re-application
                if let Curve3::Trimmed(t) = &c {
                    c = t.curve.as_ref().clone();
                }

                // OCCT L573-603: if (!S1->IsKind(Plane)).
                let s1_is_plane = matches!(&s1, Surface3::Plane(_));
                if !s1_is_plane {
                    // OCCT L575-579: C IsKind(Geom_Conic).
                    let conic_plane: Option<Plane> = match &c {
                        Curve3::Circle(cir) => Some(conic_position_plane(cir.center, cir.normal)),
                        Curve3::Ellipse(el) => Some(conic_position_plane(el.center, el.normal)),
                        Curve3::Hyperbola(h) => Some(conic_position_plane(h.center, h.normal)),
                        Curve3::Parabola(p) => Some(conic_position_plane(p.vertex, p.normal)),
                        _ => None,
                    };
                    if let Some(the_pl) = conic_plane {
                        // OCCT L578: S2 = new Geom_Plane(thePl);
                        s2 = Some(Surface3::Plane(the_pl));
                    } else if matches!(&c, Curve3::Line(_)) {
                        // OCCT L580-595: the axis of S1.
                        let axis = elementary_surface_axis(&s1);
                        // OCCT L593: they = gp_Vec(axis.Location(), C->Value(0.));
                        let they = c.point_at(0.0) - axis.location;
                        // OCCT L594: axz = axis.Direction().Crossed(they);
                        let axz = axis.direction.cross(they);
                        // OCCT L595.
                        let a3 = Ax3::from_pnt_n_vx(axis.location, axz, axis.direction);
                        s2 = Some(Surface3::Plane(ax3_plane(&a3)));
                    } else {
                        // OCCT L597-602.
                        self.bad_shape_mut().clone_from(&ek);
                        self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                        break; // leave from for
                    }
                } else {
                    // OCCT L604-685: on the plane.
                    for j in 1..=self.vmap().extent() {
                        // OCCT L607-615: find EK among the edges of Vinf
                        // (the loop leaves the member iterator at EK).
                        let more1 = {
                            let vinf = self.vmap_mut().change_from_index(j);
                            vinf.init_edge_iterator();
                            while vinf.more_edge() {
                                if vinf.edge().is_same(&ek) {
                                    break;
                                }
                                vinf.next_edge();
                            }
                            vinf.more_edge()
                        };
                        // OCCT L616: if (Vinf.MoreEdge()).
                        if more1 {
                            // OCCT L617-631: walk the edges again; break at
                            // the edge whose faces do not carry
                            // Einf.FirstFace (the iterator position then
                            // holds that edge).
                            let einf_first_face =
                                self.emap().find_from_key(&ek).first_face().clone();
                            let vinf_edges: Vec<Shape> = {
                                let vinf = self.vmap_mut().change_from_index(j);
                                vinf.init_edge_iterator();
                                let mut v = Vec::new();
                                while vinf.more_edge() {
                                    v.push(vinf.edge());
                                    vinf.next_edge();
                                }
                                v
                            };
                            let mut hit_edge: Option<Shape> = None;
                            for edg in &vinf_edges {
                                if !edg.is_same(&ek) {
                                    let ei = self.emap().find_from_key(edg);
                                    if !ei.first_face().is_same(&einf_first_face)
                                        && (ei.second_face().is_null()
                                            || !ei.second_face().is_same(&einf_first_face))
                                    {
                                        hit_edge = Some(edg.clone());
                                        break;
                                    }
                                }
                            }
                            // OCCT L631: if (Vinf.MoreEdge()).
                            if let Some(found_edge) = &hit_edge {
                                // OCCT L633-640.
                                let (c2_raw, _f2, _l2) =
                                    brep_tool_curve(found_edge)
                                        .expect("BRep_Tool::Curve(Vinf.Edge())");
                                let mut c2 = c2_raw; // identity location re-application
                                if let Curve3::Trimmed(t) = &c2 {
                                    c2 = t.curve.as_ref().clone();
                                }
                                let direc: DVec3;
                                let hcur: Curve3;
                                if matches!(&c, Curve3::Line(_)) {
                                    // OCCT L641-645.
                                    direc = match &c {
                                        Curve3::Line(l) => l.direction,
                                        _ => unreachable!(),
                                    };
                                    hcur = c2;
                                } else if matches!(&c2, Curve3::Line(_)) {
                                    // OCCT L646-650.
                                    direc = match &c2 {
                                        Curve3::Line(l) => l.direction,
                                        _ => unreachable!(),
                                    };
                                    hcur = c.clone();
                                } else {
                                    // OCCT L651-656.
                                    self.bad_shape_mut().clone_from(&ek);
                                    self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                                    break; // leave from while
                                }
                                // OCCT L657-673: the SLE classification.
                                match &hcur {
                                    Curve3::Line(l) => {
                                        if !dir_is_parallel(l.direction, direc, ANGULAR) {
                                            // GeomAbs_Plane — the SLE::Plane vehicle.
                                            let normal = direc.cross(l.direction);
                                            let a3 =
                                                Ax3::from_pnt_n_vx(l.origin, normal, l.direction);
                                            s2 = Some(Surface3::Plane(ax3_plane(&a3)));
                                        } else {
                                            s2 = Some(Surface3::LinearExtrusion(
                                                LinearExtrusionSurface {
                                                    profile: Box::new(hcur.clone()),
                                                    direction: direc,
                                                },
                                            ));
                                        }
                                    }
                                    Curve3::Circle(cir) => {
                                        if dir_is_parallel(cir.normal, direc, ANGULAR) {
                                            // GeomAbs_Cylinder — the SLE::Cylinder
                                            // vehicle.
                                            s2 = Some(Surface3::Cylinder(CylindricalSurface {
                                                origin: cir.center,
                                                axis: direc,
                                                radius: cir.radius,
                                                ref_dir: cir.x_dir,
                                                y_dir: Some(cir.y_dir),
                                            }));
                                        } else {
                                            s2 = Some(Surface3::LinearExtrusion(
                                                LinearExtrusionSurface {
                                                    profile: Box::new(hcur.clone()),
                                                    direction: direc,
                                                },
                                            ));
                                        }
                                    }
                                    _ => {
                                        // OCCT default: the extrusion surface.
                                        s2 = Some(Surface3::LinearExtrusion(
                                            LinearExtrusionSurface {
                                                profile: Box::new(hcur.clone()),
                                                direction: direc,
                                            },
                                        ));
                                    }
                                }
                            } else {
                                // OCCT L676-680.
                                self.bad_shape_mut().clone_from(&ek);
                                self.set_err_stat(DraftErrorStatus::EdgeRecomputation);
                                break; // leave from while
                            }
                            // OCCT L681: break;
                            break;
                        }
                        // OCCT L683: // j++;
                    }
                }

                // OCCT L687-699.
                if self.bad_shape().is_null() {
                    // OCCT L689: BRep_Builder B; (architecture difference #7)
                    let mut brep = BRep::new();
                    // OCCT L690-691: B.MakeFace(TheNewFace, S2, Confusion());
                    let the_new_face = brep.add_tface_tol(
                        s2.clone(),
                        builder_make_wire(),
                        vec![],
                        None,
                        None,
                        vec![],
                        true,
                        CONFUSION,
                    );
                    // OCCT L692: Einf.Add(TheNewFace);
                    self.emap_mut().change_from_index(i).add(&the_new_face);
                    // OCCT L693: Draft_FaceInfo FI(S2, false);
                    let fi = DraftFaceInfo::new_with_geometry(s2.as_ref(), false);
                    // OCCT L694: myFMap.Add(TheNewFace, FI);
                    self.fmap_mut().add(&the_new_face, fi);
                } else {
                    // OCCT L698: break; — leave from for
                    break;
                }
            }
        }
        // OCCT L703: return (badShape.IsNull());
        self.bad_shape().is_null()
    }

    /// OCCT Draft_Modification::Perform() (Draft_Modification_1.cxx
    /// L708-1776).  The single OCCT member is carried as the entry plus the
    /// OCCT-commented sections (see the module header); the sections live in
    /// this module (Calculate eventual faces) and in
    /// offset/draft_modification_1_c.rs.
    pub fn perform(&mut self) {
        // OCCT L710-713.
        if !self.bad_shape().is_null() {
            panic!("Standard_ConstructionError");
        }

        // OCCT L715-721.
        if !self.comp() {
            self.set_comp(true);
            if !self.propagate() {
                return;
            }

            // OCCT L723-824: Calculate eventual faces.
            self.perform_calculate_eventual_faces();

            // OCCT L826-1443: Calculate new edges.
            self.perform_calculate_edges();

            // OCCT L1445-1636: Calculate new vertices.
            self.perform_calculate_vertices();
        }

        // OCCT L1639-1775: the small loop of validation/protection.
        self.perform_validation();
    }

    /// OCCT Perform L723-824 — "Calculate eventual faces" (the postponed
    /// circle faces become linear extrusions).
    fn perform_calculate_eventual_faces(&mut self) {
        for i in 1..=self.fmap().extent() {
            let fk = self.fmap().find_key(i).clone();
            let finf = self.fmap_mut().change_from_index(i);
            // OCCT L729.
            if finf.new_geometry() && finf.geometry().is_none() {
                let f1 = finf.first_face().clone();
                let f2 = finf.second_face().clone();
                // OCCT L734-739.
                if f1.is_null() || f2.is_null() {
                    self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                    self.bad_shape_mut().clone_from(&fk);
                    return;
                }
                // OCCT L740-747.
                let mut s1 = self.fmap().find_from_key(&f1).geometry().cloned();
                let mut s2 = self.fmap().find_from_key(&f2).geometry().cloned();
                if s1.is_none() || s2.is_none() {
                    self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                    self.bad_shape_mut().clone_from(&fk);
                    return;
                }
                // OCCT L748-755.
                if let Some(Surface3::Trimmed(t)) = &s1 {
                    s1 = Some(t.basis.as_ref().clone());
                }
                if let Some(Surface3::Trimmed(t)) = &s2 {
                    s2 = Some(t.basis.as_ref().clone());
                }
                // OCCT L756-763.
                let p1: Option<Plane> = match &s1 {
                    Some(Surface3::Plane(p)) => Some(*p),
                    _ => None,
                };
                let p2: Option<Plane> = match &s2 {
                    Some(Surface3::Plane(p)) => Some(*p),
                    _ => None,
                };
                if p1.is_none() || p2.is_none() {
                    self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                    self.bad_shape_mut().clone_from(&fk);
                    return;
                }
                // OCCT L764-772: IntAna_QuadQuadGeo i2p(pp1, pp2, ...).
                let i2p = super::draft_modification_1_b::IntAnaQuadQuadGeo::new_plane_plane(
                    &p1.unwrap(),
                    &p2.unwrap(),
                );
                if !i2p.is_done()
                    || i2p.type_inter() != super::draft_modification_1_b::IntAnaTypeInter::Line
                {
                    self.set_err_stat(DraftErrorStatus::FaceRecomputation);
                    self.bad_shape_mut().clone_from(&fk);
                    return;
                }

                // OCCT L774: gp_Dir extrdir = i2p.Line(1).Direction();
                let mut extrdir = i2p.line(1).direction;

                // OCCT L776-792: preserve the direction of the base face.
                let mut ref_surf = brep_tool_surface(&fk);
                if let Some(Surface3::Trimmed(t)) = &ref_surf {
                    ref_surf = Some(t.basis.as_ref().clone());
                }
                let mut dir_ref = DVec3::ZERO;
                if let Some(Surface3::Cylinder(cy)) = &ref_surf {
                    dir_ref = cy.axis;
                } else if let Some(Surface3::LinearExtrusion(le)) = &ref_surf {
                    dir_ref = le.direction;
                }
                // OCCT L794-797.
                if extrdir.dot(dir_ref) < 0.0 {
                    extrdir = -extrdir;
                }

                // OCCT L802: CCir = down_cast<Geom_Circle>(Finf.Curve());
                let ccir = match self.fmap().find_from_index(i).curve() {
                    Some(Curve3::Circle(c)) => *c,
                    _ => panic!("Geom_Circle cast of the postponed circle"),
                };
                // OCCT L803: NewS = new Geom_SurfaceOfLinearExtrusion(CCir, extrdir);
                let mut new_s = Surface3::LinearExtrusion(LinearExtrusionSurface {
                    profile: Box::new(Curve3::Circle(ccir)),
                    direction: extrdir,
                });

                // OCCT L805-819.
                let bounds = brep_tools_uv_bounds(&fk);
                let (umin, umax, mut vmin, mut vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                if !is_negative_infinite_value(vmin) && !is_positive_infinite_value(vmax) {
                    let deltav = 2.0 * (vmax - vmin);
                    vmin -= deltav;
                    vmax += deltav;
                } else {
                    // very temporary
                    vmax = 300.0;
                    vmin = -300.0;
                }
                let _ = (umin, umax);

                // OCCT L821-822.
                new_s = Surface3::Trimmed(TrimmedSurface::new(
                    new_s,
                    0.0,
                    1.9 * std::f64::consts::PI,
                    vmin,
                    vmax,
                ));
                *self
                    .fmap_mut()
                    .change_from_index(i)
                    .change_geometry() = Some(new_s);
            }
        }
    }

    /// OCCT Draft_Modification::NewSurface(S, Oris, Direction, Angle,
    /// NeutralPlane) (Draft_Modification_1.cxx L1780-1959) — the private
    /// geometry mapper (the Rust name distinguishes the public
    /// BRepTools_Modification override).
    pub(crate) fn new_surface_geom(
        &mut self,
        the_s: &Surface3,
        oris: Orientation,
        direction: DVec3,
        angle: f64,
        neutral_plane: &Plane,
    ) -> Option<Surface3> {
        let mut new_s: Option<Surface3> = None;

        match the_s {
            // OCCT L1790-1806: the plane case.
            Surface3::Plane(pl) => {
                let mut axe = rcad_kernel::math::gp::Ax1::new(DVec3::ZERO, DVec3::Z);
                let mut theta = 0.0f64;
                if find_rotation(pl, oris, direction, angle, neutral_plane, &mut axe, &mut theta) {
                    if theta.abs() > ANGULAR {
                        new_s = Some(Surface3::Plane(
                            super::draft_modification_1_b::geom_plane_rotated(pl, &axe, theta),
                        ));
                    } else {
                        new_s = Some(the_s.clone());
                    }
                }
            }
            // OCCT L1807-1878: the cylinder case.
            Surface3::Cylinder(cy) => {
                // OCCT L1809-1816.
                let mut testdir = direction.dot(neutral_plane.normal);
                if testdir.abs() <= 1.0 - ANGULAR {
                    return new_s;
                }
                // OCCT L1817-1825.
                testdir = direction.dot(cy.axis);
                if testdir.abs() <= 1.0 - ANGULAR {
                    return new_s;
                }
                // OCCT L1826-1873.
                if angle.abs() > ANGULAR {
                    let mut i2s = super::draft_modification_1_b::IntAnaQuadQuadGeo::default();
                    i2s.perform_plane_cylinder(neutral_plane, cy, ANGULAR, CONFUSION);
                    let mut is_int_done = i2s.is_done();

                    // OCCT L1832-1838: the ellipse-axis ratio guard.
                    if i2s.type_inter()
                        == super::draft_modification_1_b::IntAnaTypeInter::Ellipse
                    {
                        let an_el = i2s.ellipse(1);
                        let a_major_r = an_el.major_radius;
                        let a_minor_r = an_el.minor_radius;
                        is_int_done = a_major_r < 100000.0 * a_minor_r;
                    }

                    // OCCT L1840-1846.
                    if !is_int_done
                        || i2s.type_inter() != super::draft_modification_1_b::IntAnaTypeInter::Circle
                    {
                        return new_s;
                    }
                    // OCCT L1847-1872.
                    let mut alpha = angle;
                    // OCCT L1850: direct(axcone.Direct()) — the sense of the
                    // cylinder frame is the X/Y orientation about the axis.
                    let direct = cy
                        .ref_dir
                        .cross(cy.y_axis())
                        .dot(cy.axis)
                        >= 0.0;
                    if (direct && oris == Orientation::Reversed)
                        || (!direct && oris == Orientation::Forward)
                    {
                        alpha = -alpha;
                    }

                    let center = i2s.circle(1).center;
                    if testdir < 0.0 {
                        alpha = -alpha;
                    }
                    // OCCT L1861: Z = ElCLib::LineParameter(Cy.Axis(), Center);
                    let z = (center - cy.origin).dot(cy.axis);
                    let mut rad = cy.radius + z * alpha.tan();
                    if rad < 0.0 {
                        rad = -rad;
                    } else {
                        alpha = -alpha;
                    }
                    new_s = Some(Surface3::Cone(ConicalSurface {
                        apex: cone_apex_from(cy.origin, cy.axis, alpha, rad),
                        axis: cy.axis,
                        radius: rad,
                        half_angle_rad: alpha,
                        ref_dir: cy.ref_dir,
                    }));
                } else {
                    // OCCT L1876.
                    new_s = Some(the_s.clone());
                }
            }
            // OCCT L1879-1951: the cone case.
            Surface3::Cone(co1) => {
                // OCCT L1882-1889.
                let mut testdir = direction.dot(neutral_plane.normal);
                if testdir.abs() <= 1.0 - ANGULAR {
                    return new_s;
                }
                // OCCT L1891-1900.
                testdir = direction.dot(co1.axis);
                if testdir.abs() <= 1.0 - ANGULAR {
                    return new_s;
                }

                // OCCT L1902-1910.
                let mut i2s = super::draft_modification_1_b::IntAnaQuadQuadGeo::default();
                i2s.perform_plane_cone(neutral_plane, co1, ANGULAR, CONFUSION);
                if !i2s.is_done()
                    || i2s.type_inter() != super::draft_modification_1_b::IntAnaTypeInter::Circle
                {
                    return new_s;
                }
                // OCCT L1911-1918.
                let mut alpha = angle;
                let direct = co1
                    .ref_dir
                    .cross(co1.axis.cross(co1.ref_dir))
                    .dot(co1.axis)
                    >= 0.0;
                if (direct && oris == Orientation::Reversed)
                    || (!direct && oris == Orientation::Forward)
                {
                    alpha = -alpha;
                }

                let center = i2s.circle(1).center;
                // OCCT L1921-1950.
                if angle.abs() > ANGULAR {
                    if testdir < 0.0 {
                        alpha = -alpha;
                    }
                    // OCCT L1927: Z = ElCLib::LineParameter(Co1.Axis(), Center);
                    let z = (center - co1.apex).dot(co1.axis);
                    let mut rad = i2s.circle(1).radius + z * alpha.tan();
                    if rad < 0.0 {
                        rad = -rad;
                    } else {
                        alpha = -alpha;
                    }
                    if (alpha - co1.half_angle_rad).abs() < ANGULAR {
                        new_s = Some(the_s.clone());
                    } else {
                        new_s = Some(Surface3::Cone(ConicalSurface {
                            apex: cone_apex_from(co1.apex, co1.axis, alpha, rad),
                            axis: co1.axis,
                            radius: rad,
                            half_angle_rad: alpha,
                            ref_dir: co1.ref_dir,
                        }));
                    }
                } else {
                    // OCCT L1949: the degenerate cone becomes a cylinder.
                    new_s = Some(Surface3::Cylinder(CylindricalSurface {
                        origin: center,
                        axis: co1.axis,
                        radius: i2s.circle(1).radius,
                        ref_dir: co1.ref_dir,
                        y_dir: None,
                    }));
                }
            }
            // OCCT L1952-1957: the unsupported surfaces.
            _ => {}
        }
        // OCCT L1958: return NewS;
        new_s
    }

    /// OCCT Draft_Modification::NewCurve(C, S, Oris, Direction, Angle,
    /// NeutralPlane, Flag) (Draft_Modification_1.cxx L1963-2043) — the
    /// private geometry mapper (the Flag parameter is unnamed/unused in the
    /// OCCT text).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new_curve_geom(
        &mut self,
        the_c: &Curve3,
        the_s: &Surface3,
        oris: Orientation,
        direction: DVec3,
        angle: f64,
        neutral_plane: &Plane,
        _the_flag: bool,
    ) -> Option<Curve3> {
        let mut new_c: Option<Curve3> = None;

        // OCCT L1976-1993: the plane case rotates the curve.
        if let Surface3::Plane(pl) = the_s {
            let mut axe = rcad_kernel::math::gp::Ax1::new(DVec3::ZERO, DVec3::Z);
            let mut theta = 0.0f64;
            if find_rotation(pl, oris, direction, angle, neutral_plane, &mut axe, &mut theta) {
                if theta.abs() > ANGULAR {
                    new_c = Some(super::draft_modification_1_b::geom_curve_rotated(
                        the_c, &axe, theta,
                    ));
                } else {
                    new_c = Some(the_c.clone());
                }
            }
            // OCCT L1992: return NewC;
            return new_c;
        }

        // OCCT L1995-1998.
        if !matches!(the_c, Curve3::Line(_)) {
            return new_c;
        }

        // OCCT L2000: gp_Lin lin = down_cast<Geom_Line>(C)->Lin();
        let lin = match the_c {
            Curve3::Line(l) => *l,
            _ => unreachable!(),
        };
        // OCCT L2005-2025: Norm = the surface normal at the line location.
        let mut norm = DVec3::ZERO;
        match the_s {
            Surface3::Cylinder(cy) => {
                let (u, v) = elslib_parameters_cyl(cy, lin.origin);
                let (_, d1u, d1v) =
                    rcad_kernel::geom::Surface3::Cylinder(*cy).derivatives(u, v);
                norm = d1u.cross(d1v);
            }
            Surface3::Cone(co) => {
                let (u, v) = elslib_parameters_cone(co, lin.origin);
                let (_, d1u, d1v) = rcad_kernel::geom::Surface3::Cone(*co).derivatives(u, v);
                norm = d1u.cross(d1v);
            }
            _ => {}
        }

        // OCCT L2027-2042.
        let ilipl = IntAnaIntConicQuad::new(&lin, neutral_plane, ANGULAR);
        if ilipl.is_done() && ilipl.nb_points() != 0 {
            if oris == Orientation::Reversed {
                norm = -norm;
            }
            // OCCT L2034: gp_Ax1 axrot(ilipl.Point(1), Norm.Crossed(Direction));
            let axrot = rcad_kernel::math::gp::Ax1::new(ilipl.point(1), norm.cross(direction));
            // OCCT L2035: lires = gp_Lin(Ax1(Point, Direction)).Rotated(axrot, Angle).
            let base_lin = Line3::new(ilipl.point(1), direction);
            let mut lires = rotated_line(&base_lin, &axrot, angle);
            if lires.direction.dot(lin.direction) < 0.0 {
                lires.direction = -lires.direction;
            }
            new_c = Some(Curve3::Line(lires));
        }
        // OCCT L2042: return NewC;
        new_c
    }
}

// ===========================================================================
// Pure-math helpers of the _1 module.
// ===========================================================================

/// OCCT gp_Ax3 from a Geom_Conic::Position() — the conic axis plane carrier
/// of Propagate L577 (gp_Ax3(Ax1) constructor semantics).
fn conic_position_plane(location: DVec3, normal: DVec3) -> Plane {
    // gp_Ax3(Ax1) — the Ax3 on that axis with the gp_Ax2 default reference
    // direction.
    let a2 = rcad_kernel::math::gp::Ax2::from_direction(location, normal);
    let a3 = Ax3::from_ax2(&a2);
    ax3_plane(&a3)
}

/// The rcad Plane of an OCCT gp_Ax3 (location/direction/X/Y + sense).
fn ax3_plane(a3: &Ax3) -> Plane {
    Plane {
        origin: a3.location(),
        normal: a3.direction(),
        u_dir: a3.x_direction,
        v_dir: a3.y_direction,
    }
}

/// OCCT Geom_ElementarySurface::Axis() — the surface axis (Propagate
/// L583-591).
fn elementary_surface_axis(s: &Surface3) -> rcad_kernel::math::gp::Ax1 {
    match s {
        Surface3::Cylinder(c) => rcad_kernel::math::gp::Ax1::new(c.origin, c.axis),
        Surface3::Cone(c) => rcad_kernel::math::gp::Ax1::new(c.apex, c.axis),
        Surface3::Sphere(s) => rcad_kernel::math::gp::Ax1::new(s.center, s.axis),
        Surface3::Trimmed(t) => elementary_surface_axis(t.basis.as_ref()),
        _ => panic!("Geom_ElementarySurface::Axis: unsupported surface in the Draft flows"),
    }
}

/// OCCT ElSLib::Parameters(Pl, P, U, V) (ElSLib.cxx) — the plane UV.
pub(crate) fn elslib_parameters_plane(pl: &Plane, p: DVec3) -> (f64, f64) {
    let d = p - pl.origin;
    (d.dot(pl.u_dir), d.dot(pl.v_dir))
}

/// OCCT ElSLib::Parameters(Cy, P, U, V) — the cylinder UV (u measured on the
/// X/Y frame of the rcad cylinder).
pub(crate) fn elslib_parameters_cyl(cy: &CylindricalSurface, p: DVec3) -> (f64, f64) {
    let d = p - cy.origin;
    let v = d.dot(cy.axis);
    let radial = d - cy.axis * v;
    let y_dir = cy.y_axis();
    let u = radial.dot(y_dir).atan2(radial.dot(cy.ref_dir));
    (u, v)
}

/// OCCT ElSLib::Parameters(Co, P, U, V) — the cone UV.
pub(crate) fn elslib_parameters_cone(co: &ConicalSurface, p: DVec3) -> (f64, f64) {
    let d = p - co.apex;
    let v = d.dot(co.axis);
    let radial = d - co.axis * v;
    let y_dir = co.axis.cross(co.ref_dir);
    let u = radial.dot(y_dir).atan2(radial.dot(co.ref_dir));
    (u, v)
}

/// OCCT gp_Cone apex — the apex point derived from the reference circle
/// (radius at `origin`, half angle); Vapex = -RefRadius / sin(SemiAngle)
/// from the reference circle along the axis.
fn cone_apex_from(origin: DVec3, axis: DVec3, half_angle: f64, radius: f64) -> DVec3 {
    origin - axis * (radius / half_angle.sin())
}

/// OCCT gp_Lin::Rotated(A1, Ang) — the line rotation through the rotation
/// gp_Trsf (the gp_Trsf::SetRotation + gp_Lin::Transform vehicle).
fn rotated_line(l: &Line3, axe: &rcad_kernel::math::gp::Ax1, ang: f64) -> Line3 {
    let t = super::draft_modification_1_b::gp_trsf_rotation(axe, ang);
    Line3 {
        origin: t.apply(l.origin),
        direction: t.transform_dir(l.direction),
    }
}
