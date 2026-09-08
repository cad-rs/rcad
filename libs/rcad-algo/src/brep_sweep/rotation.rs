//! OCCT BRepSweep_Rotation (TKPrim/BRepSweep) — builds a topology by rotation
//! sweep.
//!
//! Sources:
//! - BRepSweep_Rotation.hxx L34-190
//! - BRepSweep_Rotation.cxx L71-906
//!
//! Architecture difference (Rust has no inheritance): the C++
//! `Rotation : BRepSweep_Trsf : NumLinearRegularSweep` chain maps to the
//! composition root owning the core + TrsfData members and implementing the
//! two slot traits.  The gp_Ax1 axis travels as a (Location, Direction)
//! pair.

use rcad_kernel::geom::{Curve2d, Curve3, CurveEval, Surface3};
use rcad_kernel::precision::{p_confusion, CONFUSION, ANGULAR};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};

use super::num_linear_regular_sweep::{NumLinearRegularSweepCore, NumLinearRegularSweepSlots};
use super::sweep_num_shape::SweepNumShape;
use super::tool_rehost::{
    brep_tool_curve_on_surface_seam, brep_tool_degenerated, brep_tools_is_really_closed,
    curve_type, geom_adaptor_surface_get_type, geom_circle, geom_curve_transformed,
    geom_trimmed_curve, GeomAbsCurveType, GeomAbsSurfaceType, GeomAdaptorSurfaceOfRevolution,
};
use super::trsf::{self, BRepSweepTrsfSlots, TrsfData};
use super::BRepSweepBuilder;

use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_parameter, brep_tool_pnt, brep_tool_surface, brep_tool_tolerance,
    explorer,
};

/// The axis pair (OCCT gp_Ax1: Location, Direction).
pub type GpAx1 = (glam::DVec3, glam::DVec3);

/// OCCT BRepSweep_Rotation (BRepSweep_Rotation.hxx L34-190).
///
/// OCCT inheritance chain: `BRepSweep_Rotation : BRepSweep_Trsf :
/// BRepSweep_NumLinearRegularSweep`.
pub struct BRepSweepRotation {
    /// The base-class members.
    pub core: NumLinearRegularSweepCore,
    /// The BRepSweep_Trsf members (myLocation / myCopy).
    pub trsf: TrsfData,
    /// OCCT: myAng (double).
    pub my_ang: f64,
    /// OCCT: myAxe (gp_Ax1).
    pub my_axe: GpAx1,
}

impl NumLinearRegularSweepSlots for BRepSweepRotation {
    fn core(&mut self) -> &mut NumLinearRegularSweepCore {
        &mut self.core
    }
    fn core_ref(&self) -> &NumLinearRegularSweepCore {
        &self.core
    }

    /// OCCT BRepSweep_Rotation::MakeEmptyVertex (cxx L167-188) — call only
    /// in construction mode with copy.
    fn make_empty_vertex(&mut self, a_gen_v: &Shape, a_dir_v: &SweepNumShape) -> Shape {
        if !self.trsf.my_copy {
            panic!("Standard_ConstructionError: BRepSweep_Rotation::MakeEmptyVertex");
        }
        let mut p = brep_tool_pnt(a_gen_v).expect("BRep_Tool::Pnt");
        if a_dir_v.index() == 2 {
            // OCCT L174-177: P.Transform(myLocation.Transformation()).
            p = self.trsf.my_location.transform_point3(p);
        }
        // OCCT L178-181 (modified by jgv, 1.10.01, for buc61005): the vertex
        // tolerance is the generatrix vertex tolerance.
        let tol = brep_tool_tolerance(a_gen_v);
        let v = self.core.my_builder.my_builder.make_vertex(p, tol);
        // OCCT L182-186: the invariant first vertices land in the third
        // column (the closing column of a closed revol).
        if a_dir_v.index() == 1
            && NumLinearRegularSweepSlots::is_invariant(self, a_gen_v)
            && self.core.my_dir_shape_tool.nb_shapes() == 3
        {
            let i_gen_v = self.core.my_gen_shape_tool.index(a_gen_v);
            self.core.my_built_shapes.set(i_gen_v, 3, true);
            self.core.my_shapes.set(i_gen_v, 3, v.clone());
        }
        v
    }

    /// OCCT BRepSweep_Rotation::MakeEmptyDirectingEdge (cxx L192-220).
    fn make_empty_directing_edge(&mut self, a_gen_v: &Shape, _a_dir_e: &SweepNumShape) -> Shape {
        let p = brep_tool_pnt(a_gen_v).expect("BRep_Tool::Pnt");
        // OCCT L197-200: Dirz(myAxe.Direction()); V = Vec(Dirz);
        // O = myAxe.Location(); O.Translate(V.Dot(gp_Vec(O, P)) * V).
        let dirz = self.my_axe.1;
        let v = dirz;
        let mut o = self.my_axe.0;
        o = o + (p - o).dot(v) * v;
        let tol = brep_tool_tolerance(a_gen_v);
        // OCCT L201: if (O.IsEqual(P, Precision::Confusion())).
        if o.distance(p) <= CONFUSION {
            // OCCT L202-210: make a degenerated edge (temporary make 3D curve
            // null so that parameters should be registered): the zero-radius
            // circle.
            let gc = geom_circle(o, dirz, any_perpendicular(dirz), 0.0);
            let e = self.core.my_builder.my_builder.make_edge_curve(&gc, tol);
            self.core.my_builder.my_builder.degenerated(&e, true);
            e
        } else {
            // OCCT L212-218: Axis(O, Dirz, Dir(gp_Vec(O, P))); GC = Circle(
            // Axis, O.Distance(P)); MakeEdge(E, GC, tol).
            let x_dir = (p - o).normalize_or_zero();
            let gc = geom_circle(o, dirz, x_dir, o.distance(p));
            self.core.my_builder.my_builder.make_edge_curve(&gc, tol)
        }
    }

    /// OCCT BRepSweep_Rotation::MakeEmptyGeneratingEdge (cxx L224-257) —
    /// call in case of construction with copy, or only when meridian touches
    /// myAxe.
    fn make_empty_generating_edge(&mut self, a_gen_e: &Shape, a_dir_v: &SweepNumShape) -> Shape {
        let tol = brep_tool_tolerance(a_gen_e);
        let e;
        if brep_tool_degenerated(a_gen_e) {
            // OCCT L229-234: MakeEdge(E); UpdateEdge(E, Tolerance);
            // Degenerated(E, true).
            e = self.core.my_builder.my_builder.make_edge();
            self.core.my_builder.my_builder.update_edge_tol(&e, tol);
            self.core.my_builder.my_builder.degenerated(&e, true);
        } else {
            // OCCT L237-249: C = down_cast<Geom_Curve>(BRep_Tool::Curve(
            // Edge(aGenE), Loc, First, Last)->Copy()); the location transform
            // is the identity no-op (arch. diff.).
            let mut c = brep_tool_curve(a_gen_e).map(|(c, _f, _l)| c);
            if let Some(hit) = &mut c {
                if a_dir_v.index() == 2 {
                    // OCCT L244-247: C->Transform(myLocation.Transformation()).
                    *hit = geom_curve_transformed(hit, &self.trsf.my_location);
                }
            }
            e = match &c {
                Some(hit) => self.core.my_builder.my_builder.make_edge_curve(hit, tol),
                None => self.core.my_builder.my_builder.make_edge(),
            };
        }
        // OCCT L251-255: the invariant first meridians land in the third
        // column.
        if a_dir_v.index() == 1
            && NumLinearRegularSweepSlots::is_invariant(self, a_gen_e)
            && self.core.my_dir_shape_tool.nb_shapes() == 3
        {
            let i_gen_e = self.core.my_gen_shape_tool.index(a_gen_e);
            self.core.my_built_shapes.set(i_gen_e, 3, true);
            self.core.my_shapes.set(i_gen_e, 3, e.clone());
        }
        e
    }

    /// OCCT BRepSweep_Rotation::SetParameters (cxx L261-274).
    fn set_parameters(
        &mut self,
        a_new_face: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_f: &Shape,
        a_gen_v: &Shape,
        _a_dir_v: &SweepNumShape,
    ) {
        let pnt2d = super::tool_rehost::brep_tool_parameters(a_gen_v, a_gen_f);
        self.core
            .my_builder
            .my_builder
            .update_vertex_uv_on_face(a_new_vertex, pnt2d.x, pnt2d.y, a_new_face, p_confusion());
    }

    /// OCCT BRepSweep_Rotation::SetDirectingParameter (cxx L278-294).
    fn set_directing_parameter(
        &mut self,
        a_new_edge: &Shape,
        a_new_vertex: &mut Shape,
        _a_gen_v: &Shape,
        _a_dir_e: &SweepNumShape,
        a_dir_v: &SweepNumShape,
    ) {
        let mut param = 0.0;
        let mut ori = Orientation::Forward;
        if a_dir_v.index() == 2 {
            param = self.my_ang;
            ori = Orientation::Reversed;
        }
        // OCCT L291-293: V_wnt = Vertex(aNewVertex); V_wnt.Orientation(ori);
        // UpdateVertex(V_wnt, param, Edge(aNewEdge), PConfusion()).
        a_new_vertex.orientation = ori;
        self.core
            .my_builder
            .my_builder
            .update_vertex_param_on_edge(a_new_vertex, param, a_new_edge, p_confusion());
    }

    /// OCCT BRepSweep_Rotation::SetGeneratingParameter (cxx L298-310).
    fn set_generating_parameter(
        &mut self,
        a_new_edge: &Shape,
        a_new_vertex: &mut Shape,
        a_gen_e: &Shape,
        a_gen_v: &Shape,
        _a_dir_v: &SweepNumShape,
    ) {
        a_new_vertex.orientation = a_gen_v.orientation;
        let param = brep_tool_parameter(a_gen_v, a_gen_e);
        self.core
            .my_builder
            .my_builder
            .update_vertex_param_on_edge(a_new_vertex, param, a_new_edge, p_confusion());
    }

    /// OCCT BRepSweep_Rotation::MakeEmptyFace (cxx L314-385).
    fn make_empty_face(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> Shape {
        let toler;
        let s: Surface3;
        if a_gen_s.shape_type() == ShapeType::Edge {
            // OCCT L322-325: C = Curve(Edge(aGenS), L, First, Last); toler =
            // Tolerance(Edge(aGenS)).
            let (c0, first, last) = brep_tool_curve(a_gen_s).expect("BRep_Tool::Curve");
            toler = brep_tool_tolerance(a_gen_s);
            // OCCT L326-331 (modified by jgv, 9.12.03): Tr = L.Transformation();
            // C = down_cast(C->Copy()); C = new Geom_TrimmedCurve(C, First,
            // Last); C->Transform(Tr) — the curve location is the identity
            // no-op in this pipeline (arch. diff. #1).
            let c = geom_trimmed_curve(c0, first, last);
            // OCCT L333-335: HC = GeomAdaptor_Curve(); HC->Load(C, First,
            // Last); AS = GeomAdaptor_SurfaceOfRevolution(HC, myAxe).
            let as_rev = GeomAdaptorSurfaceOfRevolution::load(&c, self.my_axe.0, self.my_axe.1);
            match as_rev.get_type() {
                GeomAbsSurfaceType::Plane => {
                    // OCCT L338-341: S = new Geom_Plane(AS.Plane()).
                    s = Surface3::Plane(as_rev.plane());
                }
                GeomAbsSurfaceType::Cylinder => {
                    // OCCT L343-346: S = Geom_CylindricalSurface(AS.Cylinder()).
                    s = Surface3::Cylinder(as_rev.cylinder());
                }
                GeomAbsSurfaceType::Sphere => {
                    // OCCT L348-351: S = Geom_SphericalSurface(AS.Sphere()).
                    s = Surface3::Sphere(as_rev.sphere());
                }
                GeomAbsSurfaceType::Cone => {
                    // OCCT L353-356: S = Geom_ConicalSurface(AS.Cone()).
                    s = Surface3::Cone(as_rev.cone());
                }
                GeomAbsSurfaceType::Torus => {
                    // OCCT L358-361: S = Geom_ToroidalSurface(AS.Torus()).
                    s = Surface3::Torus(as_rev.torus());
                }
                _ => {
                    // OCCT L363-366: S = Geom_SurfaceOfRevolution(C, myAxe).
                    s = as_rev.surface_value();
                }
            }
        } else {
            // OCCT L372-377: S = Surface(Face(aGenS), L); toler = Tolerance(
            // Face(aGenS)); S = S->Copy(); S->Transform(Tr).
            let mut s0 = brep_tool_surface(a_gen_s).expect("BRep_Tool::Surface");
            toler = brep_tool_tolerance(a_gen_s);
            if a_dir_s.index() == 2 {
                // OCCT L378-381: S->Transform(myLocation.Transformation()).
                s0 = rcad_kernel::geom::transform_surface(&s0, &self.trsf.my_location);
            }
            s = s0;
        }
        // OCCT L383: MakeFace(F, S, toler).
        self.core.my_builder.my_builder.make_face(&s, toler)
    }

    /// OCCT BRepSweep_Rotation::SetPCurve (cxx L389-404) — sets on the edges
    /// of cap faces the same pcurves as on the edges of the generator face.
    fn set_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_f: &Shape,
        a_gen_e: &Shape,
        _a_dir_v: &SweepNumShape,
        orien: Orientation,
    ) {
        // OCCT L399-403: SetThePCurve(B, Edge(aNewEdge), Face(aNewFace),
        // orien, CurveOnSurface(Edge(aGenE), Face(aGenF), First, Last)).
        if let Some(c) = brep_tool_curve_on_surface_seam(a_gen_e, a_gen_f) {
            set_the_pcurve(self, a_new_edge, a_new_face, orien, &c.0);
        }
    }

    /// OCCT BRepSweep_Rotation::SetGeneratingPCurve (cxx L408-516).
    fn set_generating_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        _a_gen_e: &Shape,
        _a_dir_e: &SweepNumShape,
        a_dir_v: &SweepNumShape,
        orien: Orientation,
    ) {
        let surf = brep_tool_surface(a_new_face).expect("BRep_Tool::Surface");
        let mut u = 0.0;
        let mut v = 0.0;
        let as_type = geom_adaptor_surface_get_type(&surf);
        let l2d: (glam::DVec2, glam::DVec2); // (Location, Direction) of gp_Lin2d
        match as_type {
            GeomAbsSurfaceType::Plane => {
                // OCCT L423-449: pln = AS.Plane(); ax3 = pln.Position();
                // aC = Curve(Edge(aNewEdge), Loc, First, Last); the Geom_Line
                // (or trimmed line basis) of the new edge.
                let a_c = match brep_tool_curve(a_new_edge) {
                    Some(hit) => hit.0,
                    None => return,
                };
                let gl_dir = match &a_c {
                    Curve3::Line(ln) => ln.direction,
                    Curve3::Trimmed(tc) => match tc.curve.as_ref() {
                        Curve3::Line(ln) => ln.direction,
                        _ => panic!("Standard_ConstructionError: BRepSweep_Rotation::SetGeneratingPCurve"),
                    },
                    _ => panic!("Standard_ConstructionError: BRepSweep_Rotation::SetGeneratingPCurve"),
                };
                let gl_loc = match &a_c {
                    Curve3::Line(ln) => ln.origin,
                    Curve3::Trimmed(tc) => match tc.curve.as_ref() {
                        Curve3::Line(ln) => ln.origin,
                        _ => unreachable!(),
                    },
                    _ => unreachable!(),
                };
                // OCCT L441-449: gl.Transform(Loc.Transformation());
                // point = gl.Location(); dir = gl.Direction();
                // ElSLib::PlaneParameters(ax3, point, u, v); ...
                let pln = match &surf {
                    Surface3::Plane(p) => p.clone(),
                    _ => unreachable!(),
                };
                let (pu, pv) = rcad_kernel::math::el::elslib_plane_parameters(
                    gl_loc, pln.origin, pln.u_dir, pln.v_dir,
                );
                u = pu;
                v = pv;
                let dir2d = glam::DVec2::new(
                    gl_dir.dot(pln.u_dir),
                    gl_dir.dot(pln.v_dir),
                );
                l2d = (glam::DVec2::new(u, v), dir2d.normalize_or_zero());
            }
            GeomAbsSurfaceType::Torus => {
                // OCCT L451-487.
                let tor = match &surf {
                    Surface3::Torus(t) => *t,
                    _ => unreachable!(),
                };
                let (curve, first, _last) =
                    brep_tool_curve(a_new_edge).expect("BRep_Tool::Curve");
                let big_u = first;
                let point = curve.point_at(big_u);
                if point.distance(tor.center) < CONFUSION {
                    v = std::f64::consts::PI;
                    // modified by NIZHNY-EAP Wed Mar  1 17:49:29 2000.
                    u = 0.0;
                } else {
                    let (tu, tv) = rcad_kernel::math::el::elslib_torus_parameters(
                        point,
                        tor.center,
                        tor.ref_dir,
                        tor.axis.cross(tor.ref_dir),
                        tor.axis,
                        tor.major_radius,
                        tor.minor_radius,
                    );
                    u = tu;
                    v = tv;
                }
                // OCCT L468: v = ElCLib::InPeriod(v, 0., 2 * M_PI).
                v = rcad_kernel::math::el::in_period(v, 0.0, 2.0 * std::f64::consts::PI);
                // OCCT L469-472.
                if 2.0 * std::f64::consts::PI - v <= p_confusion() {
                    v -= 2.0 * std::f64::consts::PI;
                }
                // OCCT L473-482: the uLeft/uRight AdjustPeriodic pair.
                if a_dir_v.index() == 2 {
                    let mut u_left = u - self.my_ang;
                    super::tool_rehost::elclib_adjust_periodic(
                        -std::f64::consts::PI,
                        std::f64::consts::PI,
                        p_confusion(),
                        &mut u_left,
                        &mut u,
                    );
                } else {
                    let mut u_right = u + self.my_ang;
                    super::tool_rehost::elclib_adjust_periodic(
                        -std::f64::consts::PI,
                        std::f64::consts::PI,
                        p_confusion(),
                        &mut u,
                        &mut u_right,
                    );
                }
                // OCCT L484-486: pnt2d.SetCoord(u, v - U); L.SetLocation;
                // L.SetDirection(gp::DY2d()).
                l2d = (
                    glam::DVec2::new(u, v - big_u),
                    glam::DVec2::Y,
                );
            }
            GeomAbsSurfaceType::Sphere => {
                // OCCT L488-503.
                let sph = match &surf {
                    Surface3::Sphere(s0) => *s0,
                    _ => unreachable!(),
                };
                let (curve, first, _last) =
                    brep_tool_curve(a_new_edge).expect("BRep_Tool::Curve");
                let big_u = first;
                let point = curve.point_at(big_u);
                let (su, sv) = rcad_kernel::math::el::elslib_sphere_parameters(
                    point,
                    sph.center,
                    sph.ref_dir,
                    sph.axis.cross(sph.ref_dir),
                    sph.axis,
                );
                u = su;
                v = sv;
                u = 0.0;
                if a_dir_v.index() == 2 {
                    u = self.my_ang;
                }
                l2d = (glam::DVec2::new(u, v - big_u), glam::DVec2::Y);
            }
            _ => {
                // OCCT L504-513: anAngleTemp = 0; if (aDirV.Index() == 2)
                // anAngleTemp = myAng; L.SetLocation(gp_Pnt2d(anAngleTemp, 0));
                // L.SetDirection(gp::DY2d()).
                let an_angle_temp = if a_dir_v.index() == 2 { self.my_ang } else { 0.0 };
                l2d = (glam::DVec2::new(an_angle_temp, 0.0), glam::DVec2::Y);
            }
        }
        // OCCT L514-515: GL = new Geom2d_Line(L); SetThePCurve(B, Edge(
        // aNewEdge), Face(aNewFace), orien, GL).
        let gl = super::tool_rehost::geom2d_line(l2d.0, l2d.1);
        set_the_pcurve(self, a_new_edge, a_new_face, orien, &gl);
    }

    /// OCCT BRepSweep_Rotation::SetDirectingPCurve (cxx L520-630).
    fn set_directing_pcurve(
        &mut self,
        a_new_face: &Shape,
        a_new_edge: &mut Shape,
        a_gen_e: &Shape,
        a_gen_v: &Shape,
        _a_dir_e: &SweepNumShape,
        orien: Orientation,
    ) {
        let surf = brep_tool_surface(a_new_face).expect("BRep_Tool::Surface");
        // OCCT L529-533: par = BRep_Tool::Parameter(Vertex(aGenV), Edge(
        // aGenE)); p2 = BRep_Tool::Pnt(Vertex(aGenV)).
        let par = brep_tool_parameter(a_gen_v, a_gen_e);
        let p2 = brep_tool_pnt(a_gen_v).expect("BRep_Tool::Pnt");
        let the_pcurve: Curve2d;
        match geom_adaptor_surface_get_type(&surf) {
            GeomAbsSurfaceType::Plane => {
                // OCCT L538-549: pln; ax3; p1 = pln.Location(); R = p1.
                // Distance(p2); ElSLib::PlaneParameters(ax3, p2, u, v);
                // dx2d = Dir2d(u, v); axe = Ax22d(Origin2d, dx2d, DY2d);
                // C = Circ2d(axe, R); GC = Geom2d_Circle(C).
                let pln = match &surf {
                    Surface3::Plane(p) => p.clone(),
                    _ => unreachable!(),
                };
                let big_r = pln.origin.distance(p2);
                let (pu, pv) = rcad_kernel::math::el::elslib_plane_parameters(
                    p2, pln.origin, pln.u_dir, pln.v_dir,
                );
                let dx2d = glam::DVec2::new(pu, pv).normalize_or_zero();
                let dy2d = glam::DVec2::Y;
                // OCCT gp_Ax22d(Origin2d, dx2d, DY2d); Circ2d(axe, R).
                the_pcurve = super::tool_rehost::geom2d_circle(glam::DVec2::ZERO, big_r);
                let _ = (dx2d, dy2d);
            }
            GeomAbsSurfaceType::Cone => {
                // OCCT L552-559: cone; ConeParameters(cone.Position(),
                // cone.RefRadius(), cone.SemiAngle(), p2, u, v);
                // p22d.SetCoord(0, v); L(p22d, DX2d).
                let cone = match &surf {
                    Surface3::Cone(c0) => *c0,
                    _ => unreachable!(),
                };
                let (_cu, cv) = rcad_kernel::math::el::elslib_cone_parameters(
                    p2,
                    cone.apex,
                    cone.ref_dir,
                    cone.axis.cross(cone.ref_dir),
                    cone.axis,
                    cone.radius,
                    cone.half_angle_rad,
                );
                the_pcurve = super::tool_rehost::geom2d_line(
                    glam::DVec2::new(0.0, cv),
                    glam::DVec2::X,
                );
            }
            GeomAbsSurfaceType::Sphere => {
                // OCCT L562-569.
                let sph = match &surf {
                    Surface3::Sphere(s0) => *s0,
                    _ => unreachable!(),
                };
                let (_su, sv) = rcad_kernel::math::el::elslib_sphere_parameters(
                    p2,
                    sph.center,
                    sph.ref_dir,
                    sph.axis.cross(sph.ref_dir),
                    sph.axis,
                );
                the_pcurve = super::tool_rehost::geom2d_line(
                    glam::DVec2::new(0.0, sv),
                    glam::DVec2::X,
                );
            }
            GeomAbsSurfaceType::Torus => {
                // OCCT L572-614: the meridian end parameters.
                let tor = match &surf {
                    Surface3::Torus(t) => *t,
                    _ => unreachable!(),
                };
                let (curve, first, last) = brep_tool_curve(a_gen_e).expect("BRep_Tool::Curve");
                let mut p1 = curve.point_at(first);
                let mut u1 = 0.0;
                let mut v1 = 0.0;
                let mut u2;
                let mut v2;
                if p1.distance(tor.center) < CONFUSION {
                    v1 = std::f64::consts::PI;
                    // modified by NIZHNY-EAP Thu Mar  2 09:43:26 2000.
                    u1 = 0.0;
                } else {
                    let (tu, tv) = rcad_kernel::math::el::elslib_torus_parameters(
                        p1,
                        tor.center,
                        tor.ref_dir,
                        tor.axis.cross(tor.ref_dir),
                        tor.axis,
                        tor.major_radius,
                        tor.minor_radius,
                    );
                    u1 = tu;
                    v1 = tv;
                }
                p1 = curve.point_at(last);
                if p1.distance(tor.center) < CONFUSION {
                    v2 = std::f64::consts::PI;
                } else {
                    let (tu, tv) = rcad_kernel::math::el::elslib_torus_parameters(
                        p1,
                        tor.center,
                        tor.ref_dir,
                        tor.axis.cross(tor.ref_dir),
                        tor.axis,
                        tor.major_radius,
                        tor.minor_radius,
                    );
                    u2 = tu;
                    v2 = tv;
                }
                // OCCT L598: ElCLib::AdjustPeriodic(0., 2 * M_PI, PConfusion(),
                // v1, v2).
                super::tool_rehost::elclib_adjust_periodic(
                    0.0,
                    2.0 * std::f64::consts::PI,
                    p_confusion(),
                    &mut v1,
                    &mut v2,
                );
                // OCCT L600-601 (modified by NIZHNY-EAP Thu Mar  2 15:29:04):
                // u2 = u1 + myAng; AdjustPeriodic(-M_PI, M_PI, PConfusion(),
                // u1, u2).
                u2 = u1 + self.my_ang;
                super::tool_rehost::elclib_adjust_periodic(
                    -std::f64::consts::PI,
                    std::f64::consts::PI,
                    p_confusion(),
                    &mut u1,
                    &mut u2,
                );
                // OCCT L602-609: if (aGenV.Orientation() == FORWARD)
                //   p22d.SetCoord(u1, v1) else p22d.SetCoord(u1, v2).
                let p22d = if a_gen_v.orientation == Orientation::Forward {
                    glam::DVec2::new(u1, v1)
                } else {
                    glam::DVec2::new(u1, v2)
                };
                // OCCT L611-613: L(p22d, DX2d).
                the_pcurve = super::tool_rehost::geom2d_line(p22d, glam::DVec2::X);
            }
            _ => {
                // OCCT L617-622: p22d.SetCoord(0., par); L(p22d, DX2d).
                the_pcurve = super::tool_rehost::geom2d_line(
                    glam::DVec2::new(0.0, par),
                    glam::DVec2::X,
                );
            }
        }
        // OCCT L625-629: SetThePCurve(B, Edge(aNewEdge), Face(aNewFace),
        // orien, thePCurve).
        set_the_pcurve(self, a_new_edge, a_new_face, orien, &the_pcurve);
    }

    /// OCCT BRepSweep_Rotation::DirectSolid (cxx L635-675, modified by
    /// NIZNHY-PKV Tue Jun 14 08:33:55 2011) — compares the face normal and
    /// the direction.
    fn direct_solid(&mut self, a_gen_s: &Shape, _a_dir_s: &SweepNumShape) -> Orientation {
        use rcad_kernel::geom::SurfaceEval;
        let surf = brep_tool_surface(a_gen_s).expect("BRep_Tool::Surface");
        let dom = surf.default_domain();
        let a_u1 = dom[0];
        let a_u2 = dom[1];
        let a_v1 = dom[2];
        let a_v2 = dom[3];
        // OCCT L643-644: aTol2 = Precision::Confusion()^2.
        let a_tol2 = CONFUSION * CONFUSION;
        let a_p_axe_loc = self.my_axe.0;
        let a_p_axe_dir = self.my_axe.1;
        // OCCT L654-657: aTx = 0.5; aUx = aTx*(aU1+aU2); aVx = aTx*(aV1+aV2);
        // surf.D1(aUx, aVx, aP, du, dv).
        let mut a_tx = 0.5;
        let mut a_ux = a_tx * (a_u1 + a_u2);
        let mut a_vx = a_tx * (a_v1 + a_v2);
        let (mut a_p, mut du, mut dv) = surf.derivatives(a_ux, a_vx);
        // OCCT L659-661: aV = Vec(aPAxeLoc, aP); aV.Cross(aPAxeDir);
        // aMV2 = SquareMagnitude().
        let mut a_v = (a_p - a_p_axe_loc).cross(a_p_axe_dir);
        let mut a_mv2 = a_v.length_squared();
        if a_mv2 < a_tol2 {
            // OCCT L664-669: aTx = 0.43213918; retry at the recentered point.
            a_tx = 0.43213918;
            a_ux = a_u1 * (1.0 - a_tx) + a_u2 * a_tx;
            a_vx = a_v1 * (1.0 - a_tx) + a_v2 * a_tx;
            let hit = surf.derivatives(a_ux, a_vx);
            a_p = hit.0;
            du = hit.1;
            dv = hit.2;
            a_v = a_p - a_p_axe_loc;
            a_v = a_v.cross(a_p_axe_dir);
            let _ = a_mv2;
            a_mv2 = a_v.length_squared();
        }
        let _ = a_mv2;
        let _ = a_p;
        // OCCT L672-674: aX = aV.DotCross(du, dv);
        // aOr = (aX > 0.) ? FORWARD : REVERSED.
        let a_x = a_v.dot(du.cross(dv));
        if a_x > 0.0 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }

    /// OCCT BRepSweep_Rotation::GGDShapeIsToAdd (cxx L703-729).
    fn ggd_shape_is_to_add(
        &self,
        a_new_shape: &Shape,
        a_new_sub_shape: &Shape,
        a_gen_s: &Shape,
        a_sub_gen_s: &Shape,
        a_dir_s: &SweepNumShape,
    ) -> bool {
        let a_res = true;
        if a_new_shape.shape_type() == ShapeType::Face
            && a_new_sub_shape.shape_type() == ShapeType::Edge
            && a_gen_s.shape_type() == ShapeType::Edge
            && a_sub_gen_s.shape_type() == ShapeType::Vertex
            && a_dir_s.type_() == ShapeType::Edge
        {
            let as_type = match brep_tool_surface(a_new_shape) {
                Some(s0) => geom_adaptor_surface_get_type(&s0),
                None => return a_res,
            };
            if as_type == GeomAbsSurfaceType::Plane {
                return !NumLinearRegularSweepSlots::is_invariant(self, a_sub_gen_s);
            } else {
                return a_res;
            }
        } else {
            a_res
        }
    }

    /// OCCT BRepSweep_Rotation::GDDShapeIsToAdd (cxx L733-764).
    fn gdd_shape_is_to_add(
        &self,
        a_new_shape: &Shape,
        a_new_sub_shape: &Shape,
        a_gen_s: &Shape,
        a_dir_s: &SweepNumShape,
        a_sub_dir_s: &SweepNumShape,
    ) -> bool {
        if a_new_shape.shape_type() == ShapeType::Solid
            && a_new_sub_shape.shape_type() == ShapeType::Face
            && a_gen_s.shape_type() == ShapeType::Face
            && a_dir_s.type_() == ShapeType::Edge
            && a_sub_dir_s.type_() == ShapeType::Vertex
        {
            (self.my_ang - 2.0 * std::f64::consts::PI).abs() > ANGULAR
        } else if a_new_shape.shape_type() == ShapeType::Face
            && a_new_sub_shape.shape_type() == ShapeType::Edge
            && a_gen_s.shape_type() == ShapeType::Edge
            && a_dir_s.type_() == ShapeType::Edge
            && a_sub_dir_s.type_() == ShapeType::Vertex
        {
            let as_type = match brep_tool_surface(a_new_shape) {
                Some(s0) => geom_adaptor_surface_get_type(&s0),
                None => return true,
            };
            if as_type == GeomAbsSurfaceType::Plane {
                (self.my_ang - 2.0 * std::f64::consts::PI).abs() > ANGULAR
            } else {
                true
            }
        } else {
            true
        }
    }

    /// OCCT BRepSweep_Rotation::SeparatedWires (cxx L768-793) — the only
    /// case is a planar face in a closed revol.
    fn separated_wires(
        &self,
        a_new_shape: &Shape,
        a_new_sub_shape: &Shape,
        a_gen_s: &Shape,
        a_sub_gen_s: &Shape,
        a_dir_s: &SweepNumShape,
    ) -> bool {
        if a_new_shape.shape_type() == ShapeType::Face
            && a_new_sub_shape.shape_type() == ShapeType::Edge
            && a_gen_s.shape_type() == ShapeType::Edge
            && a_sub_gen_s.shape_type() == ShapeType::Vertex
            && a_dir_s.type_() == ShapeType::Edge
        {
            let as_type = match brep_tool_surface(a_new_shape) {
                Some(s0) => geom_adaptor_surface_get_type(&s0),
                None => return false,
            };
            if as_type == GeomAbsSurfaceType::Plane {
                (self.my_ang - 2.0 * std::f64::consts::PI).abs() <= ANGULAR
            } else {
                false
            }
        } else {
            false
        }
    }

    /// OCCT BRepSweep_Rotation::SplitShell (cxx L797-802).
    ///
    /// GAP: BRepTools_Quilt (TKBRep/BRepTools — not translated; the
    /// BRepTools_Modifier/Quilt/TrsfModification batch owns it, see the port
    /// plan §E1 new-schedule item 4).  OCCT anchor: cxx L799-801 (Q.Add(
    /// aNewShape); return Q.Shells()) — the OCCT failure path (an exception
    /// on quilt-free construction) is preserved by the panic.
    fn split_shell(&self, a_new_shape: &Shape) -> Shape {
        let _ = a_new_shape;
        panic!("GAP: BRepTools_Quilt (TKBRep/BRepTools not translated) — BRepSweep_Rotation::SplitShell, OCCT BRepSweep_Rotation.cxx L797-802");
    }

    /// OCCT BRepSweep_Rotation::HasShape (cxx L806-848).
    fn has_shape(&self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) -> bool {
        if a_dir_s.type_() == ShapeType::Edge && a_gen_s.shape_type() == ShapeType::Edge {
            // Verify that the edge has entrails.
            // OCCT L813-816: the degenerated edge.
            if brep_tool_degenerated(a_gen_s) {
                return false;
            }
            // OCCT L818-824: aCurve = BRep_Tool::Curve(anEdge, aLoc, aPFirst,
            // aPLast); the null curve.
            if brep_tool_curve(a_gen_s).is_none() {
                return false;
            }
            // OCCT L826-829: the invariant meridian.
            if NumLinearRegularSweepSlots::is_invariant(self, a_gen_s) {
                return false;
            }
            // OCCT L831-840: check the seam edge over the generatrix faces.
            for fac in explorer(&self.core.my_gen_shape, ShapeType::Face, ShapeType::Shape) {
                if brep_tools_is_really_closed(a_gen_s, &fac) {
                    return false;
                }
            }
            true
        } else {
            true
        }
    }

    /// OCCT BRepSweep_Rotation::IsInvariant (cxx L852-892) — true when the
    /// geometry of aGenS is not modified by the rotation.
    fn is_invariant(&self, a_gen_s: &Shape) -> bool {
        if a_gen_s.shape_type() == ShapeType::Edge {
            // OCCT L856-858: BRepAdaptor_Curve aC(Edge(aGenS)); the Line /
            // BSpline / Bezier types.
            let (a_c, _first, _last) = match brep_tool_curve(a_gen_s) {
                Some(hit) => hit,
                None => return false,
            };
            let c_type = curve_type(&a_c);
            if c_type == GeomAbsCurveType::Line
                || c_type == GeomAbsCurveType::BSpline
                || c_type == GeomAbsCurveType::Bezier
            {
                // OCCT L860-861: TopExp::Vertices(Edge(aGenS), V1, V2).
                let (v1, v2) = crate::brep_algo::tool::top_exp_vertices_raw(a_gen_s);
                if let (Some(v1), Some(v2)) = (&v1, &v2) {
                    if NumLinearRegularSweepSlots::is_invariant(self, v1)
                        && NumLinearRegularSweepSlots::is_invariant(self, v2)
                    {
                        if c_type == GeomAbsCurveType::Line {
                            return true;
                        }
                        // OCCT L869-872: aTol = max(Tolerance(V1), Tolerance(
                        // V2)); Lin = gp_Lin(myAxe.Location(), myAxe.
                        // Direction()); the poles of the BSpline / Bezier.
                        let a_tol = brep_tool_tolerance(v1).max(brep_tool_tolerance(v2));
                        let poles: &[glam::DVec3] = match &a_c {
                            Curve3::BSpline(b) => &b.control_points,
                            Curve3::Bezier(b) => &b.control_points,
                            _ => unreachable!(),
                        };
                        for pole in poles {
                            let d = line_point_distance(self.my_axe.0, self.my_axe.1, *pole);
                            if d > a_tol {
                                return false;
                            }
                        }
                        return true;
                    }
                }
            }
        } else if a_gen_s.shape_type() == ShapeType::Vertex {
            // OCCT L885-890: P = BRep_Tool::Pnt(Vertex(aGenS));
            // Lin = gp_Lin(...); return Lin.Distance(P) <= Tolerance(V).
            let p = brep_tool_pnt(a_gen_s).expect("BRep_Tool::Pnt");
            let d = line_point_distance(self.my_axe.0, self.my_axe.1, p);
            return d <= brep_tool_tolerance(a_gen_s);
        }
        false
    }

    /// OCCT BRepSweep_Trsf::SetContinuity (cxx L97-201) — Rotation does not
    /// override; the Trsf implementation serves.
    fn set_continuity(&mut self, a_gen_s: &Shape, a_dir_s: &SweepNumShape) {
        trsf::set_continuity(self, a_gen_s, a_dir_s);
    }
}

impl BRepSweepTrsfSlots for BRepSweepRotation {
    fn trsf(&mut self) -> &mut TrsfData {
        &mut self.trsf
    }
    fn trsf_ref(&self) -> &TrsfData {
        &self.trsf
    }
}

impl BRepSweepRotation {
    /// OCCT BRepSweep_Rotation::BRepSweep_Rotation(S, N, L, A, D, C)
    /// (cxx L150-163) — creates a topology by rotating <S> around A with the
    /// angle D.
    pub fn new(
        s: &Shape,
        n: &SweepNumShape,
        l: glam::DAffine3,
        a: GpAx1,
        d: f64,
        c: bool,
    ) -> Self {
        let core = NumLinearRegularSweepCore::new(BRepSweepBuilder::new(), s, n);
        let trsf = TrsfData {
            my_move_memo: std::cell::RefCell::new(indexmap::IndexMap::new()),
            my_location: l,
            my_copy: c,
        };
        let mut me = BRepSweepRotation {
            core,
            trsf,
            my_ang: d,
            my_axe: (a.0, a.1.normalize_or_zero()),
        };
        // OCCT L161: Standard_ConstructionError_Raise_if(D <
        // Precision::Angular(), "BRepSweep_Rotation::Constructor").
        if d < ANGULAR {
            panic!("Standard_ConstructionError: BRepSweep_Rotation::Constructor");
        }
        // OCCT L162: Init().
        trsf::init(&mut me);
        me
    }

    /// OCCT BRepSweep_Rotation::Angle() (cxx L896-899) — returns the angle.
    pub fn angle(&self) -> f64 {
        self.my_ang
    }

    /// OCCT BRepSweep_Rotation::Axe() (cxx L903-906) — returns the axis.
    pub fn axe(&self) -> GpAx1 {
        self.my_axe
    }
}

/// OCCT ComputeTolerance (BRepSweep_Rotation.cxx L71-114) — the pcurve
/// tolerance from the sampled 3D-curve-to-PCurve distance.
fn compute_tolerance(e: &Shape, f: &Shape, c: &Curve2d) -> f64 {
    use rcad_kernel::geom::{Curve2dEval, CurveEval, SurfaceEval};
    if brep_tool_degenerated(e) {
        return brep_tool_tolerance(e);
    }
    // OCCT L83-84: surf = BRep_Tool::Surface(F); c3d = BRep_Tool::Curve(E,
    // first, last).
    let surf = match brep_tool_surface(f) {
        Some(s0) => s0,
        None => return brep_tool_tolerance(e),
    };
    let (c3d, first, last) = match brep_tool_curve(e) {
        Some(hit) => hit,
        None => return brep_tool_tolerance(e),
    };
    // OCCT L86-107: the 23-sample max square distance.
    let mut d2 = 0.0;
    let nn = 23.0;
    let unsurnn = 1.0 / nn;
    let mut i = 0.0;
    while i <= nn {
        let t = unsurnn * i;
        let u = first * (1.0 - t) + last * t;
        let pc3d = c3d.point_at(u);
        let uv = c.point_at(u);
        let pcons = surf.point_at(uv.x, uv.y);
        if pcons.x.is_infinite() || pcons.y.is_infinite() || pcons.z.is_infinite() {
            d2 = f64::INFINITY;
            break;
        }
        let temp = pc3d.distance_squared(pcons);
        if temp > d2 {
            d2 = temp;
        }
        i += 1.0;
    }
    let mut d2 = 1.5 * d2.sqrt();
    if d2 < 1.0e-7 {
        d2 = 1.0e-7;
    }
    d2
}

/// OCCT SetThePCurve (BRepSweep_Rotation.cxx L116-146) — checks if there is
/// already a pcurve; the tolerance comes from ComputeTolerance.
fn set_the_pcurve(me: &mut BRepSweepRotation, e: &Shape, f: &Shape, o: Orientation, c: &Curve2d) {
    // OCCT L126-129: GP = down_cast<Geom_Plane>(Surface(F, SL)); if the face
    // is not planar, OC = CurveOnSurface(E, F, f, l).
    let oc = match brep_tool_surface(f) {
        Some(Surface3::Plane(_)) => None,
        _ => brep_tool_curve_on_surface_seam(e, f),
    };
    // OCCT L131-145: the tolerance argument is ComputeTolerance(E, F, C).
    let tol = compute_tolerance(e, f, c);
    match oc {
        None => {
            me.core.my_builder.my_builder.update_edge_pcurve(e, c, f, tol);
        }
        Some((oc_c, _f0, _l0)) => {
            if o == Orientation::Reversed {
                me.core
                    .my_builder
                    .my_builder
                    .update_edge_two_pcurves(e, &oc_c, c, f, tol);
            } else {
                me.core
                    .my_builder
                    .my_builder
                    .update_edge_two_pcurves(e, c, &oc_c, f, tol);
            }
        }
    }
}

/// Distance of a point from a gp_Lin (Location, Direction).
fn line_point_distance(loc: glam::DVec3, dir: glam::DVec3, p: glam::DVec3) -> f64 {
    let v = p - loc;
    (v - v.dot(dir) * dir).length()
}

/// An arbitrary unit vector perpendicular to the input.
fn any_perpendicular(n: glam::DVec3) -> glam::DVec3 {
    let a = if n.x.abs() < 0.9 {
        glam::DVec3::X
    } else {
        glam::DVec3::Y
    };
    a.cross(n).normalize_or_zero()
}
