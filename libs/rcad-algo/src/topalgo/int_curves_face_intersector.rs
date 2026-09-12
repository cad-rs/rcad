// OCCT IntCurvesFace_Intersector — 1:1 translation of
//   $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/IntCurvesFace/IntCurvesFace_Intersector.hxx (L41-157)
//   $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/IntCurvesFace/IntCurvesFace_Intersector.lxx (L22-87)
//   $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/IntCurvesFace/IntCurvesFace_Intersector.cxx (L44-559)
//
// OCCT inheritance chain: Standard_Transient only (no behavioural base).
//
// Architecture differences (referenced from the affected functions):
// 1. `myTopolTool` — OCCT builds `BRepTopAdaptor_TopolTool(Hsurface)`
//    directly from the adaptor handle (cxx L139).  The rcad
//    topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool resolves every
//    topological accessor through a `BRep` pool, so the face graph is adopted
//    into a standalone pool (topo_builder::brep_from_shape = the
//    TopoDS_Builder::MakeShape/Add glue) and the pool Arc is kept in
//    `my_topology` (the rcad carrier of the OCCT global TShape graph).
//    OCCT anchors: BRepTopAdaptor_TopolTool.cxx L71-94 (Initialize) and
//    L174-187 (Classify).
// 2. `Hsurface` — OCCT `BRepAdaptor_Surface` (TKBRep,
//    BRepAdaptor_Surface.cxx L28-53 Initialize).  The rcad build carrying the
//    HSurfaceTool / PSurfaceTool specialisations consumed by
//    IntCurveSurface_HInter currently lives in
//    bop::int_tools::bean_face_intersector::BRepAdaptorSurface (the bop
//    IntTools translation), so this TKTopAlgo class imports it from there
//    instead of from a TKTopAlgo BRepAdaptor.  `Initialize(Face, aRestr)`
//    maps to `BRepAdaptorSurface::with_uv_bounds(BRep_Tool::Surface(F),
//    BRepTools::UVBounds(F))` / `::new(BRep_Tool::Surface(F))`.
// 3. `IntCurveSurface_HInter` (TKGeomAlgo) — the rcad assembly
//    (geomalgo::int_curve_surface::inter_impl driven through the HInterHost
//    callbacks of bop::int_tools::hinter_adaptor) is instantiated for the
//    adaptors above as
//    bop::int_tools::bean_face_intersector::IntCurveSurfaceHInter.  Its
//    public surface exposes Perform(Curve, Surface) only; the
//    polygon / polyhedron overloads of Inter.pxx reach the same engine
//    through the `pub(crate)` inter_impl entry points, exactly as the OCCT
//    member functions Perform(Curve, Polygon, Surface, Polyhedron) and
//    Perform(Curve, Polygon, Surface, Polyhedron, BndBSB) call them.
// 4. `HICS.Point(index)` — the bop re-host's IntCurveSurfaceIntersectionPoint
//    carries (w, u, v) only, while the OCCT
//    IntCurveSurface_IntersectionPoint also carries `Pnt()` (read by the
//    OCC578 nearest-edge branch at cxx L275).  The base subobject accessor
//    `hics.base.point(index)` is used — it is the OCCT
//    IntCurveSurface_Intersection base subobject of HInter, and its
//    Accessor returns the full IntersectionPoint.
// 5. `Adaptor3d_Curve` handles — OCCT carries them as
//    `occ::handle<Adaptor3d_Curve>` created by `GeomAdaptor_Curve::Load`.
//    The rcad counterpart consumed by IntCurveSurface_HInter is
//    bop::int_tools::bean_face_intersector::BRepAdaptorCurve, built by
//    `BRepAdaptorCurve::new(Curve3)` (the `Load` equivalent, first/last
//    parameter from the curve itself).
// 6. `Bnd_Box` / `NCollection_Array1` — the rcad BndBox / `Vec<f64>` pairs.
//    OCCT's 1-based collections are emulated with a leading unused slot where
//    the OCCT index arithmetic is load bearing (ComputeSamplePars,
//    SeqPnt/mySeqState insertion order).
// 7. `TopoDS_Face` / `TopoDS_Edge` casts are identity on the rcad Shape
//    handle; `BRep_Tool::*` reads are brep_algo::tool re-hosts; the feat
//    pipeline shapes carry identity locations (the same architecture
//    difference recorded in feat::loc_ope_cs_intersector).

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{Curve3, Line3, Plane, Point3, Surface3, Vec3};
use rcad_kernel::math::bnd::{BndBox, BoundSortBox};
use rcad_kernel::math::gp::Lin;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, State, TShape, TVertexData};

use crate::bop::int_tools::bean_face_intersector::{
    BRepAdaptorCurve, BRepAdaptorSurface, IntCurveSurfaceHInter,
};
use crate::bop::int_tools::hinter_adaptor::{BRepAdaptorCurveTool, BRepAdaptorSurfaceTool};
use crate::brep_algo::tool::{brep_tool_curve, brep_tool_surface, brep_tool_tolerance, explorer};
use crate::geomalgo::geom_api_project_point_on_curve::GeomAPIProjectPointOnCurve;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};
use crate::geomalgo::int_curve_surface::{HCurveTool, HSurfaceTool, IntersectionPoint};
use crate::geomalgo::intf::IntfTool;
use crate::topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool;

/// OCCT `RealToInt` (Standard_Real.hxx L319-330) — truncation toward zero,
/// clamped to the int range (which Rust's `as i32` saturating cast performs).
fn real_to_int(the_value: f64) -> i32 {
    the_value as i32
}

/// OCCT static ComputeSamplePars (cxx L45-111) — the parameter grid of the
/// surface's C2 intervals subdivided proportionally to (nbsu, nbsv).
///
/// The returned vectors keep the OCCT 1-based indexing (`upars[i]` is
/// `UPars(i)`); index 0 is unused.
fn compute_sample_pars(
    h_surface: &BRepAdaptorSurface,
    nbsu: i32,
    nbsv: i32,
) -> (Vec<f64>, Vec<f64>) {
    // OCCT L51-56.
    let nb_u_ints = <BRepAdaptorSurfaceTool as HSurfaceTool>::nb_u_intervals(h_surface, rcad_kernel::math::GeomAbsShape::C2);
    let nb_v_ints = <BRepAdaptorSurfaceTool as HSurfaceTool>::nb_v_intervals(h_surface, rcad_kernel::math::GeomAbsShape::C2);
    let mut u_ints = vec![0.0; nb_u_ints + 2];
    let mut v_ints = vec![0.0; nb_v_ints + 2];
    <BRepAdaptorSurfaceTool as HSurfaceTool>::u_intervals(
        h_surface,
        &mut u_ints[1..=nb_u_ints + 1],
        rcad_kernel::math::GeomAbsShape::C2,
    );
    <BRepAdaptorSurfaceTool as HSurfaceTool>::v_intervals(
        h_surface,
        &mut v_ints[1..=nb_v_ints + 1],
        rcad_kernel::math::GeomAbsShape::C2,
    );
    //
    // OCCT L58-59.
    let mut nb_u_sub_ints = vec![0i32; nb_u_ints + 1];
    let mut nb_v_sub_ints = vec![0i32; nb_v_ints + 1];
    //
    // OCCT L61-71.
    let mut ind: i32;
    let nb_u: i32;
    let nb_v: i32;
    let mut t: f64;
    let mut dt: f64;
    t = u_ints[nb_u_ints + 1] - u_ints[1];
    t = 1. / t;
    let mut nb_u_acc = 0i32;
    for i in 1..=nb_u_ints {
        dt = u_ints[i + 1] - u_ints[i];
        nb_u_sub_ints[i] = real_to_int(nbsu as f64 * dt * t) + 1;
        nb_u_acc += nb_u_sub_ints[i];
    }
    nb_u = nb_u_acc;
    t = v_ints[nb_v_ints + 1] - v_ints[1];
    t = 1. / t;
    let mut nb_v_acc = 0i32;
    for i in 1..=nb_v_ints {
        dt = v_ints[i + 1] - v_ints[i];
        nb_v_sub_ints[i] = real_to_int(nbsv as f64 * dt * t) + 1;
        nb_v_acc += nb_v_sub_ints[i];
    }
    nb_v = nb_v_acc;
    //
    // OCCT L81-82: UPars/VPars arrays 1..NbU+1 / 1..NbV+1.
    let mut u_pars = vec![0.0; (nb_u + 1) as usize + 1];
    let mut v_pars = vec![0.0; (nb_v + 1) as usize + 1];
    //
    // OCCT L84-96.
    ind = 1;
    for i in 1..=nb_u_ints {
        u_pars[ind as usize] = u_ints[i];
        ind += 1;
        dt = (u_ints[i + 1] - u_ints[i]) / nb_u_sub_ints[i] as f64;
        t = u_ints[i];
        for _j in 1..nb_u_sub_ints[i] {
            t += dt;
            u_pars[ind as usize] = t;
            ind += 1;
        }
    }
    u_pars[ind as usize] = u_ints[nb_u_ints + 1];
    //
    // OCCT L98-110.
    ind = 1;
    for i in 1..=nb_v_ints {
        v_pars[ind as usize] = v_ints[i];
        ind += 1;
        dt = (v_ints[i + 1] - v_ints[i]) / nb_v_sub_ints[i] as f64;
        t = v_ints[i];
        for _j in 1..nb_v_sub_ints[i] {
            t += dt;
            v_pars[ind as usize] = t;
            ind += 1;
        }
    }
    v_pars[ind as usize] = v_ints[nb_v_ints + 1];
    (u_pars, v_pars)
}

/// OCCT `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)` as carried by the
/// rcad face data (`TFaceData.uv_domain`, the face window used by
/// topalgo::brep_adaptor::surface::BRepAdaptorSurface::initialize_face).
fn brep_tools_uv_bounds(f: &Shape) -> Option<[f64; 4]> {
    match f.data.as_ref() {
        TShape::Face(fd) => fd.uv_domain,
        _ => None,
    }
}

/// The `BRepAdaptor_Surface()` undefined-adaptor carrier (OCCT
/// BRepAdaptor_Surface.cxx L20-24 leaves the plane member unset; rcad uses an
/// XY plane at the origin — the same placeholder the TKTopAlgo
/// topalgo::brep_adaptor::surface::BRepAdaptorSurface::new() carries).
fn undefined_surface() -> Surface3 {
    Surface3::Plane(Plane::new(Point3::ZERO, Vec3::new(0.0, 0.0, 1.0)))
}

/// Adopt the face graph into a standalone `BRep` pool (architecture
/// difference #1): the `BRepTopAdaptor_TopolTool` carrier resolves every
/// topology accessor through a pool index, so the face, its wires and their
/// edges must live at their own flat indices.  This is the
/// `topo_builder::brep_from_shape` idiom of feat::loc_ope_find_edges_in_face
/// with one correction: sub-shapes built pool-free carry
/// `Shape::index == usize::MAX` (`Shape::is_null`), for which the
/// `while tshapes.len() <= index` slot stretch cannot terminate.  The
/// vertices are the only such nodes of a face graph; every consumer on this
/// path (BRepTopAdaptor_TopolTool::Initialize, BRepTopAdaptor_FClass2d,
/// BRep_Tool::CurveOnSurface) reads the wires and edges only, so pool-free
/// nodes are carried over as pool-free leaves.
fn adopt_face_graph(face: &Shape) -> BRep {
    let mut pool = BRep::new();
    let mut visited: std::collections::HashSet<u64> = std::collections::HashSet::new();
    adopt_shape(&mut pool, face, &mut visited);
    pool
}

/// One TopoDS_Builder::Add / TopoDS_Iterator step of [`adopt_face_graph`].
fn adopt_shape(pool: &mut BRep, s: &Shape, visited: &mut std::collections::HashSet<u64>) {
    if s.index == usize::MAX || !visited.insert(s.ptr_id()) {
        return;
    }
    if pool.tshapes.len() <= s.index {
        let dummy = Arc::new(TShape::Vertex(TVertexData {
            my_shapes: Vec::new(),
            flags: 0,
            point: Point3::ZERO,
            tolerance: 0.0,
            points: Vec::new(),
        }));
        while pool.tshapes.len() <= s.index {
            pool.tshapes.push(dummy.clone());
        }
    }
    pool.tshapes[s.index] = s.data.clone();
    for c in crate::brep_algo::tool::sub_shapes(s) {
        adopt_shape(pool, &c, visited);
    }
}

/// OCCT IntCurvesFace_Intersector (hxx L41-157).
pub(crate) struct IntCurvesFaceIntersector {
    /// OCCT myTopolTool (hxx L142).
    my_topol_tool: BRepTopolTool,
    /// rcad carrier of the topology the tool resolves through
    /// (architecture difference #1); OCCT has no such member because the
    /// TShape graph is global.
    my_topology: Arc<BRep>,
    /// OCCT Hsurface (hxx L143).
    h_surface: BRepAdaptorSurface,
    /// OCCT Tol (hxx L144).
    tol: f64,
    /// OCCT SeqPnt (hxx L145).
    seq_pnt: Vec<IntersectionPoint>,
    /// OCCT mySeqState (hxx L146).
    my_seq_state: Vec<i32>,
    /// OCCT done (hxx L147).
    done: bool,
    /// OCCT myReady (hxx L148).
    my_ready: bool,
    /// OCCT nbpnt (hxx L149).
    nbpnt: i32,
    /// OCCT face (hxx L150).
    face: Shape,
    /// OCCT myPolyhedron (hxx L151).
    my_polyhedron: Option<ThePolyhedronOfHInter>,
    /// OCCT myBndBounding (hxx L152).
    my_bnd_bounding: Option<BoundSortBox>,
    /// OCCT myUseBoundTol (hxx L153).
    my_use_bound_tol: bool,
    /// OCCT myIsParallel (hxx L154).
    my_is_parallel: bool,
}

impl IntCurvesFaceIntersector {
    /// OCCT IntCurvesFace_Intersector(const TopoDS_Face& F, const double aTol,
    /// const bool aRestr = true, const bool UseBToler = true) (cxx L123-224) —
    /// the rcad 2-argument form carries the OCCT defaulted parameters.
    pub(crate) fn new(f: &Shape, a_tol: f64) -> Self {
        IntCurvesFaceIntersector::new_full(f, a_tol, true, true)
    }

    /// OCCT IntCurvesFace_Intersector(const TopoDS_Face& Face,
    /// const double aTol, const bool aRestr, const bool UseBToler)
    /// (cxx L123-224) — the constructor with every OCCT parameter spelled out.
    pub(crate) fn new_full(f: &Shape, a_tol: f64, a_restr: bool, use_b_toler: bool) -> Self {
        // Architecture difference #1: the pool the TopolTool resolves through.
        let pool = adopt_face_graph(f);
        // OCCT L128-133: the member initialiser list.
        let mut a_this = IntCurvesFaceIntersector {
            my_topol_tool: BRepTopolTool::new(Arc::new(BRep::new())),
            my_topology: Arc::new(pool),
            h_surface: BRepAdaptorSurface::new(undefined_surface()),
            tol: a_tol,
            seq_pnt: Vec::new(),
            my_seq_state: Vec::new(),
            done: false,
            my_ready: false,
            nbpnt: 0,
            face: f.clone(),
            my_polyhedron: None,
            my_bnd_bounding: None,
            my_use_bound_tol: use_b_toler,
            my_is_parallel: false,
        };

        // OCCT L135-138: BRepAdaptor_Surface surface;
        // face = Face; surface.Initialize(Face, aRestr);
        // Hsurface = new BRepAdaptor_Surface(surface);
        a_this.h_surface = match brep_tool_surface(f) {
            Some(s) => match if a_restr { brep_tools_uv_bounds(f) } else { None } {
                Some(d) => BRepAdaptorSurface::with_uv_bounds(s, d),
                None => BRepAdaptorSurface::new(s),
            },
            // BRepAdaptor_Surface::Initialize leaves the adaptor undefined when
            // BRep_Tool::Surface returns a null handle (cxx L32-37).
            None => BRepAdaptorSurface::new(undefined_surface()),
        };
        // OCCT L139: myTopolTool = new BRepTopAdaptor_TopolTool(Hsurface).
        let mut tool = BRepTopolTool::new(a_this.my_topology.clone());
        tool.initialize_surface(f);
        a_this.my_topol_tool = tool;
        // OCCT L141: GeomAbs_SurfaceType SurfaceType = GetType(Hsurface).
        let surface_type = <BRepAdaptorSurfaceTool as HSurfaceTool>::get_type(&a_this.h_surface);
        // OCCT L142-145.
        if surface_type != crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane
            && surface_type != crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
            && surface_type != crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone
            && surface_type != crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere
            && surface_type != crate::geomalgo::int_patch::GeomAbsSurfaceType::Torus
        {
            // OCCT L146-151.
            let u0 = a_this.h_surface.first_u_parameter();
            let u1 = a_this.h_surface.last_u_parameter();
            let v0 = a_this.h_surface.first_v_parameter();
            let v1 = a_this.h_surface.last_v_parameter();
            //
            // OCCT L153-154.
            let mut nbsu = a_this.my_topol_tool.nb_samples_u();
            let mut nbsv = a_this.my_topol_tool.nb_samples_v();
            //
            // OCCT L156-157.
            let a_u_res = <BRepAdaptorSurfaceTool as HSurfaceTool>::u_resolution(&a_this.h_surface, 1.0);
            let a_v_res = <BRepAdaptorSurfaceTool as HSurfaceTool>::v_resolution(&a_this.h_surface, 1.0);

            // OCCT L159-165.
            let a_tresh = 100.0;
            let mut a_min_samples = 20;
            let a_max_samples = 40;
            let a_max_samples2 = a_max_samples * a_max_samples;
            let d_u = (u1 - u0) / a_u_res;
            let d_v = (v1 - v0) / a_v_res;
            // OCCT L166-181.
            if nbsu < a_min_samples {
                nbsu = a_min_samples;
            }
            if nbsv < a_min_samples {
                nbsv = a_min_samples;
            }
            if nbsu > a_max_samples {
                nbsu = a_max_samples;
            }
            if nbsv > a_max_samples {
                nbsv = a_max_samples;
            }

            // OCCT L183-204.
            if d_u > rcad_kernel::precision::CONFUSION
                && d_v > rcad_kernel::precision::CONFUSION
            {
                if d_u.max(d_v) > d_u.min(d_v) * a_tresh {
                    a_min_samples = 10;
                    nbsu = ((d_u / d_v).sqrt() * a_max_samples as f64) as i32;
                    if nbsu < a_min_samples {
                        nbsu = a_min_samples;
                    }
                    nbsv = a_max_samples2 / nbsu;
                    if nbsv < a_min_samples {
                        nbsv = a_min_samples;
                        nbsu = a_max_samples2 / a_min_samples;
                    }
                }
            } else {
                // OCCT L203: surface has no extension along one of directions.
                return a_this;
            }

            // OCCT L206-207.
            let nb_u_on_s = <BRepAdaptorSurfaceTool as HSurfaceTool>::nb_u_intervals(
                &a_this.h_surface,
                rcad_kernel::math::GeomAbsShape::C2,
            );
            let nb_v_on_s = <BRepAdaptorSurfaceTool as HSurfaceTool>::nb_v_intervals(
                &a_this.h_surface,
                rcad_kernel::math::GeomAbsShape::C2,
            );

            // OCCT L209-221.
            if nb_u_on_s > 1 || nb_v_on_s > 1 {
                let (u_pars, v_pars) = compute_sample_pars(&a_this.h_surface, nbsu, nbsv);
                a_this.my_polyhedron = Some(ThePolyhedronOfHInter::new_tool_params::<
                    BRepAdaptorSurface,
                    BRepAdaptorSurfaceTool,
                >(
                    &a_this.h_surface, &u_pars[1..], &v_pars[1..],
                ));
            } else {
                a_this.my_polyhedron = Some(ThePolyhedronOfHInter::new_tool::<
                    BRepAdaptorSurface,
                    BRepAdaptorSurfaceTool,
                >(
                    &a_this.h_surface,
                    nbsu as usize,
                    nbsv as usize,
                    u0,
                    v0,
                    u1,
                    v1,
                ));
            }
        }
        // OCCT L223.
        a_this.my_ready = true;
        a_this
    }

    /// OCCT IntCurvesFace_Intersector::SurfaceType() (cxx L116-119).
    pub(crate) fn surface_type(&self) -> crate::geomalgo::int_patch::GeomAbsSurfaceType {
        <BRepAdaptorSurfaceTool as HSurfaceTool>::get_type(&self.h_surface)
    }

    /// OCCT IntCurvesFace_Intersector::InternalCall(const
    /// IntCurveSurface_HInter& HICS, const double parinf, const double
    /// parsup) (cxx L228-363).
    fn internal_call(&mut self, hics: &IntCurveSurfaceHInter, parinf: f64, parsup: f64) {
        // OCCT L232.
        if hics.is_done() && hics.nb_points() > 0 {
            // OCCT L235-237.
            let mut mintol3d = brep_tool_tolerance(&self.face);
            let mut maxtol3d = mintol3d;
            let mintol2d;
            let maxtol2d;
            //
            // OCCT L238-244: TopExp_Explorer anExp(face, TopAbs_EDGE).
            let an_exp = explorer(&self.face, ShapeType::Edge, ShapeType::Shape);
            for e in an_exp.iter() {
                let curtol = brep_tool_tolerance(e);
                mintol3d = mintol3d.min(curtol);
                maxtol3d = maxtol3d.max(curtol);
            }
            // OCCT L245-248.
            let minres = <BRepAdaptorSurfaceTool as HSurfaceTool>::u_resolution(&self.h_surface, mintol3d)
                .max(<BRepAdaptorSurfaceTool as HSurfaceTool>::v_resolution(&self.h_surface, mintol3d));
            let maxres = <BRepAdaptorSurfaceTool as HSurfaceTool>::u_resolution(&self.h_surface, maxtol3d)
                .max(<BRepAdaptorSurfaceTool as HSurfaceTool>::v_resolution(&self.h_surface, maxtol3d));
            mintol2d = minres.max(self.tol);
            maxtol2d = maxres.max(self.tol);
            //
            // OCCT L250: occ::handle<BRepTopAdaptor_TopolTool> anAdditionalTool.
            let mut an_additional_tool: Option<BRepTopolTool> = None;
            // OCCT L251.
            for index in (1..=hics.nb_points() as i32).rev() {
                // OCCT L253-254.
                let hics_point_index = hics.base.point(index as usize);
                let puv = DVec2::new(hics_point_index.u(), hics_point_index.v());

                // OCCT L257.
                let mut currentstate = self.my_topol_tool.classify(
                    puv,
                    if !self.my_use_bound_tol { 0.0 } else { mintol2d },
                    true,
                );
                // OCCT L258-287.
                if self.my_use_bound_tol && currentstate == State::Out && maxtol2d > mintol2d {
                    if an_additional_tool.is_none() {
                        // OCCT L262: anAdditionalTool =
                        // new BRepTopAdaptor_TopolTool(Hsurface).
                        let mut tool = BRepTopolTool::new(self.my_topology.clone());
                        tool.initialize_surface(&self.face);
                        an_additional_tool = Some(tool);
                    }
                    // OCCT L264.
                    currentstate =
                        an_additional_tool.as_mut().unwrap().classify(puv, maxtol2d, true);
                    if currentstate == State::On {
                        currentstate = State::Out;
                        // OCCT L266-285: find out the nearest edge and its
                        // tolerance.
                        for an_e in an_exp.iter() {
                            let Some((a_pc, f, l)) = brep_tool_curve(an_e) else {
                                continue;
                            };
                            let a_proj = GeomAPIProjectPointOnCurve::new_point_curve_ranged(
                                hics_point_index.pnt(),
                                &a_pc,
                                f,
                                l,
                            );
                            if a_proj.nb_points() > 0 {
                                if a_proj.lower_distance() <= maxtol3d {
                                    // OCCT L280-282: nearest edge is found,
                                    // state is really ON.
                                    currentstate = State::On;
                                    break;
                                }
                            }
                        }
                    }
                }
                // OCCT L288.
                if currentstate == State::In || currentstate == State::On {
                    // OCCT L290-297.
                    let hicsw = hics_point_index.w();
                    if hicsw >= parinf && hicsw <= parsup {
                        let u = hics_point_index.u();
                        let v = hics_point_index.v();
                        let w = hicsw;
                        let mut transition = hics_point_index.transition();
                        let pnt = hics_point_index.pnt();
                        // OCCT L299: int anIntState = (currentstate == TopAbs_IN) ? 0 : 1;
                        let an_int_state = if currentstate == State::In { 0i32 } else { 1i32 };

                        // OCCT L302-312.
                        if transition != crate::geomalgo::int_curve_surface::TransitionOnCurve::Tangent
                            && self.face.orientation == Orientation::Reversed
                        {
                            if transition
                                == crate::geomalgo::int_curve_surface::TransitionOnCurve::In
                            {
                                transition =
                                    crate::geomalgo::int_curve_surface::TransitionOnCurve::Out;
                            } else {
                                transition =
                                    crate::geomalgo::int_curve_surface::TransitionOnCurve::In;
                            }
                        }
                        // OCCT L313-352: insertion of the point.
                        if self.nbpnt == 0 {
                            // OCCT L316-320.
                            let ppp = IntersectionPoint::with_values(pnt, u, v, w, transition);
                            self.seq_pnt.push(ppp);
                            self.my_seq_state.push(an_int_state);
                        } else {
                            // OCCT L324-336.
                            let mut i = 1;
                            let mut b = self.nbpnt + 1;
                            while i <= self.nbpnt {
                                let pnti = &self.seq_pnt[(i - 1) as usize];
                                let wi = pnti.w();
                                if wi >= w {
                                    b = i;
                                    i = self.nbpnt;
                                }
                                i += 1;
                            }
                            // OCCT L337.
                            let ppp = IntersectionPoint::with_values(pnt, u, v, w, transition);
                            // OCCT L341-350.
                            if b > self.nbpnt {
                                self.seq_pnt.push(ppp);
                                self.my_seq_state.push(an_int_state);
                            } else if b > 0 {
                                self.seq_pnt.insert((b - 1) as usize, ppp);
                                self.my_seq_state.insert((b - 1) as usize, an_int_state);
                            }
                        }

                        // OCCT L354.
                        self.nbpnt += 1;
                    }
                } // OCCT L356: classifier state is IN or ON.
            } // OCCT L357: loop on intersection points.
        } // OCCT L358: HICS.IsDone().
        else if hics.is_done() {
            // OCCT L361.
            self.my_is_parallel = hics.is_parallel();
        }
    }

    /// OCCT IntCurvesFace_Intersector::Perform(const gp_Lin& L,
    /// const double ParMin, const double ParMax) (cxx L367-464).
    pub(crate) fn perform_lin(&mut self, l: &Lin, par_min: f64, par_max: f64) {
        // OCCT L369-374.
        self.done = false;
        if !self.my_ready {
            return;
        }
        self.done = true;
        self.seq_pnt.clear();
        self.my_seq_state.clear();
        self.nbpnt = 0;

        // OCCT L379-384: IntCurveSurface_HInter HICS; Geom_Line geomline(L);
        // GeomAdaptor_Curve LL(geomline); HLL; parinf = ParMin;
        // parsup = ParMax.
        let mut hics = IntCurveSurfaceHInter::new();
        let hll = BRepAdaptorCurve::new(Curve3::Line(Line3::new(l.pos, l.dir)));
        let mut parinf = par_min;
        let mut parsup = par_max;
        //
        // OCCT L386.
        if self.my_polyhedron.is_none() {
            // OCCT L388.
            hics.perform(&hll, &self.h_surface);
        } else {
            // OCCT L392-394.
            let mut bnd_tool = IntfTool::new();
            let mut box_line = BndBox::new();
            let l3 = Line3::new(l.pos, l.dir);
            bnd_tool.lin_box(
                &l3,
                self.my_polyhedron.as_ref().unwrap().bounding(),
                &mut box_line,
            );
            // OCCT L395-398.
            if bnd_tool.nb_segments() == 0 {
                return;
            }
            // OCCT L399-427.
            for nbseg in 1..=bnd_tool.nb_segments() {
                let mut pinf = bnd_tool.begin_param(nbseg);
                let mut psup = bnd_tool.end_param(nbseg);
                let pppp = 0.05 * (psup - pinf);
                pinf -= pppp;
                psup += pppp;
                if (psup - pinf) < 1e-10 {
                    pinf -= 1e-10;
                    psup += 1e-10;
                }
                if nbseg == 1 {
                    parinf = pinf;
                    parsup = psup;
                } else {
                    if parinf > pinf {
                        parinf = pinf;
                    }
                    if parsup < psup {
                        parsup = psup;
                    }
                }
            }
            // OCCT L428-443.
            if parinf > par_max {
                return;
            }
            if parsup < par_min {
                return;
            }
            if parinf < par_min {
                parinf = par_min;
            }
            if parsup > par_max {
                parsup = par_max;
            }
            // OCCT L444-447.
            if parinf > (parsup - 1e-9) {
                return;
            }
            // OCCT L448: IntCurveSurface_ThePolygonOfHInter polygon(HLL,
            // parinf, parsup, 2).
            let polygon = ThePolygonOfHInter::new_tool_range::<BRepAdaptorCurve, BRepAdaptorCurveTool>(
                &hll, parinf, parsup, 2,
            );
            // OCCT L449-460 (OPTIMISATION).
            if self.my_bnd_bounding.is_none() {
                let mut bsb = BoundSortBox::new();
                let ph0 = self.my_polyhedron.as_ref().unwrap();
                bsb.initialize_with_enclosing(
                    ph0.bounding().clone(),
                    ph0.components_bounding().to_vec(),
                );
                self.my_bnd_bounding = Some(bsb);
            }
            // OCCT L457: HICS.Perform(HLL, polygon, Hsurface, *myPolyhedron,
            // *myBndBounding).
            let ph = self.my_polyhedron.as_ref().unwrap();
            crate::geomalgo::int_curve_surface::inter_impl::perform_polygon_polyhedron_bsb::<
                BRepAdaptorCurve,
                BRepAdaptorCurveTool,
                BRepAdaptorSurface,
                BRepAdaptorSurfaceTool,
                IntCurveSurfaceHInter,
            >(
                &hll,
                &polygon,
                &self.h_surface,
                ph,
                self.my_bnd_bounding.as_mut().unwrap(),
                &mut hics,
            );
        }

        // OCCT L463.
        self.internal_call(&hics, parinf, parsup);
    }

    /// OCCT IntCurvesFace_Intersector::Perform(const
    /// occ::handle<Adaptor3d_Curve>& HCu, const double ParMin,
    /// const double ParMax) (cxx L468-527).
    pub(crate) fn perform_hcurve(&mut self, hcu: &BRepAdaptorCurve, par_min: f64, par_max: f64) {
        // OCCT L470-482.
        self.done = false;
        if !self.my_ready {
            return;
        }
        self.done = true;
        self.seq_pnt.clear();
        self.my_seq_state.clear();
        self.nbpnt = 0;
        let mut hics = IntCurveSurfaceHInter::new();

        // OCCT L486-487.
        let mut parinf = par_min;
        let mut parsup = par_max;

        // OCCT L489.
        if self.my_polyhedron.is_none() {
            // OCCT L491.
            hics.perform(hcu, &self.h_surface);
        } else {
            // OCCT L495-508.
            parinf = <BRepAdaptorCurveTool as HCurveTool>::first_parameter(hcu);
            parsup = <BRepAdaptorCurveTool as HCurveTool>::last_parameter(hcu);
            if parinf < par_min {
                parinf = par_min;
            }
            if parsup > par_max {
                parsup = par_max;
            }
            if parinf > (parsup - 1e-9) {
                return;
            }
            // OCCT L509-512.
            let nbs = <BRepAdaptorCurveTool as HCurveTool>::nb_samples(hcu, parinf, parsup);
            let polygon = ThePolygonOfHInter::new_tool_range::<BRepAdaptorCurve, BRepAdaptorCurveTool>(
                hcu, parinf, parsup, nbs,
            );
            // OCCT L513-524 (OPTIMISATION).
            if self.my_bnd_bounding.is_none() {
                let mut bsb = BoundSortBox::new();
                let ph0 = self.my_polyhedron.as_ref().unwrap();
                bsb.initialize_with_enclosing(
                    ph0.bounding().clone(),
                    ph0.components_bounding().to_vec(),
                );
                self.my_bnd_bounding = Some(bsb);
            }
            // OCCT L521: HICS.Perform(HCu, polygon, Hsurface, *myPolyhedron,
            // *myBndBounding).
            let ph = self.my_polyhedron.as_ref().unwrap();
            crate::geomalgo::int_curve_surface::inter_impl::perform_polygon_polyhedron_bsb::<
                BRepAdaptorCurve,
                BRepAdaptorCurveTool,
                BRepAdaptorSurface,
                BRepAdaptorSurfaceTool,
                IntCurveSurfaceHInter,
            >(
                hcu,
                &polygon,
                &self.h_surface,
                ph,
                self.my_bnd_bounding.as_mut().unwrap(),
                &mut hics,
            );
        }
        // OCCT L526.
        self.internal_call(&hics, parinf, parsup);
    }

    /// The rcad carrier of the OCCT `GeomAdaptor_Curve HC; HC->Load(Scur(i))`
    /// step of the callers (LocOpe_CSIntersector.cxx L167-171): the curve is
    /// loaded into the adaptor, then
    /// `Perform(const occ::handle<Adaptor3d_Curve>&, PInf, PSup)` runs.
    pub(crate) fn perform_curve(&mut self, the_curve: &Curve3, par_min: f64, par_max: f64) {
        // OCCT LocOpe_CSIntersector.cxx L170: HC->Load(Scur(i)).
        let hcu = BRepAdaptorCurve::new(the_curve.clone());
        self.perform_hcurve(&hcu, par_min, par_max);
    }

    /// OCCT IntCurvesFace_Intersector::Bounding() const (cxx L530-541).
    pub(crate) fn bounding(&self) -> BndBox {
        if let Some(ph) = &self.my_polyhedron {
            ph.bounding().clone()
        } else {
            BndBox::new()
        }
    }

    /// OCCT IntCurvesFace_Intersector::ClassifyUVPoint(const gp_Pnt2d& Puv)
    /// (cxx L543-547).
    pub(crate) fn classify_uv_point(&mut self, puv: DVec2) -> State {
        self.my_topol_tool.classify(puv, 1e-7, true)
    }

    /// OCCT IntCurvesFace_Intersector::SetUseBoundToler(bool) (cxx L549-552).
    pub(crate) fn set_use_bound_toler(&mut self, use_b_toler: bool) {
        self.my_use_bound_tol = use_b_toler;
    }

    /// OCCT IntCurvesFace_Intersector::GetUseBoundToler() const (cxx L554-557).
    pub(crate) fn get_use_bound_toler(&self) -> bool {
        self.my_use_bound_tol
    }

    // ---------------------------------------------------------------------
    // Inline accessors — IntCurvesFace_Intersector.lxx L22-87.
    // ---------------------------------------------------------------------

    /// OCCT IsDone() (lxx L23-26).
    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT NbPnt() (lxx L29-32).
    pub(crate) fn nb_pnt(&self) -> i32 {
        self.nbpnt
    }

    /// OCCT Pnt(const int i) (lxx L35-38) — 1-based.
    pub(crate) fn pnt(&self, i: i32) -> DVec3 {
        self.seq_pnt[(i - 1) as usize].pnt()
    }

    /// OCCT UParameter(const int i) (lxx L41-44) — 1-based.
    pub(crate) fn u_parameter(&self, i: i32) -> f64 {
        self.seq_pnt[(i - 1) as usize].u()
    }

    /// OCCT VParameter(const int i) (lxx L47-50) — 1-based.
    pub(crate) fn v_parameter(&self, i: i32) -> f64 {
        self.seq_pnt[(i - 1) as usize].v()
    }

    /// OCCT WParameter(const int i) (lxx L53-56) — 1-based.
    pub(crate) fn w_parameter(&self, i: i32) -> f64 {
        self.seq_pnt[(i - 1) as usize].w()
    }

    /// OCCT Transition(const int i) (lxx L59-62) — 1-based.
    ///
    /// Architecture difference: rcad currently carries a second copy of
    /// `IntCurveSurface_TransitionOnCurve` under
    /// topalgo::brep_int_curve_surface::inter; the two in-domain consumers
    /// (feat::loc_ope_curve_shape_intersector, topalgo::brep_class3d::
    /// intersector3d) are typed against that copy, so the accessor converts
    /// the canonical geomalgo enumeration into it.
    pub(crate) fn transition(
        &self,
        i: i32,
    ) -> crate::topalgo::brep_int_curve_surface::inter::TransitionOnCurve {
        use crate::geomalgo::int_curve_surface::TransitionOnCurve as Src;
        use crate::topalgo::brep_int_curve_surface::inter::TransitionOnCurve as Dst;
        match self.seq_pnt[(i - 1) as usize].transition() {
            Src::In => Dst::In,
            Src::Out => Dst::Out,
            Src::Tangent => Dst::Tangent,
        }
    }

    /// OCCT State(const int i) (lxx L70-73) — 1-based.
    pub(crate) fn state(&self, i: i32) -> State {
        if self.my_seq_state[(i - 1) as usize] == 0 {
            State::In
        } else {
            State::On
        }
    }

    /// OCCT IsParallel() (lxx L77-80).
    pub(crate) fn is_parallel(&self) -> bool {
        self.my_is_parallel
    }

    /// OCCT Face() (lxx L82-85).
    pub(crate) fn face(&self) -> &Shape {
        &self.face
    }
}
