//! OCCT BRepAlgo_NormalProjection — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepAlgo/
//!         BRepAlgo_NormalProjection.cxx (L74-677) +
//!         BRepAlgo_NormalProjection.hxx (L34-119).
//!
//! Architecture differences:
//! 1. NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, ShapeMapHasher> ->
//!    HashMap<ShapeKey, (Shape, Shape)> (the entry keeps the key shape);
//!    NCollection_DataMap<TopoDS_Shape, List, ...> -> HashMap<ShapeKey,
//!    (Shape, Vec<Shape>)>.
//! 2. BRepAdaptor_Curve/BRepAdaptor_Surface are the kernel 1:1 bodies
//!    (rcad_kernel::base::proj_lib::brep_adaptor); GeomAdaptor::MakeCurve
//!    (cxx L324/L363) is the identity on the adaptor's loaded 3D curve (the
//!    rcad adaptor IS the geom value).  GAP: the curve-on-surface arm of
//!    GeomAdaptor::MakeCurve (an edge without a 3D curve) is not translated;
//!    the OCCT path derives the 3D image there, rcad skips the solution.
//! 3. GeomAbs_Shape -> rcad_kernel::topods::GeomAbsShape; GeomAbs_CurveType
//!    maps to the kernel base::proj_lib::CurveType (IsElementary).
//! 4. BRepLib_MakeWire (TKBRep connected-wire builder) — pending;
//!    IsDone() is false so BuildWire takes the OCCT not-a-wire exit.
//!    GAP: closes with the BRepLib batch.
//! 5. BRepAlgoAPI_Section (cxx L465) -> bop::brep_algo_api::SectionOp (the
//!    real BOPAlgo_Section-driven implementation).
//! 6. BRepTopAdaptor_FClass2d -> topalgo::brep_top_adaptor::fclass2d::
//!    FClass2d over FaceShapeSource (loc_ope_wires_on_shape_b.rs #8).
//! 7. BRep_Builder edits are in-place Arc::make_mut mutations (tool.rs);
//!    OCCT TShape sharing makes them visible through every handle, rcad
//!    mutates the owning copy.

use std::sync::Arc;

use crate::brep_algo::tool::{
    brep_tool_surface, builder_add_compound_shape, builder_make_compound, builder_make_edge,
    builder_make_vertex, builder_range_edge_on_face, builder_set_degenerated,
    builder_update_edge_pcurve, builder_update_vertex_point_tol, builder_update_vertex_tol,
    explorer, oriented, shape_key, top_exp_vertices_raw, ShapeKey,
};
use crate::geomalgo::approx_curve_on_surface::ApproxCurveOnSurface;
use crate::geomalgo::proj_lib_h_comp_projected_curve::{
    geom2d_adaptor_curve, new_h_comp_projected_curve, CompProjectedCurve,
};
use rcad_kernel::base::proj_lib::adaptor::Curve2dHandle;
use rcad_kernel::base::proj_lib::brep_adaptor::{BRepAdaptorCurve, BRepAdaptorSurface};
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::Adaptor3dCurveGeom;
use rcad_kernel::base::proj_lib::Adaptor3dSurface;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, TrimmedCurve2, BSplineCurve2,
};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{GeomAbsShape, Orientation, ShapeType};
use std::collections::HashMap;

use crate::bop::brep_algo_api::{Algo, SectionOp};

/// OCCT BRepAlgo_NormalProjection (BRepAlgo_NormalProjection.hxx L34-119) —
/// makes the projection of a wire on a shape.
pub struct BRepAlgoNormalProjection {
    /// OCCT: myShape.
    my_shape: Shape,
    /// OCCT: myIsDone.
    my_is_done: bool,
    /// OCCT: myTol3d.
    my_tol3d: f64,
    /// OCCT: myTol2d.
    my_tol2d: f64,
    /// OCCT: myMaxDist.
    my_max_dist: f64,
    /// OCCT: myWith3d.
    my_with3d: bool,
    /// OCCT: myContinuity.
    my_continuity: GeomAbsShape,
    /// OCCT: myMaxDegree.
    my_max_degree: usize,
    /// OCCT: myMaxSeg.
    my_max_seg: usize,
    /// OCCT: myFaceBounds.
    my_face_bounds: bool,
    /// OCCT: myToProj.
    my_to_proj: Shape,
    /// OCCT: myAncestorMap.
    my_ancestor_map: HashMap<ShapeKey, (Shape, Shape)>,
    /// OCCT: myCorresp.
    my_corresp: HashMap<ShapeKey, (Shape, Shape)>,
    /// OCCT: myDescendants.
    my_descendants: HashMap<ShapeKey, (Shape, Vec<Shape>)>,
    /// OCCT: myRes.
    my_res: Shape,
}

impl Default for BRepAlgoNormalProjection {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAlgoNormalProjection {
    /// OCCT BRepAlgo_NormalProjection::BRepAlgo_NormalProjection()
    /// (cxx L74-83).
    pub fn new() -> Self {
        let mut p = BRepAlgoNormalProjection {
            my_shape: Shape::null(),
            my_is_done: false,
            my_tol3d: 0.0,
            my_tol2d: 0.0,
            my_max_dist: -1.0,
            my_with3d: true,
            my_continuity: GeomAbsShape::C0,
            my_max_degree: 0,
            my_max_seg: 0,
            my_face_bounds: true,
            my_to_proj: Shape::null(),
            my_ancestor_map: HashMap::new(),
            my_corresp: HashMap::new(),
            my_descendants: HashMap::new(),
            my_res: Shape::null(),
        };
        // OCCT L80-81: BB.MakeCompound(TopoDS::Compound(myToProj)).
        p.my_to_proj = builder_make_compound();
        p.set_default_params();
        p
    }

    /// OCCT BRepAlgo_NormalProjection::BRepAlgo_NormalProjection(S)
    /// (cxx L87-97).
    pub fn new_with_shape(s: &Shape) -> Self {
        let mut p = Self::new();
        p.init(s);
        p
    }

    /// OCCT BRepAlgo_NormalProjection::Init(S) (cxx L101-104).
    pub fn init(&mut self, s: &Shape) {
        self.my_shape = s.clone();
    }

    /// OCCT BRepAlgo_NormalProjection::Add(ToProj) (cxx L108-112) — adds an
    /// edge or a wire to the list of shape to project.
    pub fn add(&mut self, to_proj: &Shape) {
        builder_add_compound_shape(&mut self.my_to_proj, to_proj);
    }

    /// OCCT BRepAlgo_NormalProjection::SetParams(...) (cxx L116-127) — sets
    /// the parameters used for computation.
    pub fn set_params(
        &mut self,
        tol3_d: f64,
        tol2_d: f64,
        internal_continuity: GeomAbsShape,
        max_degree: usize,
        max_seg: usize,
    ) {
        self.my_tol3d = tol3_d;
        self.my_tol2d = tol2_d;
        self.my_continuity = internal_continuity;
        self.my_max_degree = max_degree;
        self.my_max_seg = max_seg;
    }

    /// OCCT BRepAlgo_NormalProjection::SetDefaultParams() (cxx L131-138).
    pub fn set_default_params(&mut self) {
        self.my_tol3d = 1.0e-4;
        self.my_tol2d = self.my_tol3d.powf(2.0 / 3.0);
        self.my_continuity = GeomAbsShape::C1;
        self.my_max_degree = 14;
        self.my_max_seg = 16;
    }

    /// OCCT BRepAlgo_NormalProjection::SetLimit(FaceBounds) (cxx L142-145).
    pub fn set_limit(&mut self, face_bounds: bool) {
        self.my_face_bounds = face_bounds;
    }

    /// OCCT BRepAlgo_NormalProjection::SetMaxDistance(MaxDist)
    /// (cxx L149-152).
    pub fn set_max_distance(&mut self, max_dist: f64) {
        self.my_max_dist = max_dist;
    }

    /// OCCT BRepAlgo_NormalProjection::Compute3d(With3d) (cxx L156-159).
    pub fn compute3d(&mut self, with3d: bool) {
        self.my_with3d = with3d;
    }

    // NOTE: the "value assigned is never read" warnings inside build are
    // intentional — they mirror the OCCT reference flow (cxx L245 prj /
    // L272-273 Only2d-Only3d are assigned per iteration before use).
    #[allow(unused_assignments)]

    /// OCCT BRepAlgo_NormalProjection::Build() (cxx L163-584).
    pub fn build(&mut self) {
        // (OCCT_DEBUG_CHRONO blocks L165-179 and L551-583 not translated.)
        self.my_is_done = false;

        // OCCT L182-191: the Edges/Faces sequences and the locals.
        let mut edges: Vec<Shape> = Vec::new();
        let mut faces: Vec<Shape> = Vec::new();
        let mut descen_list: Vec<Shape> = Vec::new();
        let mut udeb = 0.0f64;
        let mut ufin = 0.0f64;
        let mut vertex_res: Shape;
        let mut only3d: bool;
        let mut only2d: bool;

        // for isoparametric cases (OCCT L193-200).
        let mut poles: [glam::DVec2; 2] = [glam::DVec2::ZERO; 2];
        let mut knots: [f64; 2] = [0.0; 2];
        let deg: usize = 1;
        let mults: [usize; 2] = [deg + 1, deg + 1];
        let _ = &mults; // consumed by the OCCT Geom2d_BSplineCurve ctor
                         // (the rcad BSplineCurve2 carries the expanded knots)

        // OCCT L203-211: the explorations (the NbEdges/NbFaces counters
        // increment in lockstep with the sequences).
        for e in explorer(&self.my_to_proj, ShapeType::Edge, ShapeType::Shape) {
            edges.push(e);
        }
        let nb_edges = edges.len();
        for f in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            faces.push(f);
        }
        let nb_faces = faces.len();

        // OCCT L213-216.
        let mut my_res = builder_make_compound();
        vertex_res = builder_make_compound();
        let mut ya_vertex_res = false;

        // Architecture difference #2 (module header): BRepAdaptor_Curve /
        // BRepAdaptor_Surface read BRep_Tool through the owning BRep pool,
        // and a bare rcad Shape does not reference its producing pool.  The
        // input trees are adopted into standalone pools at their own flat
        // indices (topo_builder::brep_from_shape = the TopoDS_Builder::
        // MakeShape / Add glue; the feat loc_ope_find_edges_in_face.rs
        // precedent).  The locations table of the adopting pool is empty:
        // the projection inputs carry identity locations (architecture
        // difference #1 of the feat precedent).
        let brep_edges = rcad_kernel::topo::topo_builder::brep_from_shape(&self.my_to_proj, &[]);
        let brep_faces = rcad_kernel::topo::topo_builder::brep_from_shape(&self.my_shape, &[]);

        // The rcad result pool for the created vertices/edges (the
        // BRepLib_MakeVertex/BRepLib_MakeEdge vehicles; feat precedent).
        let mut pool = rcad_kernel::topods::BRep::new();

        for i in 1..=nb_edges {
            descen_list.clear();
            let edge_shape = edges[i - 1].clone();
            // OCCT L221: hcur = new BRepAdaptor_Curve(TopoDS::Edge(
            // Edges->Value(i))).
            let hcur = BRepAdaptorCurve::with_edge(&brep_edges, &edge_shape);
            // OCCT L222: Elementary = IsElementary(*hcur).
            let elementary = is_elementary(&hcur);
            for j in 1..=nb_faces {
                let face_shape = faces[j - 1].clone();
                // OCCT L225-226: hsur = new BRepAdaptor_Surface(TopoDS::Face(
                // Faces->Value(j))) — the ctor default R = Standard_True (the
                // face UV-window restriction of Initialize).
                let hsur = BRepAdaptorSurface::with_face(&brep_faces, &face_shape, true);

                // OCCT L232-233: TolU = hsur->UResolution(myTol3d) / 20;
                // TolV = hsur->VResolution(myTol3d) / 20.
                let tol_u = hsur.u_resolution(self.my_tol3d) / 20.0;
                let tol_v = hsur.v_resolution(self.my_tol3d) / 20.0;

                // OCCT L238-239: HProjector = new ProjLib_HCompProjectedCurve(
                // hsur, hcur, TolU, TolV, myMaxDist) — the projector over the
                // BRepAdaptor_Surface / BRepAdaptor_Curve handles.
                let h_projector: Arc<CompProjectedCurve> = Arc::new(new_h_comp_projected_curve(
                    Arc::new(hsur.clone()),
                    Arc::new(hcur.clone()),
                    tol_u,
                    tol_v,
                    self.my_max_dist,
                ));

                // OCCT L245-252: the per-solution locals.
                let mut prj = Shape::null(); // OCCT L245: TopoDS_Shape prj (null).
                let mut degenerated = false; // OCCT L246: bool Degenerated = false.
                let mut p2d = glam::DVec2::ZERO;
                let mut pdeb = glam::DVec2::ZERO;
                let mut pfin = glam::DVec2::ZERO;
                let mut uiso = 0.0f64;
                let mut viso = 0.0f64;
                // OCCT L251: occ::handle<Adaptor2d_Curve2d> HPCur — either
                // the isoparametric Geom2dAdaptor_Curve(PCur2d) (L291/L309)
                // or the projector handle itself (L314).
                let mut hp_cur: Curve2dHandle;
                let mut pcur2d: Option<Curve2d> = None; // OCCT L252: only for the isoparametric projection

                for k in 1..=h_projector.nb_curves() {
                    if h_projector.is_single_pnt(k, &mut p2d) {
                        // OCCT L262-268: the punctual solution —
                        // GetSurface()->D0(P2d.X(), P2d.Y(), P).
                        let p = h_projector.get_surface().value(p2d.x, p2d.y);
                        prj = pool.add_tvertex(p);
                        descen_list.push(prj.clone());
                        builder_add_compound_shape(&mut vertex_res, &prj);
                        ya_vertex_res = true;

                        self.my_ancestor_map
                            .insert(shape_key(&prj), (prj.clone(), edge_shape.clone()));
                    } else {
                        // OCCT L272-273.
                        only2d = false;
                        only3d = false;
                        h_projector.bounds(k, &mut udeb, &mut ufin);

                        /**************************************************************/
                        // OCCT L276-315: the isoparametric branches.
                        if h_projector.is_u_iso(k, &mut uiso) {
                            // OCCT L282-283: HProjector->D0(Udeb, Pdeb) /
                            // D0(Ufin, Pfin).
                            pdeb = h_projector.d0(udeb);
                            pfin = h_projector.d0(ufin);
                            poles[0] = pdeb;
                            poles[1] = pfin;
                            knots[0] = udeb;
                            knots[1] = ufin;
                            // OCCT L288-289: BS2d = new Geom2d_BSplineCurve(
                            // Poles, Knots, Mults, Deg).
                            let bs2d = Curve2d::BSpline(BSplineCurve2 {
                                degree: deg,
                                // OCCT knots/mults expanded (Mults == Deg+1).
                                knots: vec![knots[0], knots[0], knots[1], knots[1]],
                                control_points: vec![poles[0], poles[1]],
                                weights: vec![1.0, 1.0],
                            });
                            // OCCT L290: PCur2d = new Geom2d_TrimmedCurve(
                            // BS2d, Udeb, Ufin).
                            pcur2d = Some(Curve2d::Trimmed(TrimmedCurve2 {
                                curve: Box::new(bs2d),
                                t_min: udeb,
                                t_max: ufin,
                            }));
                            // OCCT L291: HPCur = new Geom2dAdaptor_Curve(
                            // PCur2d).
                            hp_cur = geom2d_adaptor_curve(pcur2d.clone().expect("PCur2d"));
                            only3d = true;
                        } else if h_projector.is_v_iso(k, &mut viso) {
                            // OCCT L300-301: HProjector->D0(Udeb, Pdeb) /
                            // D0(Ufin, Pfin).
                            pdeb = h_projector.d0(udeb);
                            pfin = h_projector.d0(ufin);
                            poles[0] = pdeb;
                            poles[1] = pfin;
                            knots[0] = udeb;
                            knots[1] = ufin;
                            // OCCT L306-307: BS2d = new Geom2d_BSplineCurve(
                            // Poles, Knots, Mults, Deg).
                            let bs2d = Curve2d::BSpline(BSplineCurve2 {
                                degree: deg,
                                knots: vec![knots[0], knots[0], knots[1], knots[1]],
                                control_points: vec![poles[0], poles[1]],
                                weights: vec![1.0, 1.0],
                            });
                            // OCCT L308: PCur2d = new Geom2d_TrimmedCurve(
                            // BS2d, Udeb, Ufin).
                            pcur2d = Some(Curve2d::Trimmed(TrimmedCurve2 {
                                curve: Box::new(bs2d),
                                t_min: udeb,
                                t_max: ufin,
                            }));
                            // OCCT L309: HPCur = new Geom2dAdaptor_Curve(
                            // PCur2d).
                            hp_cur = geom2d_adaptor_curve(pcur2d.clone().expect("PCur2d"));
                            only3d = true;
                        } else {
                            // OCCT L314: HPCur = HProjector (the projector
                            // handle upcast to handle(Adaptor2d_Curve2d)).
                            hp_cur = h_projector.clone();
                        }

                        // OCCT L317-320.
                        if (!self.my_with3d || elementary)
                            && h_projector.max_distance(k) <= self.my_tol3d
                        {
                            only2d = true;
                        }

                        if only2d && only3d {
                            // OCCT L324: BRepLib_MakeEdge MKed(
                            // GeomAdaptor::MakeCurve(*hcur), Udeb, Ufin) —
                            // MakeCurve is the identity on the adaptor's
                            // loaded 3D curve (module header architecture
                            // difference #2).  GAP: the curve-on-surface arm
                            // of GeomAdaptor::MakeCurve (an edge without a
                            // 3D curve) is not translated; the OCCT path
                            // derives the 3D image there, rcad skips the
                            // solution.
                            if hcur.transformed.my_con_surf.is_some() {
                                continue;
                            }
                            let mked_curve = hcur.transformed.my_curve.curve.clone();
                            prj = pool.add_tedge(
                                Some(mked_curve),
                                Shape::null(),
                                Shape::null(),
                                [udeb, ufin],
                            );
                            if let Some(pc) = &pcur2d {
                                builder_update_edge_pcurve(&mut prj, pc, &face_shape, self.my_tol3d);
                            }
                            if let Some(mut v) = top_exp_vertices_raw(&prj).0 {
                                builder_update_vertex_tol(&mut v, self.my_tol3d);
                            }
                            if let Some(mut v) = top_exp_vertices_raw(&prj).1 {
                                builder_update_vertex_tol(&mut v, self.my_tol3d);
                            }
                        } else {
                            // OCCT L335-336: Approx_CurveOnSurface appr(
                            // HPCur, hsur, Udeb, Ufin, myTol3d);
                            // appr.Perform(myMaxSeg, myMaxDegree, myContinuity,
                            // Only3d, Only2d).
                            //
                            // Architecture difference: the rcad kernel carries
                            // two GeomAbsShape encodings (topods = the OCCT
                            // 7-variant order, math = the 5-variant adaptor
                            // stack order); myContinuity flows from the
                            // topods form (the set_params API) to the math
                            // form (the approx body).  The G1/G2 arms
                            // reproduce the OCCT normalize step of
                            // Approx_CurveOnSurface.cxx L373-385 (G1 -> C1,
                            // G2 -> C2); the C3/CN -> C2 restriction stays in
                            // the approx body.
                            let appr_continuity = match self.my_continuity {
                                GeomAbsShape::C0 => rcad_kernel::math::GeomAbsShape::C0,
                                GeomAbsShape::G1 => rcad_kernel::math::GeomAbsShape::C1,
                                GeomAbsShape::C1 => rcad_kernel::math::GeomAbsShape::C1,
                                GeomAbsShape::G2 => rcad_kernel::math::GeomAbsShape::C2,
                                GeomAbsShape::C2 => rcad_kernel::math::GeomAbsShape::C2,
                                GeomAbsShape::C3 => rcad_kernel::math::GeomAbsShape::C3,
                                GeomAbsShape::CN => rcad_kernel::math::GeomAbsShape::CN,
                            };
                            let mut appr = ApproxCurveOnSurface::new(
                                hp_cur,
                                Arc::new(hsur.clone()),
                                udeb,
                                ufin,
                                self.my_tol3d,
                            );
                            appr.perform(
                                self.my_max_seg as i32,
                                self.my_max_degree as i32,
                                appr_continuity,
                                only3d,
                                only2d,
                            );

                            // OCCT L338-341.
                            if appr.max_error3d() > 1.0e3 * self.my_tol3d {
                                continue;
                            }

                            // OCCT L357-360.
                            if !only3d {
                                pcur2d = appr.curve2d();
                            }
                            if only2d {
                                // OCCT L361-365: BRepLib_MakeEdge MKed(
                                // GeomAdaptor::MakeCurve(*hcur), Udeb, Ufin) —
                                // see the MakeCurve GAP note above.
                                if hcur.transformed.my_con_surf.is_some() {
                                    continue;
                                }
                                let mked_curve = hcur.transformed.my_curve.curve.clone();
                                prj = pool.add_tedge(
                                    Some(mked_curve),
                                    Shape::null(),
                                    Shape::null(),
                                    [udeb, ufin],
                                );
                            } else {
                                // OCCT L368-439: the degenerated test and
                                // the edge creation.
                                degenerated = true;
                                // OCCT L373: BS3d = appr.Curve3d() — the null
                                // handle would crash the OCCT dereference; the
                                // rcad Option keeps the skip.
                                let Some(bs3d) = appr.curve3d() else {
                                    continue;
                                };
                                let mut p1 = glam::DVec3::ZERO;
                                // OCCT BS3d->FirstParameter()/LastParameter() —
                                // the rcad curve domain.
                                let bs3d_dom = bs3d.default_domain();
                                let bs3d_first = bs3d_dom[0];
                                let bs3d_last = bs3d_dom[1];
                                // start from 3 points to reject non
                                // degenerated edges very fast (OCCT L377-397).
                                let mut nb_point = 3usize;
                                let mut d_par = (bs3d_last - bs3d_first)
                                    / (nb_point - 1) as f64;
                                for ii in 0..nb_point {
                                    let par = bs3d_first + ii as f64 * d_par;
                                    let pp = bs3d.point_at(par);
                                    p1 += pp / nb_point as f64;
                                }
                                for ii in 0..nb_point {
                                    if !degenerated {
                                        break;
                                    }
                                    let par = bs3d_first + ii as f64 * d_par;
                                    let pp = bs3d.point_at(par);
                                    let dist = p1.distance(pp);
                                    if dist > self.my_tol3d {
                                        degenerated = false;
                                        break;
                                    }
                                }
                                // if the test passes a more exact test with
                                // 10 points (OCCT L398-421).
                                if degenerated {
                                    p1 = glam::DVec3::ZERO;
                                    nb_point = 10;
                                    d_par = (bs3d_last - bs3d_first)
                                        / (nb_point - 1) as f64;
                                    for ii in 0..nb_point {
                                        let par = bs3d_first + ii as f64 * d_par;
                                        let pp = bs3d.point_at(par);
                                        p1 += pp / nb_point as f64;
                                    }
                                    for ii in 0..nb_point {
                                        if !degenerated {
                                            break;
                                        }
                                        let par = bs3d_first + ii as f64 * d_par;
                                        let pp = bs3d.point_at(par);
                                        let dist = p1.distance(pp);
                                        if dist > self.my_tol3d {
                                            degenerated = false;
                                            break;
                                        }
                                    }
                                }
                                if degenerated {
                                    // OCCT L422-435.
                                    let mut vv = builder_make_vertex();
                                    builder_update_vertex_point_tol(&mut vv, p1, self.my_tol3d);
                                    prj = builder_make_edge();
                                    builder_add_edge_vertex_fwd_rev(&mut prj, &vv);
                                    builder_set_degenerated(&mut prj, true);
                                } else {
                                    // OCCT L438: prj = BRepLib_MakeEdge(BS3d).Edge().
                                    prj = pool.add_tedge(
                                        Some(Curve3::BSpline(bs3d)),
                                        Shape::null(),
                                        Shape::null(),
                                        [bs3d_first, bs3d_last],
                                    );
                                }
                            }

                            // OCCT L442-451.
                            if let Some(pc) = &pcur2d {
                                let max_err = appr.max_error3d();
                                builder_update_edge_pcurve(&mut prj, pc, &face_shape, max_err);
                                if let Some(v) = top_exp_vertices_raw(&prj).0 {
                                    let mut v = v;
                                    builder_update_vertex_tol(&mut v, max_err);
                                }
                                if let Some(v) = top_exp_vertices_raw(&prj).1 {
                                    let mut v = v;
                                    builder_update_vertex_tol(&mut v, max_err);
                                }
                            }
                            if degenerated {
                                builder_range_edge_on_face(&mut prj, &face_shape, udeb, ufin);
                            }
                        }

                        // OCCT L454-534: the face-bounds trimming and the
                        // result bookkeeping.
                        if self.my_face_bounds {
                            if !degenerated {
                                // OCCT L465: Perform Boolean COMMON
                                // operation to get parts of projected edge
                                // inside the face.
                                let mut a_section =
                                    SectionOp::from_shapes(face_shape.clone(), prj.clone());
                                a_section.build();
                                if a_section.is_done() {
                                    let a_rc = a_section.shape().clone();
                                    for a_e in explorer(&a_rc, ShapeType::Edge, ShapeType::Shape) {
                                        builder_add_compound_shape(&mut my_res, &a_e);
                                        self.my_ancestor_map
                                            .insert(shape_key(&a_e), (a_e.clone(), edge_shape.clone()));
                                        self.my_corresp
                                            .insert(shape_key(&a_e), (a_e.clone(), face_shape.clone()));
                                    }
                                } else {
                                    // OCCT L479-498: if the common operation
                                    // has failed, try to classify the part.
                                    if self.classify_pcur_mid(&face_shape, &pcur2d) {
                                        builder_add_compound_shape(&mut my_res, &prj);
                                        descen_list.push(prj.clone());
                                        self.my_ancestor_map.insert(
                                            shape_key(&prj),
                                            (prj.clone(), edge_shape.clone()),
                                        );
                                        self.my_corresp.insert(
                                            shape_key(&prj),
                                            (prj.clone(), face_shape.clone()),
                                        );
                                    }
                                }
                            } else {
                                // OCCT L500-526: no solution — the
                                // classifier path.
                                if self.classify_pcur_mid(&face_shape, &pcur2d) {
                                    builder_add_compound_shape(&mut my_res, &prj);
                                    descen_list.push(prj.clone());
                                    self.my_ancestor_map
                                        .insert(shape_key(&prj), (prj.clone(), edge_shape.clone()));
                                    self.my_corresp
                                        .insert(shape_key(&prj), (prj.clone(), face_shape.clone()));
                                }
                            }
                        } else {
                            // OCCT L528-534.
                            builder_add_compound_shape(&mut my_res, &prj);
                            descen_list.push(prj.clone());
                            self.my_ancestor_map
                                .insert(shape_key(&prj), (prj.clone(), edge_shape.clone()));
                            self.my_corresp
                                .insert(shape_key(&prj), (prj.clone(), face_shape.clone()));
                        }
                    }
                }
            }
            // OCCT L538.
            self.my_descendants
                .insert(shape_key(&edges[i - 1]), (edges[i - 1].clone(), descen_list.clone()));
        }
        // OCCT L540-547: the eventual wire creation is reported in BuildWire.
        if ya_vertex_res {
            builder_add_compound_shape(&mut my_res, &vertex_res);
        }

        self.my_res = my_res;
        self.my_is_done = true;
    }

    /// The BRepTopAdaptor_FClass2d mid-parameter classification of
    /// Build() (cxx L481-497 / L506-521): state == TopAbs_IN ||
    /// TopAbs_ON.
    fn classify_pcur_mid(&self, face: &Shape, pcur2d: &Option<Curve2d>) -> bool {
        let Some(surf) = brep_tool_surface(face) else {
            return false;
        };
        let Some(pc) = pcur2d else {
            // OCCT would call through the null handle; rcad classifies
            // nothing (marked).
            return false;
        };
        let src = crate::topalgo::shape_source::FaceShapeSource::new(
            face,
            surf,
            &[glam::DAffine3::IDENTITY],
        );
        let classifier = crate::topalgo::brep_top_adaptor::fclass2d::FClass2d::new(
            &src,
            0,
            CONFUSION,
        );
        let f = pc.default_domain()[0];
        let l = pc.default_domain()[1];
        let pmil = (f + l) / 2.0;
        let puv = pc.point_at(pmil);
        let state = classifier.perform(&src, puv, true);
        state == rcad_kernel::topods::State::In || state == rcad_kernel::topods::State::On
    }

    /// OCCT BRepAlgo_NormalProjection::IsDone() const (cxx L588-591).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }

    /// OCCT BRepAlgo_NormalProjection::Projection() const (cxx L595-598).
    pub fn projection(&self) -> &Shape {
        &self.my_res
    }

    /// OCCT BRepAlgo_NormalProjection::Ancestor(E) const (cxx L602-605) —
    /// for a resulting edge, returns the corresponding initial edge
    /// (DataMap::Find throws when unbound).
    pub fn ancestor(&self, e: &Shape) -> &Shape {
        match self.my_ancestor_map.get(&shape_key(e)) {
            Some((_, v)) => v,
            None => panic!("BRepAlgo_NormalProjection::Ancestor"),
        }
    }

    /// OCCT BRepAlgo_NormalProjection::Couple(E) const (cxx L609-612) — for
    /// a projected edge, returns the corresponding initial face.
    pub fn couple(&self, e: &Shape) -> &Shape {
        match self.my_corresp.get(&shape_key(e)) {
            Some((_, v)) => v,
            None => panic!("BRepAlgo_NormalProjection::Couple"),
        }
    }

    /// OCCT BRepAlgo_NormalProjection::Generated(S) (cxx L616-619) — the
    /// list of shapes generated from the shape S.
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        match self.my_descendants.get(&shape_key(s)) {
            Some((_, l)) => l.clone(),
            None => panic!("BRepAlgo_NormalProjection::Generated"),
        }
    }

    /// OCCT BRepAlgo_NormalProjection::IsElementary(C) const
    /// (cxx L623-638) — Line/Circle/Ellipse/Hyperbola/Parabola are
    /// elementary (GetType over the Adaptor3d_Curve; the GeomAbs_CurveType
    /// maps to the kernel base::proj_lib::CurveType — architecture
    /// difference #3).
    pub fn is_elementary(&self, c: &dyn Adaptor3dCurveGeom) -> bool {
        is_elementary(c)
    }

    /// OCCT BRepAlgo_NormalProjection::BuildWire(ListOfWire) const
    /// (cxx L642-677) — builds the result as a list of wire if possible; a
    /// wire is returned only if there is only one wire.
    pub fn build_wire(&self, list_of_wire: &mut Vec<Shape>) -> bool {
        let mut is_wire = false;
        let exp_of_shape = explorer(&self.my_res, ShapeType::Edge, ShapeType::Shape);
        if !exp_of_shape.is_empty() {
            let mut list: Vec<Shape> = Vec::new();

            for cur_e in &exp_of_shape {
                list.push(cur_e.clone());
            }
            // OCCT L656-659: BRepLib_MakeWire MW; MW.Add(List) — pending
            // (architecture difference #4); IsDone() is false so the OCCT
            // not-a-wire exit is taken.
            let mut mw = BRepLibMakeWire::new();
            mw.add(&list);
            if mw.is_done() {
                let wire = mw.shape();
                // If the resulting wire contains the same edge as at the
                // beginning OK, otherwise the result really consists of
                // several wires (OCCT L662-673).
                let nb_edges = explorer(wire, ShapeType::Edge, ShapeType::Shape).len();
                if nb_edges == list.len() {
                    list_of_wire.push(wire.clone());
                    is_wire = true;
                }
            }
        }
        is_wire
    }
}

/// OCCT TopExp::FirstVertex/LastVertex + BRep_Builder::Add(E, V FORWARD /
/// V REVERSED) of the degenerated-edge creation (cxx L431-433).
fn builder_add_edge_vertex_fwd_rev(e: &mut Shape, v: &Shape) {
    let vf = oriented(v, Orientation::Forward);
    crate::brep_algo::tool::builder_add_edge_vertex(e, &vf);
    let vr = oriented(v, Orientation::Reversed);
    crate::brep_algo::tool::builder_add_edge_vertex(e, &vr);
}

/// OCCT BRepAlgo_NormalProjection::IsElementary — shared with the method
/// (cxx L623-638): the GetType switch over the elementary conic kinds.
fn is_elementary(c: &dyn Adaptor3dCurveGeom) -> bool {
    use rcad_kernel::base::proj_lib::CurveType;
    matches!(
        c.get_type(),
        CurveType::Line
            | CurveType::Circle
            | CurveType::Ellipse
            | CurveType::Hyperbola
            | CurveType::Parabola
    )
}

/// OCCT BRepLib_MakeWire — pending TKBRep translation (architecture
/// difference #4).  Add stores the list; IsDone() is false so the callers
/// take the OCCT not-done exit.  GAP: closes with the BRepLib batch.
pub struct BRepLibMakeWire {
    my_list: Vec<Shape>,
}

impl Default for BRepLibMakeWire {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepLibMakeWire {
    pub fn new() -> Self {
        BRepLibMakeWire { my_list: Vec::new() }
    }

    /// OCCT BRepLib_MakeWire::Add(const TopTools_ListOfShape&) — pending.
    pub fn add(&mut self, list: &[Shape]) {
        self.my_list = list.to_vec();
    }

    /// OCCT BRepLib_MakeWire::IsDone() — pending (false).
    pub fn is_done(&self) -> bool {
        false
    }

    /// OCCT BRepLib_MakeWire::Shape() — pending (null).
    pub fn shape(&self) -> &Shape {
        static NULL: std::sync::OnceLock<Shape> = std::sync::OnceLock::new();
        NULL.get_or_init(Shape::null)
    }
}

// OCCT ProjLib_HCompProjectedCurve / Approx_CurveOnSurface: the real 1:1
// bodies live in geomalgo (proj_lib_h_comp_projected_curve{,_b}.rs and
// approx_curve_on_surface.rs) — the former stubs were deleted and every
// call point in Build() is wired to them (the E3-S queue-1 follow-on).
