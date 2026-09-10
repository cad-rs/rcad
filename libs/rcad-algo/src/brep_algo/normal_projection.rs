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
//! 2. BRepAdaptor_Curve/BRepAdaptor_Surface map to the rcad geometry values
//!    carried by the Shape (Curve3 / Surface3); GeomAdaptor::MakeCurve
//!    (cxx L324/L363) is the identity (the rcad adaptor IS the geom value).
//! 3. GeomAbs_Shape -> rcad_kernel::topods::GeomAbsShape; GeomAbs_CurveType
//!    maps to the Curve3 variant discriminant (IsElementary).
//! 4. The projection machinery is GAP-pending:
//!    - ProjLib_HCompProjectedCurve (TKGeomBase/TKTopAlgo projection of a
//!      curve on a surface) — pending; nb_curves() is 0 so the OCCT
//!      per-solution loop never runs (the pending-classification behavior).
//!      GAP: closes with the ProjLib batch.
//!    - Approx_CurveOnSurface (TKGeomBase approximation) — pending; kept
//!      behind the same GAP.
//!    - BRepLib_MakeWire (TKBRep connected-wire builder) — pending;
//!      IsDone() is false so BuildWire takes the OCCT not-a-wire exit.
//!      GAP: closes with the BRepLib batch.
//! 5. BRepAlgoAPI_Section (cxx L465) -> bop::brep_algo_api::SectionOp (the
//!    real BOPAlgo_Section-driven implementation).
//! 6. BRepTopAdaptor_FClass2d -> topalgo::brep_top_adaptor::fclass2d::
//!    FClass2d over FaceShapeSource (loc_ope_wires_on_shape_b.rs #8).
//! 7. BRep_Builder edits are in-place Arc::make_mut mutations (tool.rs);
//!    OCCT TShape sharing makes them visible through every handle, rcad
//!    mutates the owning copy.

use crate::brep_algo::tool::{
    brep_tool_curve, brep_tool_surface, builder_add_compound_shape, builder_make_compound,
    builder_make_edge, builder_make_vertex, builder_range_edge_on_face,
    builder_set_degenerated, builder_update_edge_pcurve, builder_update_vertex_point_tol,
    builder_update_vertex_tol, explorer, oriented, shape_key, top_exp_vertices_raw, ShapeKey,
};
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval, TrimmedCurve2,
    BSplineCurve2,
};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{GeomAbsShape, Orientation, ShapeType};
use std::collections::HashMap;

use crate::bop::brep_algo_api::{Algo, SectionOp};

/// OCCT handle<Adaptor2d_Curve2d> HPCur (cxx L251) — either the
/// isoparametric Geom2dAdaptor_Curve(PCur2d) (cxx L291/L309) or the
/// ProjLib_HCompProjectedCurve handle itself (cxx L314).
pub enum HPCurve {
    /// Geom2dAdaptor_Curve over the trimmed isoparametric pcurve.
    Geom2d,
    /// The ProjLib_HCompProjectedCurve handle.
    Projector,
}

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

        // The rcad result pool for the created vertices/edges (the
        // BRepLib_MakeVertex/BRepLib_MakeEdge vehicles; feat precedent).
        let mut pool = rcad_kernel::topods::BRep::new();

        for i in 1..=nb_edges {
            descen_list.clear();
            let edge_shape = edges[i - 1].clone();
            // OCCT L221: hcur = new BRepAdaptor_Curve(TopoDS::Edge(...)).
            let hcur = brep_tool_curve(&edge_shape);
            let elementary = match &hcur {
                Some((c, _, _)) => is_elementary(c),
                None => false,
            };
            for j in 1..=nb_faces {
                let face_shape = faces[j - 1].clone();
                // OCCT L225-226: hsur = new BRepAdaptor_Surface(...).
                let Some(hsur) = brep_tool_surface(&face_shape) else {
                    continue;
                };

                // computation of TolU and TolV (OCCT L230-233).
                let tol_u = rcad_kernel::topo::topods::u_resolution_for_surface(&hsur, self.my_tol3d) / 20.0;
                let tol_v = rcad_kernel::topo::topods::v_resolution_for_surface(&hsur, self.my_tol3d) / 20.0;

                // OCCT L238-239: the projection tool (architecture
                // difference #4 — pending ProjLib; nb_curves() is 0 so the
                // per-solution loop below never runs).
                let h_projector = ProjLibHCompProjectedCurve::new(
                    &hsur,
                    hcur.as_ref().map(|(c, _, _)| c),
                    tol_u,
                    tol_v,
                    self.my_max_dist,
                );

                // OCCT L245-252: the per-solution locals.
                let mut prj = Shape::null(); // OCCT L245: TopoDS_Shape prj (null).
                let mut degenerated = false; // OCCT L246: bool Degenerated = false.
                let mut p2d = glam::DVec2::ZERO;
                let mut pdeb = glam::DVec2::ZERO;
                let mut pfin = glam::DVec2::ZERO;
                let mut uiso = 0.0f64;
                let mut viso = 0.0f64;
                let mut hp_cur = HPCurve::Projector;
                let mut pcur2d: Option<Curve2d> = None; // Only for isoparametric projection

                for k in 1..=h_projector.nb_curves() {
                    if h_projector.is_single_pnt(k, &mut p2d) {
                        // OCCT L262-268: the punctual solution.
                        let p = h_projector.get_surface().point_at(p2d.x, p2d.y);
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
                            h_projector.d0(udeb, &mut pdeb);
                            h_projector.d0(ufin, &mut pfin);
                            poles[0] = pdeb;
                            poles[1] = pfin;
                            knots[0] = udeb;
                            knots[1] = ufin;
                            let bs2d = Curve2d::BSpline(BSplineCurve2 {
                                degree: deg,
                                // OCCT knots/mults expanded (Mults == Deg+1).
                                knots: vec![knots[0], knots[0], knots[1], knots[1]],
                                control_points: vec![poles[0], poles[1]],
                                weights: vec![1.0, 1.0],
                            });
                            pcur2d = Some(Curve2d::Trimmed(TrimmedCurve2 {
                                curve: Box::new(bs2d),
                                t_min: udeb,
                                t_max: ufin,
                            }));
                            hp_cur = HPCurve::Geom2d;
                            only3d = true;
                        } else if h_projector.is_v_iso(k, &mut viso) {
                            h_projector.d0(udeb, &mut pdeb);
                            h_projector.d0(ufin, &mut pfin);
                            poles[0] = pdeb;
                            poles[1] = pfin;
                            knots[0] = udeb;
                            knots[1] = ufin;
                            let bs2d = Curve2d::BSpline(BSplineCurve2 {
                                degree: deg,
                                knots: vec![knots[0], knots[0], knots[1], knots[1]],
                                control_points: vec![poles[0], poles[1]],
                                weights: vec![1.0, 1.0],
                            });
                            pcur2d = Some(Curve2d::Trimmed(TrimmedCurve2 {
                                curve: Box::new(bs2d),
                                t_min: udeb,
                                t_max: ufin,
                            }));
                            hp_cur = HPCurve::Geom2d;
                            only3d = true;
                        } else {
                            // OCCT L314: HPCur = HProjector.
                            hp_cur = HPCurve::Projector;
                        }

                        // OCCT L317-320.
                        if (!self.my_with3d || elementary)
                            && h_projector.max_distance(k) <= self.my_tol3d
                        {
                            only2d = true;
                        }

                        if only2d && only3d {
                            // OCCT L322-329.
                            let Some((c, _, _)) = &hcur else {
                                continue;
                            };
                            prj = pool.add_tedge(Some(c.clone()), Shape::null(), Shape::null(), [udeb, ufin]);
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
                            // OCCT L331-341: the approximation
                            // (architecture difference #4 — pending).
                            let mut appr = ApproxCurveOnSurface::new(
                                &hp_cur,
                                &hsur,
                                udeb,
                                ufin,
                                self.my_tol3d,
                            );
                            appr.perform(
                                self.my_max_seg,
                                self.my_max_degree,
                                self.my_continuity,
                                only3d,
                                only2d,
                            );

                            if appr.max_error_3d() > 1.0e3 * self.my_tol3d {
                                continue;
                            }

                            // OCCT L357-360.
                            if !only3d {
                                pcur2d = appr.curve_2d();
                            }
                            if only2d {
                                // OCCT L361-365.
                                let Some((c, _, _)) = &hcur else {
                                    continue;
                                };
                                prj = pool.add_tedge(Some(c.clone()), Shape::null(), Shape::null(), [udeb, ufin]);
                            } else {
                                // OCCT L368-439: the degenerated test and
                                // the edge creation.
                                degenerated = true;
                                let Some(bs3d) = appr.curve_3d() else {
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
                                    prj = pool.add_tedge(Some(bs3d), Shape::null(), Shape::null(), [
                                        bs3d_first,
                                        bs3d_last,
                                    ]);
                                }
                            }

                            // OCCT L442-451.
                            if let Some(pc) = &pcur2d {
                                let max_err = appr.max_error_3d();
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
    /// elementary (the GeomAbs_CurveType maps to the rcad Curve3 variant of
    /// the basis curve; architecture difference #3).
    pub fn is_elementary(&self, c: &Curve3) -> bool {
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

/// OCCT BRepAlgo_NormalProjection::IsElementary — shared with the method.
fn is_elementary(c: &Curve3) -> bool {
    use rcad_kernel::geom::Curve3 as C3;
    // OCCT BRepAdaptor_Curve::GetType: the trimmed wrapper is transparent
    // (the basis curve type is reported).
    let basis = match c {
        C3::Trimmed(t) => t.basis_curve(),
        other => other,
    };
    matches!(
        basis,
        C3::Line(_) | C3::Circle(_) | C3::Ellipse(_) | C3::Hyperbola(_) | C3::Parabola(_)
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

/// OCCT ProjLib_HCompProjectedCurve — pending TKGeomBase/TKTopAlgo
/// translation (architecture difference #4).  The accessor surface mirrors
/// the OCCT handle; nb_curves() is 0 so the Build() per-solution loop never
/// runs (the pending-classification behavior).  GAP: closes with the ProjLib
/// batch.
// The stored constructor state (myCurve/myTolU/myTolV/myMaxDist) is the
// pending implementation working set; dead-code is allowed while the GAP
// stands.
#[allow(dead_code)]
pub struct ProjLibHCompProjectedCurve {
    my_surface: Option<Surface3>, // OCCT: GetSurface()
    my_curve: Option<Curve3>,
    my_tol_u: f64,
    my_tol_v: f64,
    my_max_dist: f64,
}

impl ProjLibHCompProjectedCurve {
    /// OCCT new ProjLib_HCompProjectedCurve(S, C, TolU, TolV, MaxDist).
    pub fn new(
        the_surface: &Surface3,
        the_curve: Option<&Curve3>,
        the_tol_u: f64,
        the_tol_v: f64,
        the_max_dist: f64,
    ) -> Self {
        ProjLibHCompProjectedCurve {
            my_surface: Some(the_surface.clone()),
            my_curve: the_curve.cloned(),
            my_tol_u: the_tol_u,
            my_tol_v: the_tol_v,
            my_max_dist: the_max_dist,
        }
    }

    /// OCCT ProjLib_HCompProjectedCurve::NbCurves() — pending (0).
    pub fn nb_curves(&self) -> usize {
        0
    }

    /// OCCT IsSinglePnt(Index, P) — pending (false).
    pub fn is_single_pnt(&self, _index: usize, _p: &mut glam::DVec2) -> bool {
        false
    }

    /// OCCT Bounds(Index, U1, U2) — pending.
    pub fn bounds(&self, _index: usize, _u1: &mut f64, _u2: &mut f64) {}

    /// OCCT IsUIso(Index, U) — pending (false).
    pub fn is_u_iso(&self, _index: usize, _u: &mut f64) -> bool {
        false
    }

    /// OCCT IsVIso(Index, V) — pending (false).
    pub fn is_v_iso(&self, _index: usize, _v: &mut f64) -> bool {
        false
    }

    /// OCCT D0(U, P) — pending.
    pub fn d0(&self, _u: f64, _p: &mut glam::DVec2) {}

    /// OCCT MaxDistance(Index) — pending (0).
    pub fn max_distance(&self, _index: usize) -> f64 {
        0.0
    }

    /// OCCT GetSurface() — the projected-on surface.
    pub fn get_surface(&self) -> &Surface3 {
        self.my_surface.as_ref().expect("GetSurface()")
    }
}

/// OCCT Approx_CurveOnSurface — pending TKGeomBase translation
/// (architecture difference #4).  Perform leaves the errors at 0 and the
/// curves null (the OCCT failure output).  GAP: closes with the Approx
/// batch.
#[allow(dead_code)]
pub struct ApproxCurveOnSurface {
    my_hp_cur: HPCurve, // OCCT: the handle<Adaptor2d_Curve2d>
    my_surf: Surface3,  // OCCT: the handle<Adaptor3d_Surface>
    my_udeb: f64,
    my_ufin: f64,
    my_tol3d: f64,
    my_max_error_3d: f64,
    my_curve2d: Option<Curve2d>,
    my_curve3d: Option<Curve3>,
}

impl ApproxCurveOnSurface {
    /// OCCT Approx_CurveOnSurface(HPC, SURF, Udeb, Ufin, Tol3d).
    pub fn new(
        the_hp_cur: &HPCurve,
        the_surf: &Surface3,
        the_udeb: f64,
        the_ufin: f64,
        the_tol3d: f64,
    ) -> Self {
        ApproxCurveOnSurface {
            my_hp_cur: match the_hp_cur {
                HPCurve::Geom2d => HPCurve::Geom2d,
                HPCurve::Projector => HPCurve::Projector,
            },
            my_surf: the_surf.clone(),
            my_udeb: the_udeb,
            my_ufin: the_ufin,
            my_tol3d: the_tol3d,
            my_max_error_3d: 0.0,
            my_curve2d: None,
            my_curve3d: None,
        }
    }

    /// OCCT Approx_CurveOnSurface::Perform(MaxSeg, MaxDegree, Continuity,
    /// Only3d, Only2d) — pending (the failure output).
    pub fn perform(
        &mut self,
        _max_seg: usize,
        _max_degree: usize,
        _continuity: GeomAbsShape,
        _only3d: bool,
        _only2d: bool,
    ) {
        self.my_max_error_3d = 0.0;
        self.my_curve2d = None;
        self.my_curve3d = None;
    }

    /// OCCT MaxError3d().
    pub fn max_error_3d(&self) -> f64 {
        self.my_max_error_3d
    }

    /// OCCT Curve2d().
    pub fn curve_2d(&self) -> Option<Curve2d> {
        self.my_curve2d.clone()
    }

    /// OCCT Curve3d().
    pub fn curve_3d(&self) -> Option<Curve3> {
        self.my_curve3d.clone()
    }
}
