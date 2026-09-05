//! BRepTopAdaptor_FClass2d (TKTopAlgo) — the 2D face classifier of the
//! BRepTopAdaptor_TopolTool.
//!
//! 1:1 translation of OCCT `BRepTopAdaptor_FClass2d.hxx` (L17-76) + `.cxx`
//! (L17-859): per-wire UV boundary polygons (adaptive sampling with the
//! area/perimeter re-discretization loop), wire orientations, and the
//! classify walk with the periodic recadre and the
//! BRepClass_FaceClassifier fallback.
//!
//! rcad data-model notes:
//! - OCCT TopoDS_Face ↔ kernel `Shape` (face) + `BRep` (the TShape graph).
//! - OCCT BRepAdaptor_Surface(F, false) ↔ `brep.face_surface(face)` (the
//!   face-local surface, consistent with the pcurve frame).
//! - OCCT BRepTools_WireExplorer ↔
//!   `fclass2d::order_wire_edges` over a [`FaceShapeSource`] view of the
//!   face (the same 1:1 WireExplorer translation already landed for the
//!   IntTools_FClass2d port).
//! - OCCT Geom2dInt_Geom2dCurveTool::NbSamples(C) single-arg overload →
//!   [`geom2dint_nb_samples`].

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::topods::{BRep, BRepTool, Orientation, Shape, State, TShape};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::topalgo::shape_source::ShapeSource;
use crate::topalgo::brep_class::face_classifier::FClassifier;
use crate::topalgo::brep_adaptor::curve2d::BRepCurve2d;
use crate::topalgo::brep_top_adaptor::class2d::{Class2d, Class2dResult};
use crate::topalgo::brep_top_adaptor::fclass2d::{order_wire_edges, surface_periodic};
use crate::topalgo::shape_source::FaceShapeSource;
use crate::topalgo::gcpnts::QuasiUniformDeflection;

/// OCCT BRepTopAdaptor_FClass2d.cxx L110: `constexpr double eps = 1.e-10;`
const EPS: f64 = 1.0e-10;

/// OCCT Geom2dInt_Geom2dCurveTool::NbSamples(C) (the single-arg overload,
/// Geom2dInt_Geom2dCurveTool.cxx) — the adaptor NbSamples refined for
/// circles with R > 1 (deflection eps*R, eps = 0.01).
pub(crate) fn geom2dint_nb_samples(c: &BRepCurve2d) -> i32 {
    let mut nbs = Curve2dAdaptor::nb_samples(c);
    if Curve2dAdaptor::get_type(c) == crate::geomalgo::geom2d_int::Curve2dType::Circle {
        let r = Curve2dAdaptor::circle(c).radius;
        if r > 1.0 {
            let angl = 0.283079; // 2.*std::acos(1. - eps);
            let n = ((Curve2dAdaptor::last_parameter(c)
                - Curve2dAdaptor::first_parameter(c))
                / angl) as i32;
            nbs = n.max(nbs);
        }
    }
    nbs
}

/// OCCT anonymous-namespace `safeIncrement` (cxx L54-60) — the next
/// representable value in the direction when + increment is absorbed.
fn safe_increment(the_value: f64, the_direction: f64, the_increment: f64) -> f64 {
    let a_next_value = the_value + the_increment;
    if a_next_value == the_value {
        // std::nextafter(theValue, theDirection).
        if the_direction > the_value {
            the_value.next_up()
        } else {
            the_value.next_down()
        }
    } else {
        a_next_value
    }
}

/// OCCT anonymous-namespace `isDegenerated` (cxx L68-85) — the 3D curve
/// collapses to a point along the parameter span.
fn is_degenerated_3d(brep: &BRep, edge: &Shape, start_param: f64, end_param: f64) -> bool {
    let Some((curve, _range)) = brep.edge_curve_world(edge) else {
        return true;
    };
    let a_parametric_step = (end_param - start_param) * 0.1;
    let a_start_point = curve.point_at(start_param);
    let mut a_curr_param = start_param;
    while a_curr_param < end_param {
        let a_current_point = curve.point_at(a_curr_param);
        if a_start_point.distance_squared(a_current_point)
            > rcad_kernel::precision::CONFUSION
        {
            return false;
        }
        a_curr_param = safe_increment(a_curr_param, end_param, a_parametric_step);
    }
    true
}

/// The 3D curve point of an edge at parameter u (the BRepAdaptor_Curve
/// Value of the FClass2d sampling loop).
fn curve3d_point_at(brep: &BRep, edge: &Shape, u: f64) -> DVec3 {
    match brep.edge_curve_world(edge) {
        Some((curve, _)) => curve.point_at(u),
        // OCCT would dereference a null curve (the No_OutOfRange guards are
        // compiled off); rcad returns the origin for the degenerate path.
        None => DVec3::ZERO,
    }
}

/// OCCT `ElCLib::Parameter(gp_Lin2d, gp_Pnt2d)` for the sag checks.
fn elclib_parameter_2d(a: DVec2, lin_dir: DVec2, p: DVec2) -> f64 {
    (p - a).dot(lin_dir)
}

/// The 3-point sag computation of the sampling loop (cxx L264-281) — the
/// deviation of the middle point from the chord of its neighbours.
fn accumulate_sag(seq: &[DVec2], ii: usize, fleche_u: &mut f64, fleche_v: &mut f64) {
    if ii < 3 {
        return;
    }
    let a = seq[ii - 3];
    let b = seq[ii - 1];
    let mid = seq[ii - 2];
    if a.distance_squared(b) == 0.0 {
        return;
    }
    let chord = b - a;
    let len = chord.length();
    if len <= 1e-15 {
        return;
    }
    let lin_dir = chord / len;
    let ul = elclib_parameter_2d(a, lin_dir, mid);
    let pp = a + ul * lin_dir;
    let d_u = (pp.x - mid.x).abs();
    let d_v = (pp.y - mid.y).abs();
    if d_u > *fleche_u {
        *fleche_u = d_u;
    }
    if d_v > *fleche_v {
        *fleche_v = d_v;
    }
}

/// Build the wire polygon (PClass) and the enclosed-area / perimeter pair
/// (cxx L304-342 and L412-430).
fn build_polygon(seq: &[DVec2], nbpnts: usize) -> (Vec<DVec2>, f64, f64) {
    let mut p_class = vec![DVec2::ZERO; nbpnts];
    let mut square = 0.0f64;
    let mut a_per = 0.0f64;
    if nbpnts < 2 {
        return (p_class, square, a_per);
    }
    let mut im1 = nbpnts - 1;
    let mut im0 = 1usize;
    p_class[nbpnts - 2] = seq[nbpnts - 2];
    p_class[nbpnts - 1] = seq[nbpnts - 1];
    for ii in 1..nbpnts {
        if im1 >= nbpnts {
            im1 = 1;
        }
        p_class[ii - 1] = seq[ii - 1];
        square += (p_class[im0 - 1].x - p_class[im1 - 1].x)
            * (p_class[im0 - 1].y + p_class[im1 - 1].y)
            * 0.5;
        a_per += (p_class[im0 - 1] - p_class[im1 - 1]).length();
        im0 += 1;
        im1 += 1;
    }
    (p_class, square, a_per)
}

/// The default plane for the FaceShapeSource fallback (an inert
/// constructor argument; never evaluated by the classifier when the face
/// carries a surface).
fn fallback_plane_surf() -> Surface3 {
    Surface3::Plane(rcad_kernel::geom::Plane {
        origin: DVec3::ZERO,
        normal: DVec3::Z,
        u_dir: DVec3::X,
        v_dir: DVec3::Y,
    })
}

/// OCCT BRepTopAdaptor_FClass2d — the face classifier.
pub struct FClass2dTopol {
    /// OCCT keeps the shapes alive through TopoDS handles (the refcounted
    /// BRep handle copy — see `brep_adaptor::curve2d`).
    brep: std::sync::Arc<BRep>,
    /// OCCT Face (FORWARD-oriented in the ctor).
    face: Shape,
    /// OCCT TabClass.
    tab_class: Vec<Class2d>,
    /// OCCT TabOrien.
    tab_orien: Vec<i32>,
    /// OCCT Toluv.
    toluv: f64,
    /// OCCT U1 / V1 / U2 / V2 — the periodic recadre windows.
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    /// OCCT Umin / Umax / Vmin / Vmax — the UV bounds of the boundary.
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
}

impl FClass2dTopol {
    /// OCCT BRepTopAdaptor_FClass2d(F, TolUV) (cxx L88-513).
    pub fn new(brep: std::sync::Arc<BRep>, a_face: &Shape, tol_uv: f64) -> Self {
        let mut r = FClass2dTopol {
            brep,
            face: a_face.clone(),
            tab_class: Vec::new(),
            tab_orien: Vec::new(),
            toluv: tol_uv,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            umin: 0.0,
            umax: 0.0,
            vmin: 0.0,
            vmax: 0.0,
        };
        r.init();
        r
    }

    /// The fallback FaceShapeSource view of the face (the
    /// BRepClass_FaceClassifier path).
    fn classifier_source(&self) -> FaceShapeSource<'_> {
        let surf = self
            .brep
            .face_surface_world(&self.face)
            .unwrap_or_else(fallback_plane_surf);
        // The kernel BRep locations table has no identity slot (transforms
        // start at index 0); the ShapeSource pcurve lookup reads the table
        // with the DS convention (slot 0 = identity) — prepend it.
        let locations_ds: Vec<glam::DAffine3> = std::iter::once(glam::DAffine3::IDENTITY)
            .chain(self.brep.locations.iter().copied())
            .collect();
        FaceShapeSource::new(&self.face, surf, &locations_ds)
    }

    /// The ctor body (cxx L96-513).
    fn init(&mut self) {
        //-- dead end on surfaces defined on more than one period
        self.face.orientation = Orientation::Forward;
        // OCCT BRepAdaptor_Surface surf; surf->Initialize(aFace, false) —
        // the face surface (local frame).
        let surf = self.brep.face_surface(&self.face);

        self.umin = 0.0;
        self.vmin = 0.0; // RealLast();
        self.umax = -self.umin;
        self.vmax = -self.vmin;

        let a_face = self.face.clone();
        // OCCT TopExp_Explorer(Face, TopAbs_WIRE).
        let wires: Vec<Shape> = match &*a_face.data {
            TShape::Face(fd) => std::iter::once(&fd.outer_wire)
                .chain(fd.inner_wires.iter())
                .cloned()
                .collect(),
            _ => Vec::new(),
        };

        // The FaceShapeSource view for the WireExplorer (index 0 = face,
        // 1..n = wire edges in traversal order).
        let surf_world = self
            .brep
            .face_surface_world(&a_face)
            .unwrap_or_else(fallback_plane_surf);
        // Kernel -> DS location table convention (slot 0 = identity).
        let locations_ds: Vec<glam::DAffine3> = std::iter::once(glam::DAffine3::IDENTITY)
            .chain(self.brep.locations.iter().copied())
            .collect();
        let source = FaceShapeSource::new(&a_face, surf_world, &locations_ds);

        let mut a_nb_e = 0usize;
        let mut an_is_bad_wire = false;

        for wire in &wires {
            let mut nbpnts = 0usize;
            let mut seq_pnt2d: Vec<DVec2> = Vec::new();
            let mut firstpoint = 1usize;
            let mut fleche_u = 0.0f64;
            let mut fleche_v = 0.0f64;
            let mut wire_is_not_empty = false;
            let mut nb_edges = 0usize;

            // OCCT TopExp_Explorer(wire, TopAbs_EDGE) count.
            if let TShape::Wire(wd) = &*wire.data {
                nb_edges = wd.edges.len();
            }
            a_nb_e = nb_edges;

            let mut ancien_pnt3d = DVec3::ZERO;
            let mut ancien_pnt3d_initialise = false;

            // OCCT BRepTools_WireExplorer(Wire, Face) — the reordered edges.
            let raw_edges: Vec<(usize, Orientation)> = match &*wire.data {
                TShape::Wire(wd) => wd
                    .edges
                    .iter()
                    .map(|e| {
                        (
                            source
                                .map_shape_index(e.ptr_id(), e.location)
                                .unwrap_or(0),
                            e.orientation,
                        )
                    })
                    .collect(),
                _ => Vec::new(),
            };
            let ordered = order_wire_edges(&source, 0, &raw_edges);

            for (ei, or) in &ordered {
                nb_edges = nb_edges.saturating_sub(1);
                let edge = source.shape_at(*ei);
                if *or != Orientation::Forward && *or != Orientation::Reversed {
                    continue;
                }
                // OCCT BRep_Tool::CurveOnSurface(edge, Face, pfbid, plbid).
                let Some((_pcurve, pfbid, plbid)) =
                    self.brep.curve_on_surface(&edge, &a_face)
                else {
                    return;
                };
                if (plbid - pfbid).abs() < 1.0e-9 {
                    continue;
                }

                let mut degenerated = self.brep.is_edge_degenerated(&edge);
                if self.brep.is_edge_closed_on_face(&edge, &a_face) {
                    degenerated = true;
                }
                let va = self.brep.first_vertex(&edge);
                let vb = self.brep.last_vertex(&edge);
                if va.is_null() || vb.is_null() {
                    degenerated = true;
                }

                let a_curve_adaptor_2d =
                    BRepCurve2d::new_edge_face(self.brep.clone(), &edge, &a_face);

                //-- Check cases when it was forgotten to code degenerated:
                //-- PRO17410 (janv 99)
                if !degenerated {
                    degenerated = is_degenerated_3d(&self.brep, &edge, pfbid, plbid);
                }

                //-- ----------------------------------------

                // int nbs = Geom2dInt_Geom2dCurveTool::NbSamples(C).
                let mut nbs = geom2dint_nb_samples(&a_curve_adaptor_2d);
                //-- Attention to rational bsplines of degree 3.
                if nbs > 2 {
                    nbs *= 4;
                }
                let mut du = (plbid - pfbid) / (nbs as f64 - 1.0);
                let mut u;
                if *or == Orientation::Forward {
                    u = pfbid;
                } else {
                    u = plbid;
                    du = -du;
                }

                //-- ------------------------------------------------------
                //-- Try to remote the first point of the current edge from
                //-- the last saved point.
                if firstpoint == 2 {
                    u += du;
                }
                let avant = nbpnts;
                for _e in firstpoint..=nbs as usize {
                    let p2d = Curve2dAdaptor::value(&a_curve_adaptor_2d, u);
                    if p2d.x < self.umin {
                        self.umin = p2d.x;
                    }
                    if p2d.x > self.umax {
                        self.umax = p2d.x;
                    }
                    if p2d.y < self.vmin {
                        self.vmin = p2d.y;
                    }
                    if p2d.y > self.vmax {
                        self.vmax = p2d.y;
                    }

                    let mut dist3dptcourant_ancienpnt = 1e+20;
                    let mut p3d = DVec3::ZERO;
                    if !degenerated {
                        p3d = curve3d_point_at(&self.brep, &edge, u);
                        if nbpnts > 1 && ancien_pnt3d_initialise {
                            dist3dptcourant_ancienpnt = p3d.distance(ancien_pnt3d);
                        }
                    }
                    // patch
                    let mut is_real_curve3d = true;
                    if dist3dptcourant_ancienpnt < rcad_kernel::precision::CONFUSION {
                        let mid_p3d = curve3d_point_at(&self.brep, &edge, u - du / 2.0);
                        if p3d.distance(mid_p3d) < rcad_kernel::precision::CONFUSION {
                            is_real_curve3d = false;
                        }
                    }
                    if is_real_curve3d {
                        if !degenerated {
                            ancien_pnt3d = p3d;
                            ancien_pnt3d_initialise = true;
                        }
                        nbpnts += 1;
                        seq_pnt2d.push(p2d);
                    }

                    u += du;
                    // Modified by Sergey KHROMOV - Fri Apr 19 09:46:12 2002:
                    // the sag of the middle point from the neighbour chord.
                    let ii = nbpnts;
                    if ii > (avant + 4) {
                        accumulate_sag(&seq_pnt2d, ii, &mut fleche_u, &mut fleche_v);
                    }
                } // for(e = firstpoint..=nbs)
                if firstpoint == 1 {
                    firstpoint = 2;
                }
                wire_is_not_empty = true;
            } //-- Edges -> for(WireExplorer)

            if nb_edges != 0 {
                //-- On compte ++ with a normal explorer and with the wire
                //-- explorer: bad wire.
                let p_class = vec![DVec2::ZERO, DVec2::ZERO];
                self.tab_class.push(Class2d::new(
                    &p_class, fleche_u, fleche_v, self.umin, self.vmin, self.umax, self.vmax,
                ));
                an_is_bad_wire = true;
                self.tab_orien.push(-1);
            } else if wire_is_not_empty {
                if nbpnts > 3 {
                    let (mut p_class, mut square, mut a_per) =
                        build_polygon(&seq_pnt2d, nbpnts);

                    let mut an_exp_thick = (2.0 * square.abs() / a_per).max(1e-7);
                    let mut a_defl = fleche_u.max(fleche_v);
                    let mut a_discr_defl = (a_defl * 0.1).min(an_exp_thick * 10.0);
                    while a_defl > an_exp_thick && a_discr_defl > 1e-7 {
                        // Deflection of the polygon is too much for this
                        // ratio of area and perimeter — discretize the wire
                        // more tightly (cxx L347-435).
                        firstpoint = 1;
                        seq_pnt2d.clear();
                        fleche_u = 0.0;
                        fleche_v = 0.0;
                        // OCCT re-runs the WireExplorer over the SAME wire.
                        for (ei2, or2) in &ordered {
                            let edge = source.shape_at(*ei2);
                            if *or2 != Orientation::Forward && *or2 != Orientation::Reversed {
                                continue;
                            }
                            // OCCT BRep_Tool::Range(edge, Face, pfbid, plbid).
                            let Some((c2d, f, l)) =
                                self.brep.curve_on_surface(&edge, &a_face)
                            else {
                                continue;
                            };
                            if (l - f).abs() < 1.0e-9 {
                                continue;
                            }
                            let a_discr =
                                QuasiUniformDeflection::new(&c2d, a_discr_defl, f, l);
                            if !a_discr.is_done() {
                                break;
                            }
                            let nbp = a_discr.nb_points() as i32;
                            let (mut i_step, mut i, mut i_end): (i32, i32, i32) =
                                (1, 1, nbp + 1);
                            if *or2 == Orientation::Reversed {
                                i_step = -1;
                                i = nbp;
                                i_end = 0;
                            }
                            if firstpoint == 2 {
                                i += i_step;
                            }
                            while i != i_end {
                                let a_p2d = rcad_kernel::geom::Curve2dEval::point_at(&c2d, a_discr.parameter(i as usize));
                                seq_pnt2d.push(a_p2d);
                                i += i_step;
                            }
                            if nbp > 2 {
                                let ii = seq_pnt2d.len();
                                accumulate_sag(&seq_pnt2d, ii, &mut fleche_u, &mut fleche_v);
                            }
                            firstpoint = 2;
                        }
                        nbpnts = seq_pnt2d.len();
                        p_class.clear();
                        p_class.resize(nbpnts, DVec2::ZERO);
                        let rebuilt = build_polygon(&seq_pnt2d, nbpnts);
                        p_class = rebuilt.0;
                        square = rebuilt.1;
                        a_per = rebuilt.2;

                        an_exp_thick = (2.0 * square.abs() / a_per).max(1e-7);
                        a_defl = fleche_u.max(fleche_v);
                        a_discr_defl = (a_discr_defl * 0.1).min(an_exp_thick * 10.0);
                    }

                    //-- FlecheU *= 10; FlecheV *= 10;
                    // OCCT: TabOrien.Append((square < 0.0) ? 1 : 0).
                    if a_nb_e == 1 && fleche_u < EPS && fleche_v < EPS && square.abs() < EPS {
                        self.tab_orien.push(1);
                    } else {
                        self.tab_orien.push(if square < 0.0 { 1 } else { 0 });
                    }

                    if fleche_u < self.toluv {
                        fleche_u = self.toluv;
                    }
                    if fleche_v < self.toluv {
                        fleche_v = self.toluv;
                    }
                    self.tab_class.push(Class2d::new(
                        &p_class,
                        fleche_u,
                        fleche_v,
                        self.umin,
                        self.vmin,
                        self.umax,
                        self.vmax,
                    ));
                } // if(nbpoints > 3)
                else {
                    an_is_bad_wire = true;
                    self.tab_orien.push(-1);
                    let x_p_class = vec![seq_pnt2d[0], seq_pnt2d[1]];
                    self.tab_class.push(Class2d::new(
                        &x_p_class,
                        fleche_u,
                        fleche_v,
                        self.umin,
                        self.vmin,
                        self.umax,
                        self.vmax,
                    ));
                }
            } // else if(WireIsNotEmpty)

            if an_is_bad_wire {
                break;
            }
        } // for(FaceExplorer)

        let nbtabclass = self.tab_class.len();
        if nbtabclass > 0 {
            //-- If an error was detected on a wire: set all TabOrien to -1.
            if an_is_bad_wire {
                self.tab_orien[0] = -1;
            }

            let (typ_u_per, typ_v_per) = surf.as_ref().map(|s| surface_periodic(s)).unwrap_or((false, false));
            let is_torus = surf
                .as_ref()
                .map(|s| matches!(s, Surface3::Torus(_)))
                .unwrap_or(false);

            if typ_u_per {
                let mut uuu = std::f64::consts::PI + std::f64::consts::PI - (self.umax - self.umin);
                if uuu < 0.0 {
                    uuu = 0.0;
                }
                let _ = uuu;
                self.u1 = 0.0; // modified by NIZHNY-OFV
                self.u2 = 2.0 * std::f64::consts::PI;
            } else {
                self.u1 = 0.0;
                self.u2 = 0.0;
            }

            if is_torus {
                // The V compensation checks the torus V periodicity — OCCT
                // keys it on the Torus type; rcad mirrors with the V flag.
                let _ = typ_v_per;
                self.v1 = 0.0; // modified by NIZHNY-OFV
                self.v2 = 2.0 * std::f64::consts::PI;
            } else {
                self.v1 = 0.0;
                self.v2 = 0.0;
            }
        }
    }

    /// OCCT PerformInfinitePoint() (cxx L515-523).
    pub fn perform_infinite_point(&self) -> State {
        if self.umax == -f64::MAX
            || self.vmax == -f64::MAX
            || self.umin == f64::MAX
            || self.vmin == f64::MAX
        {
            return State::In;
        }
        let p = DVec2::new(
            self.umin - (self.umax - self.umin),
            self.vmin - (self.vmax - self.vmin),
        );
        self.perform(p, false)
    }

    /// The classify core shared by Perform / TestOnRestriction (cxx
    /// L590-683 / L752-841): the TabClass walk with the recadre loop.
    fn classify_walk(
        &self,
        _puv: DVec2,
        tol: f64,
        recadre_on_periodic: bool,
        on_mode: bool,
    ) -> State {
        let mut dedans;
        let nbtabclass = self.tab_class.len();

        if nbtabclass == 0 {
            return State::In;
        }

        //-- U1 is the First Param and U2 in this case is U1 + Period.
        let mut u = _puv.x;
        let mut v = _puv.y;
        let mut uu = u;
        let mut vv = v;

        let surf = self.brep.face_surface(&self.face);
        let (is_u_per, is_v_per) = surf
            .as_ref()
            .map(|s| surface_periodic(s))
            .unwrap_or((false, false));
        // OCCT uperiod/vperiod from BRepAdaptor_Surface (2PI for the
        // analytic periodic surfaces rcad models).
        let uperiod = if is_u_per {
            std::f64::consts::TAU
        } else {
            0.0
        };
        let vperiod = if is_v_per {
            std::f64::consts::TAU
        } else {
            0.0
        };
        let mut a_status;
        let mut urecadre = false;
        let mut vrecadre = false;

        if recadre_on_periodic {
            if is_u_per {
                if uu < self.umin {
                    while uu < self.umin {
                        uu += uperiod;
                    }
                } else {
                    while uu >= self.umin {
                        uu -= uperiod;
                    }
                    uu += uperiod;
                }
            }
            if is_v_per {
                if vv < self.vmin {
                    while vv < self.vmin {
                        vv += vperiod;
                    }
                } else {
                    while vv >= self.vmin {
                        vv -= vperiod;
                    }
                    vv += vperiod;
                }
            }
        }

        loop {
            dedans = 1;
            let puv = DVec2::new(u, v);

            if self.tab_orien[0] != -1 {
                for n in 0..nbtabclass {
                    let cur = if on_mode {
                        self.tab_class[n].si_dans_on_mode(puv, tol)
                    } else {
                        self.tab_class[n].si_dans(puv)
                    };
                    match cur {
                        Class2dResult::Inside => {
                            if self.tab_orien[n] == 0 {
                                dedans = -1;
                                break;
                            }
                        }
                        Class2dResult::Outside => {
                            if self.tab_orien[n] == 1 {
                                dedans = -1;
                                break;
                            }
                        }
                        Class2dResult::Uncertain => {
                            dedans = 0;
                            break;
                        }
                    }
                }
                if dedans == 0 {
                    if on_mode {
                        a_status = State::On;
                    } else {
                        // BRepClass_FaceClassifier fallback; the OCCT
                        // m_Toluv clamp (Perform cxx L625).
                        let mut a_classifier = FClassifier::new();
                        let m_toluv = self.toluv.min(4.0);
                        let source = self.classifier_source();
                        a_classifier.perform(&source, 0, puv, m_toluv);
                        a_status = a_classifier.state();
                    }
                } else if dedans == 1 {
                    a_status = State::In;
                } else {
                    a_status = State::Out;
                }
            } else {
                //-- TabOrien(1) = -1: false wire.
                let mut a_classifier = FClassifier::new();
                let source = self.classifier_source();
                a_classifier.perform(&source, 0, puv, if on_mode { tol } else { self.toluv });
                a_status = a_classifier.state();
            }

            if !recadre_on_periodic || (!is_u_per && !is_v_per) {
                return a_status;
            }
            if a_status == State::In || a_status == State::On {
                return a_status;
            }

            if !urecadre {
                u = uu;
                urecadre = true;
            } else if is_u_per {
                u += uperiod;
            }
            if u > self.umax || !is_u_per {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if is_v_per {
                    v += vperiod;
                }

                u = uu;

                if v > self.vmax || !is_v_per {
                    return a_status;
                }
            }
        }
    }

    /// OCCT Perform(Puv, RecadreOnPeriodic) (cxx L525-684).
    pub fn perform(&self, puv: DVec2, recadre_on_periodic: bool) -> State {
        self.classify_walk(puv, self.toluv, recadre_on_periodic, false)
    }

    /// OCCT TestOnRestriction(Puv, Tol, RecadreOnPeriodic) (cxx L686-842).
    pub fn test_on_restriction(&self, puv: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        self.classify_walk(puv, tol, recadre_on_periodic, true)
    }

    /// OCCT Destroy() (cxx L844-847) — the TabClass clear; Rust drop.
    pub fn destroy(&mut self) {
        self.tab_class.clear();
    }

    #[cfg(test)]
    pub(crate) fn debug_state(&self) -> (usize, Vec<i32>, Vec<(i32, f64, f64, f64, f64, Vec<(DVec2, String)>)>) {
        (
            self.tab_class.len(),
            self.tab_orien.clone(),
            self.tab_class
                .iter()
                .map(|c| {
                    let probe = |x: f64, y: f64| match c.si_dans(DVec2::new(x, y)) {
                        Class2dResult::Inside => "+1".to_string(),
                        Class2dResult::Outside => "-1".to_string(),
                        Class2dResult::Uncertain => " 0".to_string(),
                    };
                    (
                        c.si_dans(DVec2::new(0.0, 0.0)) as i32, // placeholder for count
                        0.0,
                        0.0,
                        0.0,
                        0.0,
                        vec![
                            (DVec2::new(0.5, 0.5), probe(0.5, 0.5)),
                            (DVec2::new(2.0, 2.0), probe(2.0, 2.0)),
                            (DVec2::new(0.5, 0.0), probe(0.5, 0.0)),
                        ],
                    )
                })
                .collect(),
        )
    }
}

#[allow(unused_imports)]
use Curve2d as _;
#[allow(unused_imports)]
use SurfaceEval as _;
