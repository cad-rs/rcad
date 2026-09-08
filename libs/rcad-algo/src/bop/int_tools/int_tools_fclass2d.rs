//! OCCT IntTools_FClass2d (TKBO/IntTools) — strict 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBO/IntTools/IntTools_FClass2d.hxx
//! (L37-91) + IntTools_FClass2d.cxx (L1-951).
//!
//! Function map (OCCT -> rcad):
//!   IntTools_FClass2d()                     cxx L57       -> Default / new()
//!   IntTools_FClass2d(aFace, TolUV)         cxx L61-66    -> new_face(brep, aFace, TolUV)
//!   IsHole()                                cxx L70-73    -> is_hole
//!   Init(aFace, TolUV)                      cxx L77-621   -> init
//!   PerformInfinitePoint()                  cxx L625-633  -> perform_infinite_point
//!   Perform(_Puv, RecadreOnPeriodic)        cxx L637-804  -> perform
//!   TestOnRestriction(_Puv, Tol, Recadre)   cxx L808-943  -> test_on_restriction
//!   ~IntTools_FClass2d()                    cxx L947-950  -> Drop
//!
//! rcad data-model notes (architecture differences):
//! - OCCT TopoDS_Face carries the TShape graph through handles; rcad stores
//!   the owning `BRep` explicitly (`brep: Arc<BRep>`) next to the face Shape.
//! - OCCT BRepTools_WireExplorer (Init/More/Next/Current) <-> the landed 1:1
//!   translation `fclass2d::order_wire_edges` over a `FaceShapeSource` view
//!   of the face.
//! - OCCT BRepAdaptor_Curve C3d (edge 3D curve) <-> `BRep::edge_curve_world`
//!   (edge Location applied); a missing 3D curve keeps the `degenerated`
//!   state instead of raising on C3d.Value.
//! - OCCT caches `myFExplorer` across Perform calls; the rcad
//!   FClassifier::perform takes the DS view and drives its own explorer
//!   internally, so the cached `FaceExplorer` mirrors the OCCT lazy-init
//!   form only.
//! - The OCCT `#ifdef DEBUG_PCLASS_POLYGON` block (cxx L439-452) is a
//!   compile-time debug aid and is not translated.

use std::cell::RefCell;
use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2dEval, Curve3, CurveEval};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, State, TShape};
use rcad_kernel::{CONFUSION, SQUARE_CONFUSION};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::topalgo::brep_adaptor::curve2d::BRepCurve2d;
use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;
use crate::topalgo::brep_class::face_classifier::FClassifier;
use crate::topalgo::brep_class::face_explorer::FaceExplorer;
use crate::topalgo::brep_top_adaptor::class2d::{Class2d, Class2dResult};
use crate::topalgo::brep_top_adaptor::fclass2d::order_wire_edges;
use crate::topalgo::brep_top_adaptor::fclass2d_topol::geom2dint_nb_samples;
use crate::topalgo::gcpnts::QuasiUniformDeflection;
use crate::topalgo::shape_source::{FaceShapeSource, ShapeSource};

/// OCCT RealLast() for double.
const REAL_LAST: f64 = f64::MAX;

/// OCCT GeomInt::AdjustPeriodic (GeomInt.cxx L21-48) — translate a parameter
/// by whole periods so it lands inside [theParMin, theParMax]. Returns
/// `(theNewPar, theOffset)`.
fn adjust_periodic(
    the_par: f64,
    the_par_min: f64,
    the_par_max: f64,
    the_period: f64,
) -> (f64, f64) {
    let the_offset;
    let mut the_new_par = the_par;
    let b_min = the_par_min - the_par > 0.0;
    let b_max = the_par - the_par_max > 0.0;
    if b_min || b_max {
        let dp = if b_min {
            the_par_max - the_par
        } else {
            the_par_min - the_par
        };
        let a_nb_per = (dp / the_period).trunc(); // modf() integer part
        the_offset = a_nb_per * the_period;
        the_new_par += the_offset;
    } else {
        the_offset = 0.0;
    }
    (the_new_par, the_offset)
}

/// OCCT Poly::PolygonProperties (Poly.hxx L165-196) — signed area and
/// perimeter of a 2D polygon. Area is negative when bypassed clockwise.
fn polygon_properties(pts: &[DVec2]) -> (f64, f64) {
    let n = pts.len();
    if n < 2 {
        return (0.0, 0.0);
    }
    let a_ref_pnt = pts[0];
    let mut a_prev_pt = pts[1] - a_ref_pnt;
    let mut the_area = 0.0;
    let mut the_perimeter = a_prev_pt.length();
    for i in 2..n {
        let a_curr_pt = pts[i] - a_ref_pnt;
        let a_delta = a_prev_pt.x * a_curr_pt.y - a_prev_pt.y * a_curr_pt.x; // Crossed
        the_area += a_delta;
        the_perimeter += (a_prev_pt - a_curr_pt).length();
        a_prev_pt = a_curr_pt;
    }
    the_perimeter += a_prev_pt.length();
    the_area *= 0.5;
    (the_area, the_perimeter)
}

/// OCCT ElCLib::Parameter(gp_Lin2d, gp_Pnt2d) for a line `(origin, unit dir)`.
fn elclib_parameter_lin2d(origin: DVec2, dir: DVec2, p: DVec2) -> f64 {
    (p - origin).dot(dir)
}

/// OCCT ElCLib::Value(ul, gp_Lin2d).
fn elclib_value_lin2d(ul: f64, origin: DVec2, dir: DVec2) -> DVec2 {
    origin + dir * ul
}

/// OCCT IntTools_FClass2d (IntTools_FClass2d.hxx L37-91) — classify a 2d
/// point in 2d space of face using boundaries of the face.
pub struct IntToolsFClass2d {
    /// rcad: the owning BRep (OCCT reaches it through the TopoDS handles).
    brep: std::sync::Arc<BRep>,
    /// OCCT TabClass — per-wire CSLib_Class2d classifiers.
    tab_class: Vec<Class2d>,
    /// OCCT TabOrien — per-wire orientation (1 / 0 / -1).
    tab_orien: Vec<i32>,
    /// OCCT Toluv — UV tolerance.
    toluv: f64,
    /// OCCT Face (FORWARD-oriented by Init).
    face: Shape,
    /// OCCT U1 / V1 / U2 / V2 — periodic recadre windows.
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
    /// OCCT Umin / Umax / Vmin / Vmax — UV bounds of the sampled boundary.
    umin: f64,
    umax: f64,
    vmin: f64,
    vmax: f64,
    /// OCCT myIsHole.
    my_is_hole: bool,
    /// OCCT myFExplorer (mutable unique_ptr<BRepClass_FaceExplorer>) —
    /// created lazily in Perform / TestOnRestriction.
    my_f_explorer: RefCell<Option<FaceExplorer>>,
}

// OCCT IntTools_FClass2d::IntTools_FClass2d() — cxx L57 (`= default;`).
impl Default for IntToolsFClass2d {
    fn default() -> Self {
        IntToolsFClass2d {
            brep: std::sync::Arc::new(BRep::new()),
            tab_class: Vec::new(),
            tab_orien: Vec::new(),
            toluv: 0.0,
            face: Shape::null(),
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            umin: 0.0,
            umax: 0.0,
            vmin: 0.0,
            vmax: 0.0,
            my_is_hole: true,
            my_f_explorer: RefCell::new(None),
        }
    }
}

impl IntToolsFClass2d {
    /// OCCT IntTools_FClass2d::IntTools_FClass2d() — cxx L57.
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT IntTools_FClass2d::IntTools_FClass2d(aFace, TolUV) — cxx L61-66:
    /// `Toluv(TolUV), Face(aFace) { Init(Face, Toluv); }`.
    pub fn new_face(brep: std::sync::Arc<BRep>, a_face: &Shape, tol_uv: f64) -> Self {
        let mut c = IntToolsFClass2d {
            brep,
            tab_class: Vec::new(),
            tab_orien: Vec::new(),
            toluv: tol_uv,
            face: a_face.clone(),
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            umin: 0.0,
            umax: 0.0,
            vmin: 0.0,
            vmax: 0.0,
            my_is_hole: true,
            my_f_explorer: RefCell::new(None),
        };
        // OCCT L65: Init(Face, Toluv).
        let toluv = c.toluv;
        let face = c.face.clone();
        c.init(&face, toluv);
        c
    }

    /// OCCT IntTools_FClass2d::IsHole — cxx L70-73.
    pub fn is_hole(&self) -> bool {
        self.my_is_hole
    }

    /// The FaceShapeSource view of the face (the BRepClass_FaceExplorer /
    /// BRepTools_WireExplorer DS adapter). Face index 0 = the face itself.
    fn classifier_source(&self) -> FaceShapeSource<'_> {
        let surf = self.brep.face_surface_world(&self.face).unwrap_or_else(|| {
            rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            })
        });
        // Kernel -> DS location table convention (slot 0 = identity).
        let locations_ds: Vec<glam::DAffine3> = std::iter::once(glam::DAffine3::IDENTITY)
            .chain(self.brep.locations.iter().copied())
            .collect();
        FaceShapeSource::new(&self.face, surf, &locations_ds)
    }

    /// OCCT IntTools_FClass2d::Init — IntTools_FClass2d.cxx L77-621.
    pub fn init(&mut self, a_face: &Shape, tol_uv: f64) {
        // OCCT L79-97: local declarations.
        let mut firstpoint: i32;
        let mut nb_edges: i32;
        let mut i_x: i32;
        let mut a_nbs1: i32;
        let mut nbs: i32;
        let mut avant: i32;
        let mut bad_wire: i32;
        let mut u: f64;
        let mut du: f64;
        let mut tole: f64;
        let mut a_tol: f64;
        let mut fleche_u: f64;
        let mut fleche_v: f64;
        let mut tol_vertex1: f64;
        let mut tol_vertex: f64;
        let mut u_first: f64;
        let mut u_last: f64;
        let a_pr_cf: f64;
        let a_pr_cf2: f64;
        //
        let mut edge: Shape;
        let mut va: Shape;
        let mut vb: Shape;
        let mut or_: Orientation;
        let mut ancien_pnt3d: DVec3;
        let mut seq_pnt2d: Vec<DVec2>;
        let mut an_index_map: HashMap<i32, i32>;
        let mut a_d1_prev: Vec<DVec2>;
        let mut a_d1_next: Vec<DVec2>;
        //
        // OCCT L99-101.
        a_pr_cf = CONFUSION;
        a_pr_cf2 = a_pr_cf * a_pr_cf;
        self.my_is_hole = true;
        //
        // OCCT L103-105.
        self.toluv = tol_uv;
        self.face = a_face.clone();
        self.face.orientation = Orientation::Forward;
        // OCCT L106-107: surf = new BRepAdaptor_Surface(); Initialize(aFace, false).
        let adaptor_brep = self.brep.clone();
        let surf = BRepAdaptorSurface::initialize_face(&adaptor_brep, &self.face, false);
        //
        // OCCT L109-113.
        tole = 0.0;
        a_tol = 0.0;
        self.umin = REAL_LAST;
        self.vmin = REAL_LAST;
        self.umax = -self.umin;
        self.vmax = -self.vmin;
        bad_wire = 0;
        //
        // if face has several wires and one of them is bad,
        // it is necessary to process all of them for correct
        // calculation of Umin, Umax, Vmin, Vmax - ifv, 23.08.06
        //
        // OCCT L119: aExpF.Init(Face, TopAbs_WIRE) — the face wires (the
        // outer wire first, then the inner wires).
        let wires: Vec<Shape> = match &*self.face.data {
            TShape::Face(fd) => std::iter::once(&fd.outer_wire)
                .chain(fd.inner_wires.iter())
                .cloned()
                .collect(),
            _ => Vec::new(),
        };
        // OCCT L120: for (; aExpF.More(); aExpF.Next()).
        for a_w in &wires {
            //
            // OCCT L124-136.
            firstpoint = 1;
            fleche_u = 0.0;
            fleche_v = 0.0;
            tol_vertex1 = 0.0;
            tol_vertex = 0.0;
            let mut wire_is_not_empty = false;
            let mut ancien_pnt3d_initialise = false;
            ancien_pnt3d = DVec3::ZERO; // SetCoord(0., 0., 0.)
            //
            seq_pnt2d = Vec::new();
            an_index_map = HashMap::new();
            a_d1_prev = Vec::new();
            a_d1_next = Vec::new();
            //
            // NbEdges — OCCT L139-144: TopExp_Explorer(aW, TopAbs_EDGE) count.
            nb_edges = 0;
            if let TShape::Wire(wd) = &*a_w.data {
                nb_edges = wd.edges.len() as i32;
            }
            //
            // OCCT L146: aWExp.Init(aW, Face) — BRepTools_WireExplorer over
            // the FaceShapeSource view of the face (face index 0).
            let surf_world = self.brep.face_surface_world(&self.face).unwrap_or_else(|| {
                rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane {
                    origin: DVec3::ZERO,
                    normal: DVec3::Z,
                    u_dir: DVec3::X,
                    v_dir: DVec3::Y,
                })
            });
            let locations_ds: Vec<glam::DAffine3> = std::iter::once(glam::DAffine3::IDENTITY)
                .chain(self.brep.locations.iter().copied())
                .collect();
            let source = FaceShapeSource::new(&self.face, surf_world, &locations_ds);
            let raw_edges: Vec<(usize, Orientation)> = match &*a_w.data {
                TShape::Wire(wd) => wd
                    .edges
                    .iter()
                    .map(|e| {
                        (
                            source.map_shape_index(e.ptr_id(), e.location).unwrap_or(0),
                            e.orientation,
                        )
                    })
                    .collect(),
                _ => Vec::new(),
            };
            let ordered = order_wire_edges(&source, 0, &raw_edges);
            // OCCT L147: for (; aWExp.More(); aWExp.Next()).
            for (ei, edge_ori) in &ordered {
                nb_edges -= 1;
                // The oriented edge occurrence (WireExplorer::Current); the
                // seam pcurve selection keys on the orientation.
                edge = source.shape_at(*ei);
                edge.orientation = *edge_ori;
                or_ = edge.orientation;
                // OCCT L152-155.
                if or_ != Orientation::Forward && or_ != Orientation::Reversed {
                    continue;
                }
                //
                // OCCT L157-161: aC2D = BRep_Tool::CurveOnSurface(edge, Face,
                // pfbid, plbid); if (aC2D.IsNull()) return;
                let (pfbid, plbid) = match self.brep.curve_on_surface(&edge, &self.face) {
                    Some((_c2d, f, l)) => (f, l),
                    None => return,
                };
                //
                // OCCT L163-164: BRepAdaptor_Curve2d C(edge, Face);
                // BRepAdaptor_Curve C3d;
                let c = BRepCurve2d::new_edge_face(self.brep.clone(), &edge, &self.face);
                let mut c3d: Option<Curve3> = None;
                //------------------------------------------
                // OCCT L166-170.
                let mut degenerated = self.brep.is_edge_degenerated(&edge)
                    || self.brep.is_edge_closed_on_face(&edge, &self.face);
                //
                // OCCT L172: TopExp::Vertices(edge, Va, Vb).
                va = self.brep.first_vertex(&edge);
                vb = self.brep.last_vertex(&edge);
                //
                // OCCT L174-196.
                tol_vertex1 = 0.0;
                tol_vertex = 0.0;
                if va.is_null() {
                    degenerated = true;
                } else {
                    tol_vertex1 = self.brep.vertex_tolerance(&va);
                }
                if vb.is_null() {
                    degenerated = true;
                } else {
                    tol_vertex = self.brep.vertex_tolerance(&vb);
                }
                //
                if tol_vertex < tol_vertex1 {
                    tol_vertex = tol_vertex1;
                }
                //
                //-- Verification of cases when forgotten to code degenereted
                // OCCT L199-220: check that whole curve is located in vicinity
                // of its middle point (within sphere of Precision::Confusion()
                // diameter). rcad: a missing 3D curve keeps the state (OCCT
                // would raise on C3d.Value of an unloaded adaptor).
                if !degenerated {
                    c3d = self.brep.edge_curve_world(&edge).map(|(cu, _)| cu);
                    if let Some(cu) = &c3d {
                        let p3da = cu.point_at(0.5 * (pfbid + plbid));
                        du = plbid - pfbid;
                        const NBSTEPS: i32 = 10;
                        let a_prec2 = 0.25 * CONFUSION * CONFUSION;
                        let mut still_degenerated = true;
                        for i in 0..=NBSTEPS {
                            let uu = pfbid + i as f64 * du / NBSTEPS as f64;
                            let p3db = cu.point_at(uu);
                            let a_r2 = p3da.distance_squared(p3db);
                            if a_r2 > a_prec2 {
                                still_degenerated = false;
                                break;
                            }
                        }
                        degenerated = still_degenerated;
                    }
                } // if(!degenerated)
                //-- ----------------------------------------
                // OCCT L222-226.
                tole = self.brep.tolerance(&edge);
                if tole > a_tol {
                    a_tol = tole;
                }
                //
                // NbSamples +> nbs — OCCT L229-233.
                nbs = geom2dint_nb_samples(&c);
                if nbs > 2 {
                    nbs *= 4;
                }
                // OCCT L234.
                du = (plbid - pfbid) / (nbs - 1) as f64;
                //
                // OCCT L236-248.
                if or_ == Orientation::Forward {
                    u = pfbid;
                    u_first = pfbid;
                    u_last = plbid;
                } else {
                    u = plbid;
                    u_first = plbid;
                    u_last = pfbid;
                    du = -du;
                }
                //
                // aPrms — OCCT L251-270.
                a_nbs1 = nbs + 1;
                let mut a_prms: Vec<f64> = vec![0.0; (a_nbs1 + 1) as usize]; // 1-based filler at 0
                if nbs == 2 {
                    let a_coef = 0.0025;
                    a_prms[1] = u_first;
                    a_prms[2] = u_first + a_coef * (u_last - u_first);
                    a_prms[3] = u_last;
                } else if nbs > 2 {
                    a_nbs1 = nbs;
                    a_prms[1] = u_first;
                    i_x = 2;
                    while i_x < a_nbs1 {
                        a_prms[i_x as usize] = u + (i_x - 1) as f64 * du;
                        i_x += 1;
                    }
                    a_prms[a_nbs1 as usize] = u_last;
                } else {
                    // rcad: OCCT leaves aPrms uninitialized for nbs < 2 (the
                    // entries are never filled); Rust requires initialized
                    // memory — fill with uFirst (unreachable in practice:
                    // nbs >= 2 by construction).
                    a_nbs1 = 1;
                    a_prms[1] = u_first;
                }
                //
                //-- ------------------------------------------------------------
                //-- Check distance uv between the start point of the edge
                //-- and the last point saved in SeqPnt2d
                //-- To to set the first point of the current
                //-- afar from the last saved point
                // OCCT L277.
                avant = seq_pnt2d.len() as i32;
                // OCCT L278: for (iX = firstpoint; iX <= aNbs1; iX++).
                i_x = firstpoint;
                while i_x <= a_nbs1 {
                    //
                    // OCCT L280-284.
                    let ii: usize;
                    let mut a_dst_x: f64;
                    let p2d: DVec2;
                    let p3d: Option<DVec3>;
                    //
                    // OCCT L286-303.
                    u = a_prms[i_x as usize];
                    p2d = c.value(u);
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
                    //
                    // OCCT L305-316.
                    a_dst_x = REAL_LAST;
                    p3d = if !degenerated {
                        c3d.as_ref().map(|cu| cu.point_at(u))
                    } else {
                        None
                    };
                    if let Some(p3) = &p3d {
                        if !seq_pnt2d.is_empty() {
                            if ancien_pnt3d_initialise {
                                a_dst_x = p3.distance_squared(ancien_pnt3d);
                            }
                        }
                    }
                    //
                    // OCCT L318-333.
                    let mut is_real_curve3d = true;
                    if a_dst_x < a_pr_cf2 {
                        if i_x > 1 {
                            if let Some(p3) = &p3d {
                                if let Some(cu) = &c3d {
                                    let a_dst_x1: f64;
                                    let mid_p3d =
                                        cu.point_at(0.5 * (u + a_prms[(i_x - 1) as usize]));
                                    a_dst_x1 = p3.distance_squared(mid_p3d);
                                    if a_dst_x1 < a_pr_cf2 {
                                        is_real_curve3d = false;
                                    }
                                }
                            }
                        }
                    }
                    //
                    // OCCT L335-343.
                    if is_real_curve3d {
                        if !degenerated {
                            if let Some(p3) = &p3d {
                                ancien_pnt3d = *p3;
                                ancien_pnt3d_initialise = true;
                            }
                        }
                        seq_pnt2d.push(p2d);
                    }
                    //
                    // OCCT L345-364.
                    ii = seq_pnt2d.len();
                    if ii as i32 > (avant + 4) {
                        // gp_Lin2d Lin(SeqPnt2d(ii-2), Dir2d(Vec2d(SeqPnt2d(ii-2), SeqPnt2d(ii))))
                        let a = seq_pnt2d[ii - 3];
                        let b = seq_pnt2d[ii - 1];
                        let chord = b - a;
                        let len = chord.length();
                        // rcad: OCCT gp_Dir2d raises on the null vector; skip.
                        if len > 0.0 {
                            let lin_dir = chord / len;
                            let ul = elclib_parameter_lin2d(a, lin_dir, seq_pnt2d[ii - 2]);
                            let pp = elclib_value_lin2d(ul, a, lin_dir);
                            let d_u = (pp.x - seq_pnt2d[ii - 2].x).abs();
                            let d_v = (pp.y - seq_pnt2d[ii - 2].y).abs();
                            if d_u > fleche_u {
                                fleche_u = d_u;
                            }
                            if d_v > fleche_v {
                                fleche_v = d_v;
                            }
                        }
                    }
                    i_x += 1;
                } // for(iX=firstpoint; iX<=aNbs1; iX++) {
                //
                // OCCT L367-372.
                if bad_wire != 0 {
                    continue; // if face has several wires and one of them is bad,
                              // it is necessary to process all of them for correct
                              // calculation of Umin, Umax, Vmin, Vmax - ifv, 23.08.06
                }
                //
                // OCCT L374-377.
                if firstpoint == 1 {
                    firstpoint = 2;
                }
                // OCCT L378.
                wire_is_not_empty = true;
                // Append the derivative of the first parameter. (OCCT L379-391)
                let a_u = a_prms[1];
                let (_a_p, mut a_v) = c.d1(a_u);
                if or_ == Orientation::Reversed {
                    a_v = -a_v; // gp_Vec2d::Reverse
                }
                a_d1_next.push(a_v);
                // Append the derivative of the last parameter. (OCCT L393-409)
                let a_u = a_prms[a_nbs1 as usize];
                let (_a_p, mut a_v) = c.d1(a_u);
                if or_ == Orientation::Reversed {
                    a_v = -a_v;
                }
                if nb_edges > 0 {
                    a_d1_prev.push(a_v);
                } else {
                    a_d1_prev.insert(0, a_v);
                }
                // Fill the map anIndexMap. (OCCT L411-419)
                if avant > 0 {
                    an_index_map.insert(avant, a_d1_next.len() as i32);
                } else {
                    an_index_map.insert(1, a_d1_next.len() as i32);
                }
            } // for(;aWExp.More(); aWExp.Next()) {
            //
            // OCCT L423-433.
            if nb_edges != 0 {
                //-- count ++ with normal explorer and -- with Wire Explorer
                // OCCT L426-429: NCollection_Array1<gp_Pnt2d> PClass(1, 2);
                // PClass.Init(anInitPnt(0., 0.)).
                let p_class = vec![DVec2::ZERO, DVec2::ZERO];
                self.tab_class.push(Class2d::new(
                    &p_class,
                    fleche_u,
                    fleche_v,
                    self.umin,
                    self.vmin,
                    self.umax,
                    self.vmax,
                ));
                bad_wire = 1;
                self.tab_orien.push(-1);
            }
            // OCCT L435-575: else if (WireIsNotEmpty).
            else if wire_is_not_empty {
                // OCCT L437: if (SeqPnt2d.Length() > 3).
                if seq_pnt2d.len() > 3 {
                    // OCCT L439-452: #ifdef DEBUG_PCLASS_POLYGON debug draw —
                    // a compile-time debug aid, not translated.
                    //
                    // OCCT L454-456: aS = 0.; aPer = 0.;
                    // Poly::PolygonProperties(SeqPnt2d, aS, aPer).
                    let mut a_s = 0.0;
                    let mut a_per = 0.0;
                    let (s, p) = polygon_properties(&seq_pnt2d);
                    a_s = s;
                    a_per = p;
                    //
                    // OCCT L458-461.
                    let mut an_exp_thick = (2.0 * a_s.abs() / a_per).max(1e-7);
                    let mut a_defl = fleche_u.max(fleche_v);
                    let mut a_discr_defl = (a_defl * 0.1).min(an_exp_thick * 10.0);
                    let mut is_changed = false;
                    // OCCT L462: while (aDefl > anExpThick && aDiscrDefl > 1e-7).
                    while a_defl > an_exp_thick && a_discr_defl > 1e-7 {
                        // Deflection of the polygon is too much for this ratio of area and perimeter,
                        // and this might lead to self-intersections.
                        // Discretize the wire more tightly to eliminate the error.
                        // OCCT L467-471.
                        firstpoint = 1;
                        is_changed = true;
                        seq_pnt2d.clear();
                        fleche_u = 0.0;
                        fleche_v = 0.0;
                        // OCCT L472: aWExp.Init(TopoDS::Wire(aExpF.Current()),
                        // Face) — the same wire as the first walk (rcad: the
                        // reordered edge list).
                        for (ei2, ori2) in &ordered {
                            edge = source.shape_at(*ei2);
                            edge.orientation = *ori2;
                            or_ = edge.orientation;
                            // OCCT L476.
                            if or_ == Orientation::Forward || or_ == Orientation::Reversed {
                                // OCCT L478: BRep_Tool::Range(edge, Face, pfbid, plbid).
                                let (pfbid, plbid) =
                                    match self.brep.curve_on_surface(&edge, &self.face) {
                                        Some((_pc, f, l)) => (f, l),
                                        // rcad: Range is always present here (the
                                        // pcurve existed in the first walk); skip.
                                        None => continue,
                                    };
                                // OCCT L479-482.
                                if (plbid - pfbid).abs() < 1.0e-9 {
                                    continue;
                                }
                                // OCCT L483: BRepAdaptor_Curve2d C(edge, Face).
                                let (pc, _f, _l) =
                                    match self.brep.curve_on_surface(&edge, &self.face) {
                                        Some(x) => x,
                                        None => continue,
                                    };
                                // OCCT L484.
                                let a_discr =
                                    QuasiUniformDeflection::new(&pc, a_discr_defl, pfbid, plbid);
                                // OCCT L485-488.
                                if !a_discr.is_done() {
                                    break;
                                }
                                let nbp = a_discr.nb_points() as i32;
                                let (mut i_step, mut i, mut i_end): (i32, i32, i32) =
                                    (1, 1, nbp + 1);
                                // OCCT L491-496.
                                if or_ == Orientation::Reversed {
                                    i_step = -1;
                                    i = nbp;
                                    i_end = 0;
                                }
                                // OCCT L497-500.
                                if firstpoint == 2 {
                                    i += i_step;
                                }
                                // OCCT L501-505.
                                while i != i_end {
                                    let a_p2d = pc.point_at(a_discr.parameter(i as usize));
                                    seq_pnt2d.push(a_p2d);
                                    i += i_step;
                                }
                                // OCCT L506-522.
                                if nbp > 2 {
                                    let ii = seq_pnt2d.len();
                                    let a = seq_pnt2d[ii - 3];
                                    let b = seq_pnt2d[ii - 1];
                                    let chord = b - a;
                                    let len = chord.length();
                                    // rcad: OCCT gp_Dir2d raises on the null vector; skip.
                                    if len > 0.0 {
                                        let lin_dir = chord / len;
                                        let ul =
                                            elclib_parameter_lin2d(a, lin_dir, seq_pnt2d[ii - 2]);
                                        let pp = elclib_value_lin2d(ul, a, lin_dir);
                                        let d_u = (pp.x - seq_pnt2d[ii - 2].x).abs();
                                        let d_v = (pp.y - seq_pnt2d[ii - 2].y).abs();
                                        if d_u > fleche_u {
                                            fleche_u = d_u;
                                        }
                                        if d_v > fleche_v {
                                            fleche_v = d_v;
                                        }
                                    }
                                }
                                // OCCT L523.
                                firstpoint = 2;
                            }
                        }
                        // OCCT L526-528.
                        an_exp_thick = (2.0 * a_s.abs() / a_per).max(1e-7);
                        a_defl = fleche_u.max(fleche_v);
                        a_discr_defl = (a_discr_defl * 0.1).min(an_exp_thick * 10.0);
                    }
                    //
                    // OCCT L531-534.
                    if is_changed {
                        let (s, _p) = polygon_properties(&seq_pnt2d);
                        a_s = s;
                    }
                    //
                    // OCCT L536-544.
                    if fleche_u < self.toluv {
                        fleche_u = self.toluv;
                    }
                    if fleche_v < self.toluv {
                        fleche_v = self.toluv;
                    }
                    // OCCT L546.
                    self.tab_class.push(Class2d::new(
                        &seq_pnt2d,
                        fleche_u,
                        fleche_v,
                        self.umin,
                        self.vmin,
                        self.umax,
                        self.vmax,
                    ));
                    //
                    // OCCT L548-565.
                    if a_s.abs() < SQUARE_CONFUSION {
                        bad_wire = 1;
                        self.tab_orien.push(-1);
                    } else {
                        if a_s > 0.0 {
                            self.my_is_hole = false;
                            self.tab_orien.push(1);
                        } else {
                            self.my_is_hole = true;
                            self.tab_orien.push(0);
                        }
                    }
                }
                // OCCT L567-574: else.
                else {
                    bad_wire = 1;
                    self.tab_orien.push(-1);
                    // OCCT L571: NCollection_Array1<gp_Pnt2d> PPClass(1, 2) —
                    // created but unused.
                    seq_pnt2d.clear();
                    // OCCT L573.
                    self.tab_class.push(Class2d::new(
                        &seq_pnt2d,
                        fleche_u,
                        fleche_v,
                        self.umin,
                        self.vmin,
                        self.umax,
                        self.vmax,
                    ));
                }
            } // else if(WireIsNotEmpty)
        } // for(; aExpF.More();  aExpF.Next()) {
        //
        // OCCT L578.
        let nbtabclass = self.tab_class.len();
        //
        // OCCT L580-620.
        if nbtabclass > 0 {
            //-- if an error on a wire was detected : all TabOrien set to -1
            // OCCT L583-586.
            if bad_wire != 0 {
                self.tab_orien[0] = -1;
            }
            // OCCT L588-599.
            let st = surf.get_type();
            if st == GeomAbsSurfaceType::Cone
                || st == GeomAbsSurfaceType::Cylinder
                || st == GeomAbsSurfaceType::Torus
                || st == GeomAbsSurfaceType::Sphere
                || st == GeomAbsSurfaceType::SurfaceOfRevolution
            {
                let mut uuu = std::f64::consts::PI + std::f64::consts::PI - (self.umax - self.umin);
                if uuu < 0.0 {
                    uuu = 0.0;
                }
                self.u1 = self.umin - uuu * 0.5;
                self.u2 = self.u1 + std::f64::consts::PI + std::f64::consts::PI;
            } else {
                self.u1 = 0.0;
                self.u2 = 0.0;
            }
            // OCCT L605-619.
            if st == GeomAbsSurfaceType::Torus {
                let mut uuu = std::f64::consts::PI + std::f64::consts::PI - (self.vmax - self.vmin);
                if uuu < 0.0 {
                    uuu = 0.0;
                }
                self.v1 = self.vmin - uuu * 0.5;
                self.v2 = self.v1 + std::f64::consts::PI + std::f64::consts::PI;
            } else {
                self.v1 = 0.0;
                self.v2 = 0.0;
            }
        }
    }

    /// OCCT IntTools_FClass2d::PerformInfinitePoint — cxx L625-633.
    pub fn perform_infinite_point(&self) -> State {
        if self.umax == -REAL_LAST
            || self.vmax == -REAL_LAST
            || self.umin == REAL_LAST
            || self.vmin == REAL_LAST
        {
            return State::In;
        }
        let p = DVec2::new(
            self.umin - (self.umax - self.umin),
            self.vmin - (self.vmax - self.vmin),
        );
        self.perform(p, false)
    }

    /// OCCT IntTools_FClass2d::Perform — cxx L637-804.
    pub fn perform(&self, _puv: DVec2, recadre_on_periodic: bool) -> State {
        // OCCT L639-643.
        let nbtabclass = self.tab_class.len();
        if nbtabclass == 0 {
            return State::In;
        }
        //
        //-- U1 is the First Param and U2 is in this case U1+Period
        // OCCT L646-650.
        let mut u = _puv.x;
        let mut v = _puv.y;
        let mut uu = u;
        let mut vv = v;
        let mut a_status = State::Unknown;
        //
        // OCCT L652-653: surf = new BRepAdaptor_Surface(); Initialize(Face, false).
        let adaptor_brep = self.brep.clone();
        let surf = BRepAdaptorSurface::initialize_face(&adaptor_brep, &self.face, false);
        //
        // OCCT L655-658.
        let is_u_per = surf.is_u_periodic();
        let is_v_per = surf.is_v_periodic();
        let uperiod = if is_u_per { surf.u_period() } else { 0.0 };
        let vperiod = if is_v_per { surf.v_period() } else { 0.0 };
        //
        // OCCT L660-664.
        let mut urecadre: bool;
        let mut vrecadre: bool;
        let mut b_use_classifier: bool;
        let mut dedans: i32 = 1;
        //
        urecadre = false;
        vrecadre = false;
        //
        // OCCT L666-678.
        if recadre_on_periodic {
            if is_u_per {
                let (new_uu, _du) = adjust_periodic(uu, self.umin, self.umax, uperiod);
                uu = new_uu;
            } // if (IsUPer) {
            //
            if is_v_per {
                let (new_vv, _dv) = adjust_periodic(vv, self.vmin, self.vmax, vperiod);
                vv = new_vv;
            } // if (IsVPer) {
        }
        //
        // OCCT L680: for (;;).
        loop {
            // OCCT L682-683.
            dedans = 1;
            let puv = DVec2::new(u, v);
            // OCCT L684.
            b_use_classifier = self.tab_orien[0] == -1;
            // OCCT L685-724.
            if !b_use_classifier {
                let mut n: usize = 1;
                let mut tab_orien_n: i32;
                let mut cur: Class2dResult;
                while n <= nbtabclass {
                    cur = self.tab_class[n - 1].si_dans(puv);
                    tab_orien_n = self.tab_orien[n - 1];
                    //
                    if cur == Class2dResult::Inside {
                        if tab_orien_n == 0 {
                            dedans = -1;
                            break;
                        }
                    } else if cur == Class2dResult::Outside {
                        if tab_orien_n == 1 {
                            dedans = -1;
                            break;
                        }
                    } else {
                        dedans = 0;
                        break;
                    }
                    n += 1;
                } // for(n=1; n<=nbtabclass; n++)
                //
                // OCCT L716-723.
                if dedans == 0 {
                    b_use_classifier = true;
                } else {
                    a_status = if dedans == 1 { State::In } else { State::Out };
                }
            } // if(TabOrien(1)!=-1) {
            //
            // compute state of the point using face classifier (OCCT L726-756)
            if b_use_classifier {
                // compute tolerance to use in face classifier (OCCT L728-745)
                let a_u_res: f64;
                let a_v_res: f64;
                let a_fc_tol: f64;
                let b_u_in: bool;
                let b_v_in: bool;
                //
                a_u_res = surf.u_resolution(self.toluv);
                a_v_res = surf.v_resolution(self.toluv);
                //
                b_u_in = u >= self.umin && u <= self.umax;
                b_v_in = v >= self.vmin && v <= self.vmax;
                //
                if b_u_in == b_v_in {
                    a_fc_tol = a_u_res.min(a_v_res);
                } else {
                    a_fc_tol = if !b_u_in { a_u_res } else { a_v_res };
                }
                //
                // OCCT L748-755: the lazily created face explorer, then
                // BRepClass_FClassifier::Perform(*myFExplorer, Puv, aFCTol).
                // rcad FClassifier::perform drives its own explorer over the
                // DS view (architecture note).
                {
                    let source = self.classifier_source();
                    if self.my_f_explorer.borrow().is_none() {
                        *self.my_f_explorer.borrow_mut() = Some(FaceExplorer::new(&source, 0));
                    }
                    let mut a_classifier = FClassifier::new();
                    a_classifier.perform(&source, 0, puv, a_fc_tol);
                    a_status = a_classifier.state();
                }
            }
            //
            // OCCT L758-761.
            if !recadre_on_periodic || (!is_u_per && !is_v_per) {
                return a_status;
            }
            // OCCT L763-766.
            if a_status == State::In || a_status == State::On {
                return a_status;
            }
            //
            // OCCT L768-779.
            if !urecadre {
                u = uu;
                urecadre = true;
            } else {
                if is_u_per {
                    u += uperiod;
                }
            }
            //
            // OCCT L781-802.
            if u > self.umax || !is_u_per {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else {
                    if is_v_per {
                        v += vperiod;
                    }
                }
                //
                u = uu;
                //
                if v > self.vmax || !is_v_per {
                    return a_status;
                }
            }
        } // while (1)
    }

    /// OCCT IntTools_FClass2d::TestOnRestriction — cxx L808-943.
    pub fn test_on_restriction(&self, _puv: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        // OCCT L812-816.
        let nbtabclass = self.tab_class.len();
        if nbtabclass == 0 {
            return State::In;
        }
        //
        //-- U1 is the First Param and U2 in this case is U1+Period
        // OCCT L819-821.
        let mut u = _puv.x;
        let mut v = _puv.y;
        let mut uu = u;
        let mut vv = v;
        //
        // OCCT L823-828.
        let adaptor_brep = self.brep.clone();
        let surf = BRepAdaptorSurface::initialize_face(&adaptor_brep, &self.face, false);
        let is_u_per = surf.is_u_periodic();
        let is_v_per = surf.is_v_periodic();
        let uperiod = if is_u_per { surf.u_period() } else { 0.0 };
        let vperiod = if is_v_per { surf.v_period() } else { 0.0 };
        // OCCT L829-831.
        let mut a_status = State::Unknown;
        let mut urecadre = false;
        let mut vrecadre = false;
        let mut dedans: i32 = 1;
        //
        // OCCT L833-845.
        if recadre_on_periodic {
            if is_u_per {
                let (new_uu, _du) = adjust_periodic(uu, self.umin, self.umax, uperiod);
                uu = new_uu;
            } // if (IsUPer) {
            //
            if is_v_per {
                let (new_vv, _dv) = adjust_periodic(vv, self.vmin, self.vmax, vperiod);
                vv = new_vv;
            } // if (IsVPer) {
        }
        //
        // OCCT L847: for (;;).
        loop {
            // OCCT L849-850.
            dedans = 1;
            let puv = DVec2::new(u, v);
            //
            // OCCT L852-891.
            if self.tab_orien[0] != -1 {
                let mut n: i32 = 1;
                while n <= nbtabclass as i32 {
                    let cur = self.tab_class[(n - 1) as usize].si_dans_on_mode(puv, tol);
                    if cur == Class2dResult::Inside {
                        if self.tab_orien[(n - 1) as usize] == 0 {
                            dedans = -1;
                            break;
                        }
                    } else if cur == Class2dResult::Outside {
                        if self.tab_orien[(n - 1) as usize] == 1 {
                            dedans = -1;
                            break;
                        }
                    } else {
                        dedans = 0;
                        break;
                    }
                    n += 1;
                }
                // OCCT L879-890.
                if dedans == 0 {
                    a_status = State::On;
                }
                if dedans == 1 {
                    a_status = State::In;
                }
                if dedans == -1 {
                    a_status = State::Out;
                }
            }
            // OCCT L892-903: else — TabOrien(1) = -1, wrong wire.
            else {
                {
                    let source = self.classifier_source();
                    if self.my_f_explorer.borrow().is_none() {
                        *self.my_f_explorer.borrow_mut() = Some(FaceExplorer::new(&source, 0));
                    }
                    // OCCT L900-902: aClassifier.Perform(*myFExplorer, Puv, Tol).
                    let mut a_classifier = FClassifier::new();
                    a_classifier.perform(&source, 0, puv, tol);
                    a_status = a_classifier.state();
                }
            }
            //
            // OCCT L905-908.
            if !recadre_on_periodic || (!is_u_per && !is_v_per) {
                return a_status;
            }
            // OCCT L909-912.
            if a_status == State::In || a_status == State::On {
                return a_status;
            }
            //
            // OCCT L914-922.
            if !urecadre {
                u = uu;
                urecadre = true;
            } else if is_u_per {
                u += uperiod;
            }
            // OCCT L923-941.
            if u > self.umax || !is_u_per {
                if !vrecadre {
                    v = vv;
                    vrecadre = true;
                } else if is_v_per {
                    v += vperiod;
                }
                //
                u = uu;
                //
                if v > self.vmax || !is_v_per {
                    return a_status;
                }
            }
        } // for (;;)
    }
}

// OCCT IntTools_FClass2d::~IntTools_FClass2d — cxx L947-950:
// `TabClass.Clear();`.
impl Drop for IntToolsFClass2d {
    fn drop(&mut self) {
        self.tab_class.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Curve2d, Line2d, Line3, Plane, Surface3};
    use rcad_kernel::topods::BRepBuilder;

    /// The planar unit square face [0,1]x[0,1] at z = 0 with CCW pcurves
    /// (the fixture pattern of the brep_top_adaptor FClass2d tests).
    fn ccw_square_face() -> (std::sync::Arc<BRep>, Shape) {
        let mut brep = BRep::new();
        let mut b = BRepBuilder::new();
        let p = |x: f64, y: f64| DVec3::new(x, y, 0.0);
        let pts = [p(0.0, 0.0), p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0)];
        let pairs = [(0, 1), (1, 2), (2, 3), (3, 0)];
        let vs: Vec<Shape> = pts.iter().map(|q| b.add_vertex(&mut brep, *q, 1e-7)).collect();
        let mut edges = Vec::new();
        for &(i, j) in pairs.iter() {
            let a = pts[i];
            let c = pts[j];
            let curve = Curve3::Line(Line3::new(a, (c - a).normalize()));
            let range = [0.0, (c - a).length()];
            let v_first = vs[i].clone();
            let mut v_last = vs[j].clone();
            v_last.orientation = Orientation::Reversed;
            edges.push(b.add_edge(&mut brep, Some(curve), v_first, v_last, range));
        }
        let wire = brep.add_twire(edges.clone());
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            })),
            wire,
            Vec::new(),
            Some(DVec3::new(0.5, 0.5, 0.0)),
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );
        for (k, &(i, j)) in pairs.iter().enumerate() {
            let a2 = DVec2::new(pts[i].x, pts[i].y);
            let b2 = DVec2::new(pts[j].x, pts[j].y);
            let len = (b2 - a2).length();
            let pc = Curve2d::Line(Line2d {
                origin: a2,
                direction: (b2 - a2) / len,
            });
            b.add_pcurve(&mut brep, edges[k].clone(), face.clone(), pc, 0.0, len);
        }
        (std::sync::Arc::new(brep), face)
    }

    /// OCCT anchor (cxx L555-564): the CCW outer wire gives aS > 0 —
    /// myIsHole = false, TabOrien = 1; Perform classifies the center IN and
    /// the far point OUT through the polygon walk (cxx L685-724).
    #[test]
    fn ccw_square_center_in_far_out_outer_wire_not_hole() {
        let (brep, face) = ccw_square_face();
        let clsf = IntToolsFClass2d::new_face(brep, &face, 1e-6);
        assert!(!clsf.is_hole());
        assert_eq!(clsf.perform(DVec2::new(0.5, 0.5), true), State::In);
        assert_eq!(clsf.perform(DVec2::new(5.0, 5.0), true), State::Out);
        assert_eq!(clsf.perform_infinite_point(), State::Out);
    }

    /// OCCT anchor (cxx L639-643): an uninitialized classifier (the default
    /// ctor) has an empty TabClass — Perform returns TopAbs_IN.
    #[test]
    fn empty_classifier_performs_in() {
        let clsf = IntToolsFClass2d::new();
        assert_eq!(clsf.perform(DVec2::new(0.5, 0.5), true), State::In);
        assert_eq!(clsf.perform_infinite_point(), State::In);
    }

    /// OCCT anchor (cxx L808-816): TestOnRestriction with an empty TabClass
    /// returns TopAbs_IN.
    #[test]
    fn empty_classifier_test_on_restriction_in() {
        let clsf = IntToolsFClass2d::default();
        assert_eq!(
            clsf.test_on_restriction(DVec2::new(0.5, 0.5), 1e-6, true),
            State::In
        );
    }
}
