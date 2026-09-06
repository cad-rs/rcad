// OCCT HLRBRep_Data.cxx L1233-2683 (NextInterference / RejectedInterference /
// AboveInterference / LocalLEGeometry2D / LocalFEGeometry2D / EdgeState /
// HidingStartLevel / Compare / OrientOutLine / OrientOthEdge / Classify /
// SimplClassify / RejectedPoint / SameVertex / IsBadFace) + the file statics
// AdjustParameter (L498-518) and REJECT1 (L1950-1989) + the interference
// lxx inline MoreInterference.  Impl block of super::Data (proxy M).
//
// InitInterference (L1223-1232) belongs to [update] (proxy L):
// update.rs = L539-1232, classify.rs = L1233-2683.
//
// The OCCT element pointers (myLEData / myFEData / myLEGeom / myFEGeom /
// iFaceGeom / iFaceMinMax / myReject) are dereferenced in small unsafe faces
// mirroring the OCCT pointer statements (the HLR methodological exception).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::Line3;
use rcad_kernel::topods::{Orientation, State};

use crate::bop::int_tools::bean_face_intersector::GeomAbsCurveType;
use crate::geomalgo::int_patch::elclib;
use crate::geomalgo::int_res2d::{IntersectionPoint, Position, Situation, TypeTrans};
use crate::hlr::algo::edges_block::MinMaxIndices;
use crate::hlr::algo::hlr_algo::HLRAlgo;
use crate::hlr::algo::interference::Interference;
use crate::hlr::algo::wires_block::WiresBlock;
use crate::hlr::brep::curve::Curve;
use crate::hlr::brep::edge_data::EdgeData;
use crate::hlr::brep::face_data::FaceData;

use super::tableau_rejection::TableauRejection;
use super::{CUT_BIG, CUT_LAR, DERIVEE_PREMIERE_NULLE, counters};

/// OCCT Standard_Real RealLast() (Standard_Real.hxx).
const REAL_LAST: f64 = f64::MAX;
/// OCCT Epsilon(1.) — the double machine epsilon.
const EPSILON_1: f64 = f64::EPSILON;
/// OCCT Precision::Infinite() (Precision.hxx) = 2.e+100.
const PRECISION_INFINITE: f64 = 2.0e100;
/// OCCT Precision::IsInfinite(R) = |R| >= 0.5 * Infinite() (Precision.hxx).
fn precision_is_infinite(r: f64) -> bool {
    r.abs() >= 0.5 * PRECISION_INFINITE
}
/// OCCT Precision::PConfusion() = Confusion() * 0.01 (Precision.hxx).
const PRECISION_P_CONFUSION: f64 = 1.0e-9;
/// OCCT gp::Resolution() (gp.hxx).
const GP_RESOLUTION: f64 = 1.0e-12;
/// OCCT `& 0x80008000` — the sign-bit pair of the packed rejection test
/// (cxx L1247 etc.); the bit pattern kept as i32 (the OCCT int mask).
const REJECT_MASK: i32 = 0x80008000u32 as i32;

/// OCCT file static AdjustParameter (HLRBRep_Data.cxx L498-518).
fn adjust_parameter(e: &mut EdgeData, h: bool, p: &mut f64, t: &mut f32) {
    let mut p1: f64 = 0.0; // OCCT uninitialized locals; neutral default 0.
    let mut p2: f64 = 0.0;
    let mut t1: f32 = 0.0;
    let mut t2: f32 = 0.0;
    if h {
        // OCCT L504: E->Status().Bounds(p, t, p2, t2);
        let (b_start, b_tol_start, b_end, b_tol_end) = e.status_ref().bounds();
        *p = b_start;
        *t = b_tol_start;
        p2 = b_end;
        t2 = b_tol_end;
        if e.ver_at_sta() {
            *p = *p + (p2 - *p) * CUT_BIG;
        }
    } else {
        // OCCT L512: E->Status().Bounds(p1, t1, p, t);
        let (b_start, b_tol_start, b_end, b_tol_end) = e.status_ref().bounds();
        p1 = b_start;
        t1 = b_tol_start;
        *p = b_end;
        *t = b_tol_end;
        if e.ver_at_end() {
            *p = *p - (*p - p1) * CUT_BIG;
        }
    }
    let _ = (p1, p2, t1, t2);
}

/// OCCT file static REJECT1 (HLRBRep_Data.cxx L1950-1989, anonymous
/// namespace) — the 16-direction integer decalage box of a point box.
fn reject1(
    the_deca: &[f64; 16],
    the_tot_min: &[f64; 16],
    the_tot_max: &[f64; 16],
    the_sur_d: &[f64; 16],
    the_vert_min: &mut MinMaxIndices,
    the_vert_max: &mut MinMaxIndices,
) {
    the_vert_min.min[0] = ((the_deca[0] + the_tot_min[0]) * the_sur_d[0]) as i32;
    the_vert_max.min[0] = ((the_deca[0] + the_tot_max[0]) * the_sur_d[0]) as i32;
    the_vert_min.min[1] = ((the_deca[1] + the_tot_min[1]) * the_sur_d[1]) as i32;
    the_vert_max.min[1] = ((the_deca[1] + the_tot_max[1]) * the_sur_d[1]) as i32;
    the_vert_min.min[2] = ((the_deca[2] + the_tot_min[2]) * the_sur_d[2]) as i32;
    the_vert_max.min[2] = ((the_deca[2] + the_tot_max[2]) * the_sur_d[2]) as i32;
    the_vert_min.min[3] = ((the_deca[3] + the_tot_min[3]) * the_sur_d[3]) as i32;
    the_vert_max.min[3] = ((the_deca[3] + the_tot_max[3]) * the_sur_d[3]) as i32;
    the_vert_min.min[4] = ((the_deca[4] + the_tot_min[4]) * the_sur_d[4]) as i32;
    the_vert_max.min[4] = ((the_deca[4] + the_tot_max[4]) * the_sur_d[4]) as i32;
    the_vert_min.min[5] = ((the_deca[5] + the_tot_min[5]) * the_sur_d[5]) as i32;
    the_vert_max.min[5] = ((the_deca[5] + the_tot_max[5]) * the_sur_d[5]) as i32;
    the_vert_min.min[6] = ((the_deca[6] + the_tot_min[6]) * the_sur_d[6]) as i32;
    the_vert_max.min[6] = ((the_deca[6] + the_tot_max[6]) * the_sur_d[6]) as i32;
    the_vert_min.min[7] = ((the_deca[7] + the_tot_min[7]) * the_sur_d[7]) as i32;
    the_vert_max.min[7] = ((the_deca[7] + the_tot_max[7]) * the_sur_d[7]) as i32;
    the_vert_min.max[0] = ((the_deca[8] + the_tot_min[8]) * the_sur_d[8]) as i32;
    the_vert_max.max[0] = ((the_deca[8] + the_tot_max[8]) * the_sur_d[8]) as i32;
    the_vert_min.max[1] = ((the_deca[9] + the_tot_min[9]) * the_sur_d[9]) as i32;
    the_vert_max.max[1] = ((the_deca[9] + the_tot_max[9]) * the_sur_d[9]) as i32;
    the_vert_min.max[2] = ((the_deca[10] + the_tot_min[10]) * the_sur_d[10]) as i32;
    the_vert_max.max[2] = ((the_deca[10] + the_tot_max[10]) * the_sur_d[10]) as i32;
    the_vert_min.max[3] = ((the_deca[11] + the_tot_min[11]) * the_sur_d[11]) as i32;
    the_vert_max.max[3] = ((the_deca[11] + the_tot_max[11]) * the_sur_d[11]) as i32;
    the_vert_min.max[4] = ((the_deca[12] + the_tot_min[12]) * the_sur_d[12]) as i32;
    the_vert_max.max[4] = ((the_deca[12] + the_tot_max[12]) * the_sur_d[12]) as i32;
    the_vert_min.max[5] = ((the_deca[13] + the_tot_min[13]) * the_sur_d[13]) as i32;
    the_vert_max.max[5] = ((the_deca[13] + the_tot_max[13]) * the_sur_d[13]) as i32;
    the_vert_min.max[6] = ((the_deca[14] + the_tot_min[14]) * the_sur_d[14]) as i32;
    the_vert_max.max[6] = ((the_deca[14] + the_tot_max[14]) * the_sur_d[14]) as i32;
    the_vert_min.max[7] = ((the_deca[15] + the_tot_min[15]) * the_sur_d[15]) as i32;
    the_vert_max.max[7] = ((the_deca[15] + the_tot_max[15]) * the_sur_d[15]) as i32;
}

/// OCCT GeomInt::AdjustPeriodic (TKGeomAlgo/GeomInt/GeomInt.cxx L21-47) —
/// adjusts the parameter to the range [theParMin, theParMax]; the OCCT
/// modf(dp / thePeriod, &aNbPer) keeps the integer part.
fn geom_int_adjust_periodic(
    the_par: f64,
    the_par_min: f64,
    the_par_max: f64,
    the_period: f64,
    the_new_par: &mut f64,
    the_offset: &mut f64,
    the_eps: f64,
) -> bool {
    let b_min: bool;
    let b_max: bool;
    //
    *the_offset = 0.0;
    *the_new_par = the_par;
    b_min = the_par_min - the_par > the_eps;
    b_max = the_par - the_par_max > the_eps;
    //
    if b_min || b_max {
        let dp: f64;
        let a_nb_per: f64;
        //
        dp = if b_min {
            the_par_max - the_par
        } else {
            the_par_min - the_par
        };
        a_nb_per = (dp / the_period).trunc();
        //
        *the_offset = a_nb_per * the_period;
        *the_new_par += *the_offset;
    }
    //
    *the_offset > 0.0
}

/// OCCT TopAbs::Complement (TopAbs.hxx) — FORWARD<->REVERSED,
/// INTERNAL<->EXTERNAL.
fn top_abs_complement(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::External,
        Orientation::External => Orientation::Internal,
    }
}

/// OCCT TopAbs::Reverse (TopAbs.hxx) — FORWARD<->REVERSED only.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// OCCT HLRBRep_Curve::GetType (hxx) returns GeomAbs_CurveType; the rcad
/// Curve carries the kernel CurveType — the value-preserving type bridge
/// (the GeomAbs_CurveType / ProjLib CurveType split is an rcad artifact).
fn curve_type_to_geom_abs(t: rcad_kernel::base::proj_lib::CurveType) -> GeomAbsCurveType {
    use rcad_kernel::base::proj_lib::CurveType as CT;
    match t {
        CT::Line => GeomAbsCurveType::Line,
        CT::Circle => GeomAbsCurveType::Circle,
        CT::Ellipse => GeomAbsCurveType::Ellipse,
        CT::Hyperbola => GeomAbsCurveType::Hyperbola,
        CT::Parabola => GeomAbsCurveType::Parabola,
        CT::Bezier => GeomAbsCurveType::BezierCurve,
        CT::BSpline => GeomAbsCurveType::BSplineCurve,
        CT::Other => GeomAbsCurveType::OtherCurve,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::Data;
    use crate::hlr::brep::cl_props::CLProps;
    use crate::hlr::brep::intersector::Intersector;
    use crate::hlr::brep::sl_props::SLProps;
    use crate::topalgo::brep_top_adaptor::topol_tool_brep::BRepTopolTool;
    use glam::DVec3;
    use rcad_kernel::base::proj_lib::CurveType;
    use rcad_kernel::geom::{Plane, Point3, Surface3, Vec3};
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::topo::topods::BRepBuilder;
    use crate::geomalgo::int_patch::GeomAbsSurfaceType;
    use crate::hlr::algo::coincidence::Coincidence;
    use crate::hlr::algo::intersection::Intersection;
    use crate::hlr::algo::projector::Projector;
    use crate::hlr::algo::wires_block::WiresBlock;
    use crate::hlr::brep::b_curve_tool::CurveView;
    use crate::hlr::brep::face_iterator::FaceIterator;
    use crate::hlr::brep::surface::Surface;

    use crate::hlr::algo::interference::Interference as AlgoInterference;

    /// The shared leaked top-view projector for the test curves (the Curve
    /// keeps the OCCT `const HLRAlgo_Projector*` raw pointer).
    fn leaked_projector() -> &'static Projector {
        Box::leak(Box::new(top_view_projector()))
    }

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double).
    struct SegEdge {
        origin: Point3,
        dir: Vec3,
        len: f64,
    }

    impl CurveView for SegEdge {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.len
        }
        fn d0(&self, u: f64) -> Point3 {
            self.origin + self.dir * u
        }
        fn d1(&self, u: f64) -> (Point3, Vec3) {
            (self.origin + self.dir * u, self.dir)
        }
        fn d2(&self, u: f64) -> (Point3, Vec3, Vec3) {
            (self.origin + self.dir * u, self.dir, Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> rcad_kernel::geom::Line3 {
            rcad_kernel::geom::Line3 {
                origin: self.origin,
                direction: self.dir,
            }
        }
        fn circle(&self) -> rcad_kernel::geom::Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn degree(&self) -> i32 {
            0
        }
        fn nb_poles(&self) -> i32 {
            0
        }
        fn nb_knots(&self) -> i32 {
            0
        }
        fn is_closed(&self) -> bool {
            false
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            0.0
        }
        fn resolution(&self, r3d: f64) -> f64 {
            r3d
        }
        fn parameter_3d(&self, p2d: f64) -> f64 {
            p2d
        }
        fn poles(&self) -> Vec<Point3> {
            Vec::new()
        }
    }

    /// The top-view projector (identity transform, type 1).
    fn top_view_projector() -> Projector {
        Projector::from_ax2(&Ax2::new(
            DVec3::ZERO,
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(1.0, 0.0, 0.0),
        ))
    }

    /// The OCCT Data Set path: the loaded projected curve inside a fresh
    /// EdgeData with the full-parameter domain.
    fn line_edge_data(proj: *const Projector, edge: &'static SegEdge) -> EdgeData<'static> {
        let mut c = Curve::new();
        c.projector(proj);
        c.load(edge);
        c.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);
        let mut e = EdgeData::new();
        e.set(
            true, false, c, 0.0, 0, 0, false, false, false, false, 0.0, 0.0, edge.len, 0.0,
        );
        e
    }

    /// The plane z = 0 square face with the UV bounds [0, 2] x [0, 2]; the
    /// four wire edges come back aligned with the Data my_e_map slots.
    fn plane_brep() -> (
        &'static rcad_kernel::BRep,
        rcad_kernel::topods::Shape,
        Vec<rcad_kernel::topods::Shape>,
    ) {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = BRepBuilder::new();
        let vs = [
            b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 0.0), 1e-7),
            b.add_vertex(&mut brep, DVec3::new(2.0, 0.0, 0.0), 1e-7),
            b.add_vertex(&mut brep, DVec3::new(2.0, 2.0, 0.0), 1e-7),
            b.add_vertex(&mut brep, DVec3::new(0.0, 2.0, 0.0), 1e-7),
        ];
        let corners = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(2.0, 2.0, 0.0),
            DVec3::new(0.0, 2.0, 0.0),
        ];
        let mut edges = Vec::new();
        for k in 0..4 {
            let a = corners[k];
            let c = corners[(k + 1) % 4];
            let e = b.add_edge(
                &mut brep,
                Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: a,
                    direction: (c - a).normalize(),
                })),
                vs[k].clone(),
                vs[(k + 1) % 4].clone(),
                [0.0, 2.0],
            );
            edges.push(e);
        }
        let wire = brep.add_twire(edges.clone());
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::new(0.0, 0.0, 0.0),
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 2.0]),
            Vec::new(),
            true,
        );
        let brep: &'static rcad_kernel::BRep = Box::leak(Box::new(brep));
        (brep, face, edges)
    }

    /// The minimal HLRBRep_Data object graph for the classify anchors.  The
    /// fields mirror the OCCT ctor defaults (cxx L522-537) with the
    /// TableauRejection SetDim of the ctor.
    fn data_fixture() -> Data<'static> {
        let wires1: &'static mut WiresBlock = Box::leak(Box::new(WiresBlock::new(0)));
        let wires2: &'static mut WiresBlock = Box::leak(Box::new(WiresBlock::new(0)));
        let mut data = Data {
            my_nb_vertices: 0,
            my_nb_edges: 4,
            my_nb_faces: 1,
            my_e_map: Vec::new(),
            my_f_map: Vec::new(),
            // The kernel brep context (session-10 ruling); the uv-driven
            // tests set it together with the aligned my_e_map slots.
            my_brep: None,
            my_e_data: (0..4).map(|_| EdgeData::new()).collect(),
            my_f_data: Vec::new(),
            my_edge_indices: Vec::new(),
            my_toler: 1.0e-5,
            my_proj: top_view_projector(),
            my_l_l_props: CLProps::new(2, f64::EPSILON),
            my_f_l_props: CLProps::new(2, f64::EPSILON),
            my_s_l_props: SLProps::new(2, f64::EPSILON),
            my_big_size: 0.0,
            my_face_itr1: FaceIterator::new(wires1),
            my_face_itr2: FaceIterator::new(wires2),
            i_face: 1,
            i_face_data: std::ptr::null_mut(),
            i_face_geom: std::ptr::null_mut(),
            i_face_min_max: std::ptr::null_mut(),
            i_face_type: GeomAbsSurfaceType::Plane,
            i_face_back: false,
            i_face_simp: true,
            i_face_smpl: true,
            i_face_test: false,
            my_hide_count: 0,
            my_deca: [0.0; 16],
            my_sur_d: [1.0; 16],
            my_cur_sort_ed: 0,
            my_nbr_sort_ed: 0,
            my_le: 1,
            my_le_out_line: false,
            my_le_internal: false,
            my_le_double: false,
            my_le_iso_line: false,
            my_le_data: std::ptr::null_mut(),
            my_le_geom: std::ptr::null(),
            my_le_min_max: std::ptr::null_mut(),
            my_le_type: GeomAbsCurveType::Line,
            my_le_tol: 0.5,
            my_fe: 2,
            my_fe_ori: Orientation::Forward,
            my_fe_out_line: false,
            my_fe_internal: false,
            my_fe_double: false,
            my_fe_data: std::ptr::null_mut(),
            my_fe_geom: std::ptr::null_mut(),
            my_fe_type: GeomAbsCurveType::Line,
            my_fe_tol: 0.5,
            my_intersector: Intersector::new(),
            my_classifier: BRepTopolTool::new(std::sync::Arc::new(
                rcad_kernel::BRep::new(),
            )),
            my_same_vertex: false,
            my_intersected: false,
            my_nb_points: 0,
            my_nb_segments: 0,
            i_interf: 0,
            my_intf: Interference::new(),
            my_above_intf: false,
            my_reject: Box::into_raw(Box::new(TableauRejection::new())),
        };
        // OCCT L536: ((TableauRejection*)myReject)->SetDim(myNbEdges);
        unsafe { &mut *data.my_reject }.set_dim(data.my_nb_edges as i32);
        data
    }

    /// The interference fixture of the HidingStartLevel anchors.
    fn interference(param: f64, level: i32, transition: Orientation) -> AlgoInterference {
        let mut inters = Intersection::new();
        inters.set_parameter(param);
        inters.set_level(level);
        let mut bound = Coincidence::new();
        bound.set_2d(2, param);
        let mut i = AlgoInterference::new();
        i.set_intersection(inters);
        i.set_boundary(bound);
        i.set_transition(transition);
        i
    }

    /// The RejectedPoint fixture: LE on z = 0 (X), FE on z = 2 (Y), the
    /// TolZ of 1.0 and the FORWARD FE orientation.
    fn rejected_point_fixture() -> Data<'static> {
        let mut data = data_fixture();
        let le: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 2.0,
        }));
        let fe: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(1.0, 0.0, 2.0),
            dir: Vec3::new(0.0, 1.0, 0.0),
            len: 2.0,
        }));
        data.my_e_data[0] = line_edge_data(leaked_projector(), le);
        data.my_e_data[1] = line_edge_data(leaked_projector(), fe);
        data.my_le = 1;
        data.my_fe = 2;
        data.my_le_data = &mut data.my_e_data[0];
        data.my_le_geom = data.my_e_data[0].curve();
        data.my_fe_data = &mut data.my_e_data[1];
        data.my_fe_geom = data.my_e_data[1].curve();
        // TolZ = myBigSize * 0.00001 = 1.0 — dz = 0 - 2 = -2 <= -TolZ -> IN.
        data.my_big_size = 100000.0;
        data.my_fe_ori = Orientation::Forward;
        data.my_same_vertex = false;
        data.i_face_test = false;
        data.i_face_back = false;
        data
    }

    /// OCCT anchor: RejectedPoint (cxx L2330-2593) — the dz < -TolZ branch
    /// (state IN, the boundary under the face): the interference records the
    /// transition computed from Tr1, Or2 = INTERNAL, and the boundary
    /// Set2D(myFE, p2); no edge is rejected.
    #[test]
    fn rejected_point_under_face_state_in() {
        let mut data = rejected_point_fixture();

        let mut p_inter = IntersectionPoint::empty();
        p_inter.set_values(
            DVec2::new(1.0, 0.5),
            1.0,
            0.5,
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::In,
            ),
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::Out,
            ),
            false,
        );

        let rejected = data.rejected_point(&p_inter, Orientation::Forward, 0);
        assert!(!rejected);
        assert!(!data.my_above_intf);

        let inter = data.my_intf.intersection();
        // Ori comes from Tr1 PositionOnCurve = Middle -> INTERNAL; with an
        // INTERNAL Ori the decal stays 1 (st == IN but Ori != FORWARD).
        assert_eq!(inter.orientation(), Orientation::Internal);
        assert_eq!(inter.level(), 1, "decal: IN + INTERNAL Ori");
        assert_eq!(inter.seg_index(), 0);
        assert_eq!(inter.index(), 0);
        assert!((inter.parameter() - 1.0).abs() < 1e-12);
        assert_eq!(inter.tolerance(), 0.5);
        assert_eq!(inter.state(), State::In);
        assert_eq!(data.my_intf.orientation(), Orientation::Internal);
        assert_eq!(data.my_intf.transition(), Orientation::Forward);
        assert_eq!(data.my_intf.boundary_transition(), Orientation::Forward);
        let (fe_2d, param_2d) = data.my_intf.boundary().value_2d();
        assert_eq!(fe_2d, 2);
        assert!((param_2d - 0.5).abs() < 1e-12);

        // The Tr1 Head variant: Ori = FORWARD, Orie = FORWARD, st = IN ->
        // decal 0.
        let mut data2 = rejected_point_fixture();
        let mut p_inter2 = IntersectionPoint::empty();
        p_inter2.set_values(
            DVec2::new(1.0, 0.5),
            1.0,
            0.5,
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Head,
                TypeTrans::In,
            ),
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::Out,
            ),
            false,
        );
        assert!(!data2.rejected_point(&p_inter2, Orientation::Forward, 3));
        assert_eq!(data2.my_intf.intersection().orientation(), Orientation::Forward);
        assert_eq!(data2.my_intf.intersection().level(), 0, "decal 0");
        assert_eq!(data2.my_intf.intersection().seg_index(), 3);
    }

    /// OCCT anchor: RejectedPoint (cxx L2362-2366) — the dz >= TolZ branch
    /// rejects the point and raises AboveInterference.
    #[test]
    fn rejected_point_above_face_rejected() {
        let mut data = data_fixture();
        let le: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(0.0, 0.0, 2.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 2.0,
        }));
        let fe: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(1.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 1.0, 0.0),
            len: 2.0,
        }));
        data.my_e_data[0] = line_edge_data(leaked_projector(), le);
        data.my_e_data[1] = line_edge_data(leaked_projector(), fe);
        data.my_le = 1;
        data.my_fe = 2;
        data.my_le_data = &mut data.my_e_data[0];
        data.my_le_geom = data.my_e_data[0].curve();
        data.my_fe_data = &mut data.my_e_data[1];
        data.my_fe_geom = data.my_e_data[1].curve();
        data.my_big_size = 100000.0;

        let mut p_inter = IntersectionPoint::empty();
        p_inter.set_values(
            DVec2::ZERO,
            1.0,
            1.0,
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::In,
            ),
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::Out,
            ),
            false,
        );

        let rejected = data.rejected_point(&p_inter, Orientation::Forward, 0);
        assert!(rejected);
        assert!(data.my_above_intf);
    }

    /// OCCT anchor: RejectedPoint reject paths (cxx L2394-2399, L2424-2429)
    /// — the same-vertex end intersection and the Undecided transition.
    #[test]
    fn rejected_point_reject_paths() {
        let mut data = data_fixture();
        let le: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 2.0,
        }));
        let fe: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(1.0, 0.0, 2.0),
            dir: Vec3::new(0.0, 1.0, 0.0),
            len: 2.0,
        }));
        data.my_e_data[0] = line_edge_data(leaked_projector(), le);
        data.my_e_data[1] = line_edge_data(leaked_projector(), fe);
        data.my_le = 1;
        data.my_fe = 2;
        data.my_le_data = &mut data.my_e_data[0];
        data.my_le_geom = data.my_e_data[0].curve();
        data.my_fe_data = &mut data.my_e_data[1];
        data.my_fe_geom = data.my_e_data[1].curve();
        data.my_big_size = 100000.0;

        // mySameVertex + Tr2 not Middle -> rejected (L2455).
        data.my_same_vertex = true;
        let mut p_inter = IntersectionPoint::empty();
        p_inter.set_values(
            DVec2::ZERO,
            1.0,
            1.0,
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::In,
            ),
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Head,
                TypeTrans::In,
            ),
            false,
        );
        assert!(data.rejected_point(&p_inter, Orientation::Forward, 0));

        // Tr1 UNDECIDED -> rejected (L2428-2429); the flag stays untouched.
        data.my_same_vertex = false;
        let mut p_inter2 = IntersectionPoint::empty();
        p_inter2.set_values(
            DVec2::ZERO,
            1.0,
            1.0,
            crate::geomalgo::int_res2d::Transition::undecided(Position::Middle),
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::Out,
            ),
            false,
        );
        assert!(data.rejected_point(&p_inter2, Orientation::Forward, 0));

        // Tr1 TOUCH with an UNKNOWN situation -> rejected (L2424-2425).
        let mut p_inter3 = IntersectionPoint::empty();
        p_inter3.set_values(
            DVec2::ZERO,
            1.0,
            1.0,
            crate::geomalgo::int_res2d::Transition::touch(
                false,
                Position::Middle,
                Situation::Unknown,
                false,
            ),
            crate::geomalgo::int_res2d::Transition::in_out(
                false,
                Position::Middle,
                TypeTrans::Out,
            ),
            false,
        );
        assert!(data.rejected_point(&p_inter3, Orientation::Forward, 0));
    }

    /// OCCT anchor: SameVertex (cxx L2597-2651) — the VSta/VEnd equality and
    /// the OutLV/CutAt flag effects on myIntersected.
    #[test]
    fn same_vertex_flag_semantics() {
        let mut data = data_fixture();
        let le: &'static EdgeData = {
            let mut e = EdgeData::new();
            e.set_v_sta(7);
            e.set_v_end(8);
            e.set_out_lv_sta(true);
            e.set_cut_at_sta(false);
            Box::leak(Box::new(e))
        };
        let fe: &'static EdgeData = {
            let mut e = EdgeData::new();
            e.set_v_sta(7);
            e.set_v_end(9);
            Box::leak(Box::new(e))
        };
        data.my_le_data = le as *const EdgeData as *mut EdgeData;
        data.my_fe_data = fe as *const EdgeData as *mut EdgeData;
        data.my_le_type = GeomAbsCurveType::Line;
        data.my_fe_type = GeomAbsCurveType::Line;
        data.i_face_test = false;
        data.my_le_internal = false;

        // h1/h2 at the starts: v1 == v2 -> same vertex; Line/Line -> no
        // other intersection.
        assert!(data.same_vertex(true, true));
        assert!(!data.my_intersected);

        // OutLVSta without iFaceTest/internal keeps otherCase, no CutAtSta
        // -> the intersections stay wanted.
        data.my_le_type = GeomAbsCurveType::BSplineCurve;
        assert!(data.same_vertex(true, true));
        assert!(data.my_intersected);

        // CutAtSta on the LE start kills the intersection (two connected
        // OutLines do not intersect themselves).
        unsafe { &mut *(data.my_le_data) }.set_cut_at_sta(true);
        assert!(data.same_vertex(true, true));
        assert!(!data.my_intersected);

        // iFaceTest with OutLVSta takes the else-if branch: otherCase =
        // false, the CutAtSta check is skipped.
        unsafe { &mut *(data.my_le_data) }.set_cut_at_sta(false);
        data.i_face_test = true;
        assert!(data.same_vertex(true, true));
        assert!(data.my_intersected);
        data.i_face_test = false;

        // Different vertices (LE end vs FE start) -> not the same vertex.
        assert!(!data.same_vertex(false, true));
    }

    /// OCCT anchor: HidingStartLevel (cxx L1658-1733) — the two-pass level
    /// accumulation over a small IL with the far hiding-face box rejecting
    /// the inner Classify (level contribution 0), then
    /// Reversed += 2 / Forward -= 3 below param - tolpar.
    #[test]
    fn hiding_start_level_small_il() {
        let mut data = data_fixture();
        let proj: &'static Projector = Box::leak(Box::new(top_view_projector()));
        let le: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 2.0,
        }));
        let ed: &'static EdgeData = Box::leak(Box::new(line_edge_data(proj, le)));
        // A far-away hiding-face box: every encoded comparison rejects.
        data.i_face_min_max = Box::into_raw(Box::new(MinMaxIndices {
            min: [100000; 8],
            max: [100001; 8],
        }));

        // sta = 0, end = 2; 0.5 and 0.6 snap sta (closer to sta each time):
        // param = 0.5 * (0.6 + 2) = 1.3, tolpar = 0.02.
        let il = vec![
            interference(0.5, 2, Orientation::Reversed),
            interference(0.6, 3, Orientation::Forward),
        ];
        let level = data.hiding_start_level(1, ed, &il);
        assert_eq!(level, -1, "Reversed += 2, Forward -= 3");

        // An interference above the end stops the first pass (param from
        // sta/end only) and the second pass (level stays 0).
        let il2 = vec![interference(5.0, 4, Orientation::Reversed)];
        let level2 = data.hiding_start_level(1, ed, &il2);
        assert_eq!(level2, 0);
    }

    /// OCCT anchor: EdgeState (cxx L1580-1654) — on the z = 0 hiding plane:
    /// the coplanar tangent gives scal 0 (ON/ON); the transversal tangent
    /// along +Z against the reversed plane normal gives scal < 0 (OUT/IN).
    #[test]
    fn edge_state_states_around_plane() {
        let mut data = data_fixture();
        let (brep, face, edges) = plane_brep();
        let mut hsurf = Surface::new();
        hsurf.load(brep, &face);
        let hsurf: &'static Surface<'static> = Box::leak(Box::new(hsurf));
        data.i_face_geom = hsurf as *const Surface<'static> as *mut Surface<'static>;
        // OCCT L1037 (InitEdge): mySLProps.SetSurface(iFaceGeom) before the
        // exploration; the test mirrors it on the leaked surface.
        data.my_s_l_props.set_surface(hsurf);
        // The kernel context: the brep and the my_e_map slots aligned with
        // my_e_data (session-10 ruling).  myFE = 2 reads slot 1.
        data.my_brep = Some(brep);
        data.my_e_map = edges;

        // The coplanar edge along X on the plane (tangent (1, 0, 0), the
        // reversed plane normal (0, 0, -1): scal = 0 -> ON/ON).
        let le: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 2.0,
        }));
        let mut ed = line_edge_data(leaked_projector(), le);
        let mut c = Curve::new();
        c.projector(leaked_projector());
        c.load(le);
        c.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);
        data.my_fe_geom = &mut c as *mut Curve<'static>;
        data.my_le_geom = &mut c as *const Curve<'static>;
        let _ = &mut ed;

        let mut stbef = State::In;
        let mut staft = State::In;
        data.edge_state(1.0, 1.0, &mut stbef, &mut staft);
        assert_eq!(stbef, State::On, "coplanar: scal ~ 0");
        assert_eq!(staft, State::On, "coplanar: scal ~ 0");

        // The transversal edge along +Z at (1, 1): tangent (0, 0, 1) against
        // the reversed normal (0, 0, -1): scal = -1 -> OUT/IN.
        let vert: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(1.0, 1.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
            len: 2.0,
        }));
        let mut c2 = Curve::new();
        c2.projector(leaked_projector());
        c2.load(vert);
        c2.update(&mut [0.0f64; 16], &mut [0.0f64; 16]);
        data.my_le_geom = &mut c2 as *const Curve<'static>;

        let mut stbef2 = State::In;
        let mut staft2 = State::In;
        data.edge_state(1.0, 1.0, &mut stbef2, &mut staft2);
        assert_eq!(stbef2, State::Out, "transversal: scal = -1");
        assert_eq!(staft2, State::In, "transversal: scal = -1");
    }

    /// The plane hiding-face FaceData with one outline edge (edge index 1)
    /// oriented as requested, plus the Data with the matching iFaceGeom.
    fn orient_out_line_fixture(
        block_ori: Orientation,
        closed: bool,
    ) -> (Data<'static>, FaceData<'static>) {
        let (brep, face, edges) = plane_brep();
        let mut fd = FaceData::new();
        fd.set(brep, &face, Orientation::Forward, closed, 1);
        fd.set_wire(1, 1);
        fd.set_w_edge(1, 1, 1, block_ori, true, false, false, false);

        let mut data = data_fixture();
        let mut hsurf = Surface::new();
        hsurf.load(brep, &face);
        let hsurf: &'static Surface<'static> = Box::leak(Box::new(hsurf));
        data.i_face_geom = hsurf as *const Surface<'static> as *mut Surface<'static>;
        // OCCT L1037 (InitEdge): mySLProps.SetSurface(iFaceGeom).
        data.my_s_l_props.set_surface(hsurf);
        data.my_brep = Some(brep);
        data.my_e_map = edges;
        let proj = leaked_projector();
        let le: &'static SegEdge = Box::leak(Box::new(SegEdge {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(1.0, 0.0, 0.0),
            len: 2.0,
        }));
        data.my_e_data[0] = line_edge_data(proj, le);
        (data, fd)
    }

    /// OCCT anchor: OrientOutLine (cxx L1746-1877) — the outline edge on the
    /// z = 0 plane gets the REVERSED orientation (Nm x Tg = (0,-1,0), r = 0
    /// -> not > 0), written back into the wire block; with a matching block
    /// orientation and an open face there is no inversion; with a mismatched
    /// block orientation on a closed face the face flips (inverted = true).
    #[test]
    fn orient_out_line_orientation_and_inversion() {
        // Sub-case A: REVERSED block on an open face — no inversion.
        let (mut data, mut fd) = orient_out_line_fixture(Orientation::Reversed, false);
        let inverted = data.orient_out_line(1, &mut fd);
        assert!(!inverted);
        assert_eq!(data.my_fe_ori, Orientation::Reversed, "r = 0 -> REVERSED");
        let wb_a: *mut WiresBlock = match fd.change_wires().as_mut() {
            Some(w) => std::sync::Arc::as_ptr(w) as *const WiresBlock as *mut WiresBlock,
            None => std::ptr::null_mut(),
        };
        assert_eq!(unsafe { &mut *wb_a }.wire(1).orientation(1), Orientation::Reversed);

        // Sub-case B: FORWARD block on a closed face — the face flips.
        let (mut data2, mut fd2) = orient_out_line_fixture(Orientation::Forward, true);
        let inverted2 = data2.orient_out_line(1, &mut fd2);
        assert!(inverted2);
        assert_eq!(data2.my_fe_ori, Orientation::Reversed);
        let wb_b: *mut WiresBlock = match fd2.change_wires().as_mut() {
            Some(w) => std::sync::Arc::as_ptr(w) as *const WiresBlock as *mut WiresBlock,
            None => std::ptr::null_mut(),
        };
        assert_eq!(unsafe { &mut *wb_b }.wire(1).orientation(1), Orientation::Reversed);
    }
}

impl<'a> super::Data<'a> {
    // TEMPORARY diagnostic (to be removed): RCAD_HLR_TRACE=1 enables a
    // run trace of orient_out_line and next_interference.
    pub(crate) fn trace_enabled() -> bool {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        *ON.get_or_init(|| std::env::var("RCAD_HLR_TRACE").is_ok())
    }

    pub(crate) fn trace_i_face(&self) -> usize {
        self.i_face
    }

    pub(crate) fn trace_my_fe(&self) -> usize {
        self.my_fe
    }

    pub(crate) fn trace_edge_ends(&self, e: usize) -> (DVec3, DVec3) {
        let ed = &self.my_e_data[e - 1];
        let ec = ed.geometry();
        let sta = ec.parameter_3d(ec.first_parameter());
        let end = ec.parameter_3d(ec.last_parameter());
        (ec.value_3d(sta), ec.value_3d(end))
    }

    /// OCCT MoreInterference (HLRBRep_Data.lxx L103-107) — the interference
    /// exploration query (part of the interference group of [classify]).
    pub fn more_interference(&self) -> bool {
        self.i_interf <= self.my_nb_points + 2 * self.my_nb_segments
    }

    /// OCCT NextInterference (cxx L1233-1482).
    pub fn next_interference(&mut self) {
        // are there more intersections on the current edge
        self.i_interf += 1;
        //  int miniWire1,miniWire2;
        //  int maxiWire1,maxiWire2,maxiWire3,maxiWire4;

        while !self.more_interference() && self.my_face_itr1.more_edge() {
            // rejection of current wire
            if self.my_face_itr1.beginning_of_wire() {
                // OCCT L1246: MinMaxWire = myFaceItr1.Wire()->MinMax()
                let min_max_wire = *self.my_face_itr1.wire().min_max();
                let my_le_min_max = unsafe { &*self.my_le_min_max };
                if ((min_max_wire.max[0].wrapping_sub(my_le_min_max.min[0])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[0].wrapping_sub(min_max_wire.min[0])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[1].wrapping_sub(my_le_min_max.min[1])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[1].wrapping_sub(min_max_wire.min[1])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[2].wrapping_sub(my_le_min_max.min[2])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[2].wrapping_sub(min_max_wire.min[2])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[3].wrapping_sub(my_le_min_max.min[3])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[3].wrapping_sub(min_max_wire.min[3])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[4].wrapping_sub(my_le_min_max.min[4])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[4].wrapping_sub(min_max_wire.min[4])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[5].wrapping_sub(my_le_min_max.min[5])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[5].wrapping_sub(min_max_wire.min[5])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[6].wrapping_sub(my_le_min_max.min[6])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[6].wrapping_sub(min_max_wire.min[6])) & REJECT_MASK) != 0
                    || ((min_max_wire.max[7].wrapping_sub(my_le_min_max.min[7])) & REJECT_MASK) != 0
                    || ((my_le_min_max.max[7].wrapping_sub(min_max_wire.min[7])) & REJECT_MASK) != 0
                {
                    //-- Rejection en Z
                    self.my_face_itr1.skip_wire();
                    continue;
                }
            }
            self.my_fe = self.my_face_itr1.edge() as usize;
            self.my_fe_ori = self.my_face_itr1.orientation();
            self.my_fe_out_line = self.my_face_itr1.out_line();
            self.my_fe_internal = self.my_face_itr1.internal();
            self.my_fe_double = self.my_face_itr1.double();
            if Self::trace_enabled() {
                eprintln!(
                    "[TRACE] next_interf: iFace={} LE={} FE={} ori={:?} out={} int={} dbl={} test={} cur={} nbr={}",
                    self.i_face,
                    self.my_le,
                    self.my_fe,
                    self.my_fe_ori,
                    self.my_fe_out_line,
                    self.my_fe_internal,
                    self.my_fe_double,
                    self.i_face_test,
                    self.my_cur_sort_ed,
                    self.my_nbr_sort_ed,
                );
            }
            // OCCT L1272-1275: myFEData / myFEGeom / myFETol / myFEType.
            self.my_fe_data = &mut self.my_e_data[(self.my_fe - 1) as usize];
            self.my_fe_geom =
                unsafe { &mut *self.my_fe_data }.change_geometry() as *mut Curve<'a>;
            self.my_fe_tol = unsafe { &*self.my_fe_data }.tolerance();
            self.my_fe_type = curve_type_to_geom_abs(unsafe { &*self.my_fe_geom }.get_type());

            if self.my_fe_ori == Orientation::Forward || self.my_fe_ori == Orientation::Reversed {
                // Edge from the boundary
                // OCCT L1280: not a vertical edge and not a double Edge.
                if !unsafe { &*self.my_fe_data }.vertical()
                    && (!self.my_fe_double || self.my_fe_out_line)
                {
                    // OCCT L1283: MinMaxFEdg = &((HLRBRep_EdgeData*)myFEData)->MinMax();
                    let min_max_fedg: &MinMaxIndices =
                        unsafe { &mut *self.my_fe_data }.min_max();
                    let my_le_min_max = unsafe { &*self.my_le_min_max };
                    //-- -----------------------------------------------------------------------
                    //-- Max - Min doit etre positif pour toutes les directions
                    //--
                    //-- Rejection 1   (FEMax-LEMin)& REJECT_MASK  !=0
                    //--
                    //--                   FE Min ...........  FE Max
                    //--                                                LE Min ....   LE Max
                    //--
                    //-- Rejection 2   (LEMax-FEMin)& REJECT_MASK  !=0
                    //--                            FE Min ...........  FE Max
                    //--     LE Min ....   LE Max
                    //-- ----------------------------------------------------------------------
                    if !unsafe { &mut *self.my_reject }.no_intersection(
                        self.my_le as i32,
                        self.my_fe as i32,
                    ) {
                        if ((min_max_fedg.max[0].wrapping_sub(my_le_min_max.min[0]))
                            & REJECT_MASK)
                            == 0
                            && ((my_le_min_max.max[0].wrapping_sub(min_max_fedg.min[0]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[1].wrapping_sub(my_le_min_max.min[1]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[1].wrapping_sub(min_max_fedg.min[1]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[2].wrapping_sub(my_le_min_max.min[2]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[2].wrapping_sub(min_max_fedg.min[2]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[3].wrapping_sub(my_le_min_max.min[3]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[3].wrapping_sub(min_max_fedg.min[3]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[4].wrapping_sub(my_le_min_max.min[4]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[4].wrapping_sub(min_max_fedg.min[4]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[5].wrapping_sub(my_le_min_max.min[5]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[5].wrapping_sub(min_max_fedg.min[5]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[6].wrapping_sub(my_le_min_max.min[6]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[6].wrapping_sub(min_max_fedg.min[6]))
                                & REJECT_MASK)
                                == 0
                            && ((min_max_fedg.max[7].wrapping_sub(my_le_min_max.min[7]))
                                & REJECT_MASK)
                                == 0
                            && ((my_le_min_max.max[7].wrapping_sub(min_max_fedg.min[7]))
                                & REJECT_MASK)
                                == 0
                        {
                            //-- Rejection en Z
                            // not rejected perform intersection
                            let mut rej = false;
                            if self.my_le == self.my_fe {
                                // test if an auto-intersection is not useful
                                // OCCT L1320-1327: the AutoIntersectionDone
                                // read re-sets the flag before Simple().
                                if unsafe { &*self.my_le_data }.auto_intersection_done() {
                                    unsafe { &mut *self.my_le_data }
                                        .set_auto_intersection_done(true);
                                    if unsafe { &*self.my_le_data }.simple() {
                                        rej = true;
                                    }
                                }
                            }
                            if !rej {
                                counters::NB_CAL1_INTERSECTION
                                    .with(|c| c.set(c.get() + 1));
                                let mut h1 = false;
                                let mut e1 = false;
                                let mut h2 = false;
                                let mut e2 = false;
                                self.my_same_vertex = false;

                                if self.my_le == self.my_fe {
                                    self.my_intersected = true;
                                    self.my_same_vertex = false;
                                } else {
                                    self.my_intersected = true;
                                    if self.same_vertex(true, true) {
                                        self.my_same_vertex = true;
                                        h1 = true;
                                        h2 = true;
                                    }
                                    if self.same_vertex(true, false) {
                                        self.my_same_vertex = true;
                                        h1 = true;
                                        e2 = true;
                                    }
                                    if self.same_vertex(false, true) {
                                        self.my_same_vertex = true;
                                        e1 = true;
                                        h2 = true;
                                    }
                                    if self.same_vertex(false, false) {
                                        self.my_same_vertex = true;
                                        e1 = true;
                                        e2 = true;
                                    }
                                }

                                self.my_nb_points = 0;
                                self.my_nb_segments = 0;
                                self.i_interf = 1;

                                if self.my_intersected {
                                    // compute real intersection
                                    counters::NB_CAL2_INTERSECTION
                                        .with(|c| c.set(c.get() + 1));

                                    let mut da1: f64 = 0.0;
                                    let mut db1: f64 = 0.0;
                                    let mut da2: f64 = 0.0;
                                    let mut db2: f64 = 0.0;

                                    if self.my_same_vertex || self.my_le == self.my_fe {
                                        if h1 {
                                            da1 = CUT_LAR;
                                        }
                                        if e1 {
                                            db1 = CUT_LAR;
                                        }
                                        if h2 {
                                            da2 = CUT_LAR;
                                        }
                                        if e2 {
                                            db2 = CUT_LAR;
                                        }
                                    }
                                    let mut no_inter: i32 = 0;
                                    if self.my_le == self.my_fe {
                                        self.my_intersector.perform_auto(
                                            unsafe { &mut *self.my_le_data },
                                            da1,
                                            db1,
                                        );
                                    } else {
                                        let mut su: f64 = 0.0;
                                        let mut sv: f64 = 0.0;
                                        unsafe { &mut *self.my_reject }.get_single_intersection(
                                            self.my_le as i32,
                                            self.my_fe as i32,
                                            &mut su,
                                            &mut sv,
                                        );
                                        if su != REAL_LAST {
                                            self.my_intersector.simulate_one_point(
                                                unsafe { &*self.my_le_data },
                                                su,
                                                unsafe { &*self.my_fe_data },
                                                sv,
                                            );
                                            //-- std::cout<<"p";
                                        } else {
                                            self.my_intersector.perform(
                                                self.my_le as i32,
                                                unsafe { &*self.my_le_data },
                                                da1,
                                                db1,
                                                self.my_fe as i32,
                                                unsafe { &*self.my_fe_data },
                                                da2,
                                                db2,
                                                self.my_same_vertex,
                                            );
                                            if self.my_intersector.is_done() {
                                                if self.my_intersector.nb_points() == 1
                                                    && self.my_intersector.nb_segments() == 0
                                                {
                                                    let p1 =
                                                        self.my_intersector.point(1).clone();
                                                    unsafe { &mut *self.my_reject }
                                                        .set_intersection(
                                                            self.my_le as i32,
                                                            self.my_fe as i32,
                                                            &p1,
                                                        );
                                                }
                                            }
                                        }
                                        no_inter = 0;
                                    }
                                    if no_inter != 0 {
                                        self.my_nb_points = 0;
                                        self.my_nb_segments = 0;
                                    } else {
                                        if self.my_intersector.is_done() {
                                            self.my_nb_points = self.my_intersector.nb_points();
                                            self.my_nb_segments =
                                                self.my_intersector.nb_segments();
                                            if (self.my_nb_segments + self.my_nb_points) > 0 {
                                                counters::NB_OK_INTERSECTION
                                                    .with(|c| c.set(c.get() + 1));
                                            } else {
                                                unsafe { &mut *self.my_reject }
                                                    .set_no_intersection(
                                                        self.my_le as i32,
                                                        self.my_fe as i32,
                                                    );
                                            }
                                        } else {
                                            self.my_nb_points = 0;
                                            self.my_nb_segments = 0;
                                            // OCCT L1454-1461: the
                                            // #ifdef OCCT_DEBUG "Intersection
                                            // not done" print is not translated.
                                        }
                                    }
                                }
                                counters::NB_PT_INTERSECTION
                                    .with(|c| c.set(c.get() + self.my_nb_points as i32));
                                counters::NB_SEG_INTERSECTION
                                    .with(|c| c.set(c.get() + self.my_nb_segments as i32));
                                if Self::trace_enabled() {
                                    eprintln!(
                                        "[TRACE]   -> LE={}({:?}) FE={}({:?}) done={} pts={} segs={}",
                                        self.my_le,
                                        unsafe { &*self.my_le_geom }.get_type(),
                                        self.my_fe,
                                        unsafe { &*self.my_fe_geom }.get_type(),
                                        self.my_intersector.is_done(),
                                        self.my_nb_points,
                                        self.my_nb_segments,
                                    );
                                    if self.my_nb_points == 0 && self.my_nb_segments == 0 {
                                        let le = unsafe { &*self.my_le_geom };
                                        let fe = unsafe { &*self.my_fe_geom };
                                        let samp = |c: &crate::hlr::brep::curve::Curve| {
                                            let (u1, u2) = (c.first_parameter(), c.last_parameter());
                                            let mut s = String::new();
                                            for k in 0..=4 {
                                                let u = u1 + (u2 - u1) * (k as f64) / 4.0;
                                                let p = c.value(u);
                                                s.push_str(&format!("({:.4},{:.4})", p.x, p.y));
                                            }
                                            format!("[{:.4},{:.4}] {}", u1, u2, s)
                                        };
                                        eprintln!(
                                            "[TRACE]      LE{} {:?} {}",
                                            self.my_le,
                                            le.get_type(),
                                            samp(le)
                                        );
                                        eprintln!(
                                            "[TRACE]      FE{} {:?} {}",
                                            self.my_fe,
                                            fe.get_type(),
                                            samp(fe)
                                        );
                                    }
                                }
                            }
                        } else {
                        }
                    } else {
                        //-- std::cout<<"+";
                    }
                }
            }
            // next edge in face
            self.my_face_itr1.next_edge();
        }
    }

    /// OCCT RejectedInterference (cxx L1486-1521).
    pub fn rejected_interference(&mut self) -> bool {
        if self.i_interf <= self.my_nb_points {
            // OCCT L1490: RejectedPoint(myIntersector.Point(iInterf), ...);
            // the rcad IntersectionPoint is cloned (OCCT passes a const ref;
            // the point is only read).
            let p = self.my_intersector.point(self.i_interf).clone();
            return self.rejected_point(&p, Orientation::External, 0);
        } else {
            let n = (self.i_interf - self.my_nb_points) as i32;
            let mut first_point = (n & 1) != 0;
            let mut nseg = n >> 1;
            if first_point {
                nseg += 1;
            }
            let my_le_geom = unsafe { &*self.my_le_geom };
            let pf = my_le_geom.parameter_3d(
                self.my_intersector
                    .segment(nseg as usize)
                    .first_point()
                    .param_on_first(),
            );
            let pl = my_le_geom.parameter_3d(
                self.my_intersector
                    .segment(nseg as usize)
                    .last_point()
                    .param_on_first(),
            );
            if pf > pl {
                first_point = !first_point;
            }

            if first_point {
                let pt = self
                    .my_intersector
                    .segment(nseg as usize)
                    .first_point()
                    .clone();
                let ret1 = self.rejected_point(&pt, Orientation::Forward, nseg);
                return ret1;
            } else {
                let pt = self
                    .my_intersector
                    .segment(nseg as usize)
                    .last_point()
                    .clone();
                let ret2 = self.rejected_point(&pt, Orientation::Reversed, -nseg);
                return ret2;
            }
        }
    }

    /// OCCT AboveInterference (cxx L1525-1528).
    pub fn above_interference(&self) -> bool {
        self.my_above_intf
    }

    /// OCCT LocalLEGeometry2D (cxx L1532-1549).
    pub fn local_le_geometry_2d(&mut self, param: f64, tg: &mut DVec2, nm: &mut DVec2, cu: &mut f64) {
        self.my_l_l_props.set_parameter(param);
        if !self.my_l_l_props.is_tangent_defined() {
            // OCCT L1537: throw Standard_Failure("HLRBRep_Data::LocalGeometry2D");
            panic!("HLRBRep_Data::LocalGeometry2D");
        }
        *tg = self
            .my_l_l_props
            .tangent()
            .expect("HLRBRep_Data::LocalGeometry2D");
        *cu = self.my_l_l_props.curvature();
        if *cu > EPSILON_1 && !precision_is_infinite(*cu) {
            *nm = self
                .my_l_l_props
                .normal()
                .expect("HLRBRep_Data::LocalGeometry2D");
        } else {
            *nm = DVec2::new(-tg.y, tg.x);
        }
    }

    /// OCCT LocalFEGeometry2D (cxx L1553-1576).
    pub fn local_fe_geometry_2d(
        &mut self,
        fe: i32,
        param: f64,
        tg: &mut DVec2,
        nm: &mut DVec2,
        cu: &mut f64,
    ) {
        // OCCT L1559: const HLRBRep_Curve* aCurve = &myEData(FE).ChangeGeometry();
        // (the OCCT pointer into the myEData array; the rcad raw-pointer
        // deref reproduces it — the ClProps curve reference keeps the 'a
        // view lifetime like the OCCT HLRBRep_CurvePtr).
        let a_curve_ptr: *const Curve<'a> =
            &*self.my_e_data[(fe - 1) as usize].geometry();
        let a_curve: &'a Curve<'a> = unsafe { &*a_curve_ptr };
        self.my_f_l_props.set_curve(a_curve);
        self.my_f_l_props.set_parameter(param);
        if !self.my_f_l_props.is_tangent_defined() {
            // OCCT L1564: throw Standard_Failure("HLRBRep_Data::LocalGeometry2D");
            panic!("HLRBRep_Data::LocalGeometry2D");
        }
        *tg = self
            .my_f_l_props
            .tangent()
            .expect("HLRBRep_Data::LocalGeometry2D");
        *cu = self.my_f_l_props.curvature();
        if *cu > EPSILON_1 && !precision_is_infinite(*cu) {
            *nm = self
                .my_f_l_props
                .normal()
                .expect("HLRBRep_Data::LocalGeometry2D");
        } else {
            *nm = DVec2::new(-tg.y, tg.x);
        }
    }

    /// OCCT EdgeState (cxx L1580-1654).
    pub fn edge_state(&mut self, p1: f64, p2: f64, stbef: &mut State, staft: &mut State) {
        // compute the state of The Edge near the Intersection
        // this method should give the states before and after
        // it should get the parameters on the surface

        let mut pu: f64 = 0.0; // OCCT uninitialized out params; neutral default 0.
        let mut pv: f64 = 0.0;
        // OCCT L1590: HLRBRep_EdgeFaceTool::UVPoint(p2, myFEGeom, iFaceGeom, pu, pv)
        // — the OCCT tool reads myFEGeom->Curve().Edge() internally; the rcad
        // call passes the Data brep and the aligned my_e_map slot (the
        // session-10 contract ruling).
        if crate::hlr::brep::edge_face_tool::uv_point(
            self.my_brep.expect("kernel brep context"),
            p2,
            self.my_fe_geom,
            &self.my_e_map[(self.my_fe - 1) as usize],
            self.i_face_geom,
            &mut pu,
            &mut pv,
        ) {
            self.my_s_l_props.set_parameters(pu, pv);
            if self.my_s_l_props.is_normal_defined() {
                // OCCT L1595: gp_Dir NrmFace = mySLProps.Normal();
                let mut nrm_face = self
                    .my_s_l_props
                    .normal()
                    .expect("HLRBRep_Data::EdgeState");

                // OCCT L1597-1599: gp_Pnt Pbid; gp_Vec TngEdge; D1 (3D, lxx L59-63).
                let (pbid, tng_edge) = unsafe { &*self.my_le_geom }.view().d1(p1);
                let _ = pbid;

                // OCCT L1601: const gp_Trsf& TI = myProj.InvertedTransformation();
                let ti = *self.my_proj.inverted_transformation();
                let v = if self.my_proj.perspective() {
                    let mut p2d = DVec2::ZERO;
                    self.my_proj.project_pnt(pbid, &mut p2d);
                    DVec3::new(p2d.x, p2d.y, -self.my_proj.focus())
                } else {
                    DVec3::new(0.0, 0.0, 1.0) // gp_Dir(gp_Dir::D::NZ)
                };
                let v = ti.transform_dir(v); // V.Transform(TI);
                if nrm_face.dot(v) > 0.0 {
                    nrm_face = -nrm_face; // NrmFace.Reverse();
                }

                // OCCT L1619: gp_Dir(TngEdge) normalizes the tangent.
                let scal = if tng_edge.length_squared() > 1.0e-10 {
                    nrm_face.dot(tng_edge.normalize())
                } else {
                    0.0
                };

                if scal > self.my_toler as f64 * 10.0 {
                    *stbef = State::In;
                    *staft = State::Out;
                } else if scal < -(self.my_toler as f64) * 10.0 {
                    *stbef = State::Out;
                    *staft = State::In;
                } else {
                    *stbef = State::On;
                    *staft = State::On;
                }
            } else {
                *stbef = State::Out;
                *staft = State::Out;
                // OCCT L1641-1643: the #ifdef OCCT_DEBUG "undefined" print
                // is not translated.
            }
        } else {
            *stbef = State::Out;
            *staft = State::Out;
            // OCCT L1650-1652: the #ifdef OCCT_DEBUG "undefined" print
            // is not translated.
        }
    }

    /// OCCT HidingStartLevel (cxx L1658-1733) — IL is the OCCT
    /// NCollection_List<HLRAlgo_Interference> as a slice (the List → Vec
    /// mapping).
    pub fn hiding_start_level(&mut self, e: i32, ed: &EdgeData<'a>, il: &[Interference]) -> i32 {
        let mut loop_: bool;
        // OCCT L1663: NCollection_List iterator → index walk.
        let ec = ed.geometry();
        let mut sta = ec.parameter_3d(ec.first_parameter());
        let mut end = ec.parameter_3d(ec.last_parameter());
        let tolpar = (end - sta) * 0.01;
        let mut param: f64;
        loop_ = true;
        let mut it = 0usize;

        while it < il.len() && loop_ {
            param = il[it].intersection().parameter();
            if param > end {
                loop_ = false;
            } else {
                if (param - sta).abs() > (param - end).abs() {
                    end = param;
                } else {
                    sta = param;
                }
            }
            it += 1;
        }
        param = 0.5 * (sta + end);
        let mut level: i32 = 0;
        /*TopAbs_State st = */ let _ = self.classify(e, ed, true, &mut level, param);
        loop_ = true;
        it = 0;

        while it < il.len() && loop_ {
            let int_ = &il[it];
            let p = int_.intersection().parameter();
            if p < param - tolpar {
                match int_.transition() {
                    Orientation::Forward => {
                        level -= int_.intersection().level();
                    }
                    Orientation::Reversed => {
                        level += int_.intersection().level();
                    }
                    Orientation::External | Orientation::Internal => {}
                }
            } else if p > param + tolpar {
                loop_ = false;
            } else {
                // OCCT L1725-1728: the #ifdef OCCT_DEBUG "Bad Parameter."
                // print is not translated.
            }
            it += 1;
        }
        if Self::trace_enabled() {
            eprintln!(
                "[TRACE] hiding_start_level: E={} iFace={} -> level={}",
                e, self.i_face, level
            );
        }
        level
    }

    /// OCCT Compare (cxx L1737-1742).
    pub fn compare(&mut self, e: i32, ed: &EdgeData<'a>) -> State {
        let mut level: i32 = 0;
        let parbid: f64 = 0.0;
        self.classify(e, ed, false, &mut level, parbid)
    }

    /// OCCT OrientOutLine (cxx L1746-1877).
    pub fn orient_out_line(&mut self, i: i32, fd: &mut FaceData<'a>) -> bool {
        let _ = i; // avoid compiler warning (OCCT L1748)

        // OCCT L1750: const occ::handle<HLRAlgo_WiresBlock>& wb = FD.Wires();
        // the OCCT handle aliases the FaceData-owned block across the whole
        // exploration (the FD scalar reads interleave with the wb
        // mutations); rcad dereferences the Arc through its raw pointer
        // (the shared-handle mutation follows the FaceData Arc::get_mut
        // precedent; the HLR Data is single-threaded).
        let wb: *mut WiresBlock = match fd.change_wires().as_mut() {
            Some(w) => std::sync::Arc::as_ptr(w) as *const WiresBlock as *mut WiresBlock,
            None => std::ptr::null_mut(),
        };
        let nw = unsafe { &*wb }.nb_wires() as i32;
        let mut iw1: i32;
        let mut ie1: i32;
        let mut ne1: i32;
        let t = *self.my_proj.transformation();
        let ti = *self.my_proj.inverted_transformation();
        let mut inverted = false;
        let mut first_inversion = true;

        iw1 = 1;
        while iw1 <= nw {
            // OCCT L1760: eb1 = wb->Wire(iw1);
            let eb1 = unsafe { &mut *wb }.wire(iw1 as usize);
            ne1 = eb1.nb_edges() as i32;

            ie1 = 1;
            while ie1 <= ne1 {
                self.my_fe = eb1.edge(ie1 as usize) as usize;
                // OCCT L1766: HLRBRep_EdgeData& ed1 = myEData(myFE);
                let ed1: &mut EdgeData<'a> = &mut self.my_e_data[(self.my_fe - 1) as usize];
                if eb1.double(ie1 as usize) || eb1.iso_line(ie1 as usize) || ed1.vertical() {
                    ed1.set_used(true);
                } else {
                    ed1.set_used(false);
                }
                if (eb1.out_line(ie1 as usize) || eb1.internal(ie1 as usize)) && !ed1.vertical() {
                    if Self::trace_enabled() {
                        eprintln!(
                            "[TRACE] orient_out_line: face={} edge={} blockOri={:?} out={} int={}",
                            i,
                            self.my_fe,
                            eb1.orientation(ie1 as usize),
                            eb1.out_line(ie1 as usize),
                            eb1.internal(ie1 as usize),
                        );
                    }
                    let mut p: f64;
                    let mut pu: f64 = 0.0; // OCCT uninitialized; neutral default 0.
                    let mut pv: f64 = 0.0;
                    let r: f64;
                    // OCCT L1778: myFEGeom = &(ed1.ChangeGeometry());
                    self.my_fe_geom = ed1.change_geometry() as *mut Curve<'a>;
                    let ec: &Curve<'a> = ed1.geometry();
                    let vsta = ed1.v_sta();
                    let vend = ed1.v_end();
                    if vsta == 0 && vend == 0 {
                        p = 0.0;
                    } else if vsta == 0 {
                        p = ec.parameter_3d(ec.last_parameter());
                    } else if vend == 0 {
                        p = ec.parameter_3d(ec.first_parameter());
                    } else {
                        p = ec.parameter_3d((ec.last_parameter() + ec.first_parameter()) / 2.0);
                    }
                    // OCCT L1798: HLRBRep_EdgeFaceTool::UVPoint(p, myFEGeom,
                    // iFaceGeom, pu, pv) — the brep/edge context travels
                    // through my_brep and the aligned my_e_map slot (the
                    // session-10 contract ruling).
                    if crate::hlr::brep::edge_face_tool::uv_point(
                        self.my_brep.expect("kernel brep context"),
                        p,
                        self.my_fe_geom,
                        &self.my_e_map[(self.my_fe - 1) as usize],
                        self.i_face_geom,
                        &mut pu,
                        &mut pv,
                    ) {
                        // OCCT L1802-1803: mySLProps.SetParameters(pu, pv); EC.D1(p, Pt, Tg);
                        self.my_s_l_props.set_parameters(pu, pv);
                        let (mut pt, mut tg) = ec.view().d1(p);
                        let v = if self.my_proj.perspective() {
                            let mut p2d = DVec2::ZERO;
                            self.my_proj.project_pnt(pt, &mut p2d);
                            DVec3::new(p2d.x, p2d.y, -self.my_proj.focus())
                        } else {
                            DVec3::new(0.0, 0.0, 1.0) // gp_Dir(gp_Dir::D::NZ)
                        };
                        let v = ti.transform_dir(v);
                        if self.my_s_l_props.is_normal_defined() {
                            let curv = crate::hlr::brep::edge_face_tool::curvature_value(
                                self.i_face_geom,
                                pu,
                                pv,
                                v,
                            );
                            // OCCT L1819: gp_Vec Nm = mySLProps.Normal();
                            let mut nm = self
                                .my_s_l_props
                                .normal()
                                .expect("HLRBRep_Data::OrientOutLine");
                            if curv == 0.0 {
                                // OCCT L1822-1826: the #ifdef OCCT_DEBUG
                                // print is not translated.
                            }
                            if curv > 0.0 {
                                nm = -nm; // Nm.Reverse();
                            }
                            tg = t.transform_vec(tg); // Tg.Transform(T);
                            pt = t.apply(pt); // Pt.Transform(T);
                            nm = t.transform_vec(nm); // Nm.Transform(T);
                            nm = nm.cross(tg); // Nm.Cross(Tg);
                            if tg.length() < GP_RESOLUTION {
                                // OCCT L1838-1842: the #ifdef OCCT_DEBUG
                                // print is not translated.
                            }
                            if self.my_proj.perspective() {
                                r = nm.z * self.my_proj.focus()
                                    - (nm.x * pt.x + nm.y * pt.y + nm.z * pt.z);
                            } else {
                                r = nm.z;
                            }
                            self.my_fe_ori = if r > 0.0 {
                                Orientation::Forward
                            } else {
                                Orientation::Reversed
                            };
                            if !fd.cut() && fd.closed() && first_inversion {
                                if (eb1.orientation(ie1 as usize) == self.my_fe_ori)
                                    != (fd.orientation() == Orientation::Forward)
                                {
                                    first_inversion = false;
                                    inverted = true;
                                }
                            }
                            eb1.set_orientation(ie1 as usize, self.my_fe_ori);
                            if Self::trace_enabled() {
                                eprintln!(
                                    "[TRACE]   -> face={} edge={} setOri={:?} (uv ok, normal ok)",
                                    i, self.my_fe, self.my_fe_ori
                                );
                            }
                        }
                    } else {
                        // OCCT L1866-1870: the #ifdef OCCT_DEBUG "UVPoint
                        // not found, OutLine not Oriented" print is not
                        // translated.
                    }
                    ed1.set_used(true);
                }
                ie1 += 1;
            }
            iw1 += 1;
        }
        inverted
    }

    /// OCCT OrientOthEdge (cxx L1881-1943).
    pub fn orient_oth_edge(&mut self, i: i32, fd: &mut FaceData<'a>) {
        let _ = i; // OCCT L1938: (void)I; in the #else branch.
        let mut p: f64;
        let mut pu: f64 = 0.0; // OCCT uninitialized; neutral default 0.
        let mut pv: f64 = 0.0;
        // OCCT L1884: wb = FD.Wires() (the raw Arc aliasing, see
        // orient_out_line).
        let wb: *mut WiresBlock = match fd.change_wires().as_mut() {
            Some(w) => std::sync::Arc::as_ptr(w) as *const WiresBlock as *mut WiresBlock,
            None => std::ptr::null_mut(),
        };
        let nw = unsafe { &*wb }.nb_wires() as i32;
        let mut iw1: i32;
        let mut ie1: i32;
        let mut ne1: i32;
        let t = *self.my_proj.transformation();

        iw1 = 1;
        while iw1 <= nw {
            // OCCT L1891: eb1 = wb->Wire(iw1);
            let eb1 = unsafe { &mut *wb }.wire(iw1 as usize);
            ne1 = eb1.nb_edges() as i32;

            ie1 = 1;
            while ie1 <= ne1 {
                self.my_fe = eb1.edge(ie1 as usize) as usize;
                self.my_fe_ori = eb1.orientation(ie1 as usize);
                // OCCT L1898: HLRBRep_EdgeData& ed1 = myEData(myFE);
                let ed1: &mut EdgeData<'a> = &mut self.my_e_data[(self.my_fe - 1) as usize];

                if !ed1.used() {
                    ed1.set_used(true);
                    // OCCT L1903: myFEGeom = &(ed1.ChangeGeometry());
                    self.my_fe_geom = ed1.change_geometry() as *mut Curve<'a>;
                    let ec: &Curve<'a> = ed1.geometry();
                    p = ec.parameter_3d((ec.last_parameter() + ec.first_parameter()) / 2.0);
                    // OCCT L1906: HLRBRep_EdgeFaceTool::UVPoint(p, myFEGeom,
                    // iFaceGeom, pu, pv) — the brep/edge context travels
                    // through my_brep and the aligned my_e_map slot (the
                    // session-10 contract ruling).
                    if crate::hlr::brep::edge_face_tool::uv_point(
                        self.my_brep.expect("kernel brep context"),
                        p,
                        self.my_fe_geom,
                        &self.my_e_map[(self.my_fe - 1) as usize],
                        self.i_face_geom,
                        &mut pu,
                        &mut pv,
                    ) {
                        let mut pt = ec.value_3d(p); // OCCT L1908: gp_Pnt Pt = EC.Value3D(p);
                        self.my_s_l_props.set_parameters(pu, pv);
                        if self.my_s_l_props.is_normal_defined() {
                            // OCCT L1912: gp_Vec Nm = mySLProps.Normal();
                            let mut nm = self
                                .my_s_l_props
                                .normal()
                                .expect("HLRBRep_Data::OrientOthEdge");
                            pt = t.apply(pt); // Pt.Transform(T);
                            nm = t.transform_vec(nm); // Nm.Transform(T);
                            // OCCT L1883 declares r at the method top; the
                            // rcad deferred binding lives per edge (the
                            // fresh-uninitialized value of each OCCT pass).
                            let r: f64;
                            if self.my_proj.perspective() {
                                r = nm.z * self.my_proj.focus()
                                    - (nm.x * pt.x + nm.y * pt.y + nm.z * pt.z);
                            } else {
                                r = nm.z;
                            }
                            if r < 0.0 {
                                self.my_fe_ori = top_abs_reverse(self.my_fe_ori); // TopAbs::Reverse
                                eb1.set_orientation(ie1 as usize, self.my_fe_ori);
                            }
                        }
                    }
                    // OCCT L1930-1936: the #ifdef OCCT_DEBUG "UVPoint not
                    // found, Edge not Oriented" print is not translated.
                }
                ie1 += 1;
            }
            iw1 += 1;
        }
    }

    /// OCCT Classify (cxx L1993-2268).
    pub fn classify(
        &mut self,
        e: i32,
        ed: &EdgeData<'a>,
        level_flag: bool,
        level: &mut i32,
        param: f64,
    ) -> State {
        let _ = e; // avoid compiler warning (OCCT L1999)

        counters::NB_CLASSIFICATION.with(|c| c.set(c.get() + 1));
        let mut vert_min = MinMaxIndices::default();
        let mut vert_max = MinMaxIndices::default();
        let mut min_max_vert = MinMaxIndices::default();
        let mut tot_min = [0.0f64; 16]; // OCCT uninitialized; InitMinMax fills before use.
        let mut tot_max = [0.0f64; 16];

        let mut i: i32;
        *level = 0;
        let mut state = State::Out;
        //  bool rej = false;
        let ec: &Curve<'a> = ed.geometry();
        let mut sta: f64;
        let mut xsta: f64 = 0.0; // OCCT uninitialized; neutral default 0.
        let mut ysta: f64 = 0.0;
        let mut zsta: f64 = 0.0;
        let mut end: f64 = 0.0;
        let mut xend: f64 = 0.0;
        let mut yend: f64 = 0.0;
        let mut zend: f64 = 0.0;
        let tol = ed.tolerance() as f64;

        if level_flag {
            sta = param;
            self.my_proj
                .project_xyz(ec.value_3d(sta), &mut xsta, &mut ysta, &mut zsta);

            //-- les rejections sont faites dans l intersecteur a moindre frais
            //-- puisque la surface sera chargee
            HLRAlgo::init_min_max(PRECISION_INFINITE, &mut tot_min, &mut tot_max);
            HLRAlgo::update_min_max(xsta, ysta, zsta, &mut tot_min, &mut tot_max);
            HLRAlgo::enlarge_min_max(tol, &mut tot_min, &mut tot_max);
            reject1(
                &self.my_deca,
                &tot_min,
                &tot_max,
                &self.my_sur_d,
                &mut vert_min,
                &mut vert_max,
            );

            HLRAlgo::encode_min_max(&vert_min, &vert_max, &mut min_max_vert);
            let i_face_min_max = unsafe { &*self.i_face_min_max };
            {
                // TEMPORARY diagnostic: find the rejecting dimension.
                let mut rej_k: i32 = -1;
                for k in 0..8usize {
                    if ((i_face_min_max.max[k].wrapping_sub(min_max_vert.min[k])) & REJECT_MASK)
                        != 0
                        || ((min_max_vert.max[k].wrapping_sub(i_face_min_max.min[k])) & REJECT_MASK)
                            != 0
                    {
                        rej_k = k as i32;
                        break;
                    }
                }
                if Self::trace_enabled() && rej_k >= 0 {
                    let p3 = ec.value_3d(sta);
                    eprintln!(
                        "[TRACE] classify REJECT level_flag: E={} iFace={} param={:.6} k={} p3d=({:.4},{:.4},{:.4}) projz={:.4} faceMin={:?} faceMax={:?} ptMin={:?} ptMax={:?}",
                        e,
                        self.i_face,
                        param,
                        rej_k,
                        p3.x,
                        p3.y,
                        p3.z,
                        zsta,
                        &i_face_min_max.min[..],
                        &i_face_min_max.max[..],
                        &min_max_vert.min[..],
                        &min_max_vert.max[..],
                    );
                }
            }
            if ((i_face_min_max.max[0].wrapping_sub(min_max_vert.min[0])) & REJECT_MASK) != 0
                || ((min_max_vert.max[0].wrapping_sub(i_face_min_max.min[0])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[1].wrapping_sub(min_max_vert.min[1])) & REJECT_MASK) != 0
                || ((min_max_vert.max[1].wrapping_sub(i_face_min_max.min[1])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[2].wrapping_sub(min_max_vert.min[2])) & REJECT_MASK) != 0
                || ((min_max_vert.max[2].wrapping_sub(i_face_min_max.min[2])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[3].wrapping_sub(min_max_vert.min[3])) & REJECT_MASK) != 0
                || ((min_max_vert.max[3].wrapping_sub(i_face_min_max.min[3])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[4].wrapping_sub(min_max_vert.min[4])) & REJECT_MASK) != 0
                || ((min_max_vert.max[4].wrapping_sub(i_face_min_max.min[4])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[5].wrapping_sub(min_max_vert.min[5])) & REJECT_MASK) != 0
                || ((min_max_vert.max[5].wrapping_sub(i_face_min_max.min[5])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[6].wrapping_sub(min_max_vert.min[6])) & REJECT_MASK) != 0
                || ((min_max_vert.max[6].wrapping_sub(i_face_min_max.min[6])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[7].wrapping_sub(min_max_vert.min[7])) & REJECT_MASK) != 0
                || ((min_max_vert.max[7].wrapping_sub(i_face_min_max.min[7])) & REJECT_MASK) != 0
            {
                //-- Rejection en Z
                return state;
            }
        } else {
            sta = ec.parameter_3d(ec.first_parameter());
            self.my_proj
                .project_xyz(ec.value_3d(sta), &mut xsta, &mut ysta, &mut zsta);

            //-- les rejections sont faites dans l intersecteur a moindre frais
            //-- puisque la surface sera chargee
            HLRAlgo::init_min_max(PRECISION_INFINITE, &mut tot_min, &mut tot_max);
            HLRAlgo::update_min_max(xsta, ysta, zsta, &mut tot_min, &mut tot_max);
            HLRAlgo::enlarge_min_max(tol, &mut tot_min, &mut tot_max);

            reject1(
                &self.my_deca,
                &tot_min,
                &tot_max,
                &self.my_sur_d,
                &mut vert_min,
                &mut vert_max,
            );

            HLRAlgo::encode_min_max(&vert_min, &vert_max, &mut min_max_vert);
            let i_face_min_max = unsafe { &*self.i_face_min_max };
            if ((i_face_min_max.max[0].wrapping_sub(min_max_vert.min[0])) & REJECT_MASK) != 0
                || ((min_max_vert.max[0].wrapping_sub(i_face_min_max.min[0])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[1].wrapping_sub(min_max_vert.min[1])) & REJECT_MASK) != 0
                || ((min_max_vert.max[1].wrapping_sub(i_face_min_max.min[1])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[2].wrapping_sub(min_max_vert.min[2])) & REJECT_MASK) != 0
                || ((min_max_vert.max[2].wrapping_sub(i_face_min_max.min[2])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[3].wrapping_sub(min_max_vert.min[3])) & REJECT_MASK) != 0
                || ((min_max_vert.max[3].wrapping_sub(i_face_min_max.min[3])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[4].wrapping_sub(min_max_vert.min[4])) & REJECT_MASK) != 0
                || ((min_max_vert.max[4].wrapping_sub(i_face_min_max.min[4])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[5].wrapping_sub(min_max_vert.min[5])) & REJECT_MASK) != 0
                || ((min_max_vert.max[5].wrapping_sub(i_face_min_max.min[5])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[6].wrapping_sub(min_max_vert.min[6])) & REJECT_MASK) != 0
                || ((min_max_vert.max[6].wrapping_sub(i_face_min_max.min[6])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[7].wrapping_sub(min_max_vert.min[7])) & REJECT_MASK) != 0
                || ((min_max_vert.max[7].wrapping_sub(i_face_min_max.min[7])) & REJECT_MASK) != 0
            {
                //-- Rejection en Z
                return state;
            }
            end = ec.parameter_3d(ec.last_parameter());
            self.my_proj
                .project_xyz(ec.value_3d(end), &mut xend, &mut yend, &mut zend);

            HLRAlgo::init_min_max(PRECISION_INFINITE, &mut tot_min, &mut tot_max);
            HLRAlgo::update_min_max(xend, yend, zend, &mut tot_min, &mut tot_max);
            HLRAlgo::enlarge_min_max(tol, &mut tot_min, &mut tot_max);

            reject1(
                &self.my_deca,
                &tot_min,
                &tot_max,
                &self.my_sur_d,
                &mut vert_min,
                &mut vert_max,
            );

            HLRAlgo::encode_min_max(&vert_min, &vert_max, &mut min_max_vert);
            let i_face_min_max = unsafe { &*self.i_face_min_max };
            if ((i_face_min_max.max[0].wrapping_sub(min_max_vert.min[0])) & REJECT_MASK) != 0
                || ((min_max_vert.max[0].wrapping_sub(i_face_min_max.min[0])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[1].wrapping_sub(min_max_vert.min[1])) & REJECT_MASK) != 0
                || ((min_max_vert.max[1].wrapping_sub(i_face_min_max.min[1])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[2].wrapping_sub(min_max_vert.min[2])) & REJECT_MASK) != 0
                || ((min_max_vert.max[2].wrapping_sub(i_face_min_max.min[2])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[3].wrapping_sub(min_max_vert.min[3])) & REJECT_MASK) != 0
                || ((min_max_vert.max[3].wrapping_sub(i_face_min_max.min[3])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[4].wrapping_sub(min_max_vert.min[4])) & REJECT_MASK) != 0
                || ((min_max_vert.max[4].wrapping_sub(i_face_min_max.min[4])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[5].wrapping_sub(min_max_vert.min[5])) & REJECT_MASK) != 0
                || ((min_max_vert.max[5].wrapping_sub(i_face_min_max.min[5])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[6].wrapping_sub(min_max_vert.min[6])) & REJECT_MASK) != 0
                || ((min_max_vert.max[6].wrapping_sub(i_face_min_max.min[6])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[7].wrapping_sub(min_max_vert.min[7])) & REJECT_MASK) != 0
                || ((min_max_vert.max[7].wrapping_sub(i_face_min_max.min[7])) & REJECT_MASK) != 0
            {
                //-- Rejection en Z
                return state;
            }
            sta = 0.4 * sta + 0.6 * end; // dangerous if it is the middle
            self.my_proj
                .project_xyz(ec.value_3d(sta), &mut xsta, &mut ysta, &mut zsta);

            //-- les rejections sont faites dans l intersecteur a moindre frais
            //-- puisque la surface sera chargee
            HLRAlgo::init_min_max(PRECISION_INFINITE, &mut tot_min, &mut tot_max);
            HLRAlgo::update_min_max(xsta, ysta, zsta, &mut tot_min, &mut tot_max);
            HLRAlgo::enlarge_min_max(tol, &mut tot_min, &mut tot_max);
            reject1(
                &self.my_deca,
                &tot_min,
                &tot_max,
                &self.my_sur_d,
                &mut vert_min,
                &mut vert_max,
            );

            HLRAlgo::encode_min_max(&vert_min, &vert_max, &mut min_max_vert);
            // OCCT L2116-2157: the commented-out OCCT_DEBUG printf block is
            // an OCCT comment and is not translated.
            let i_face_min_max = unsafe { &*self.i_face_min_max };
            if ((i_face_min_max.max[0].wrapping_sub(min_max_vert.min[0])) & REJECT_MASK) != 0
                || ((min_max_vert.max[0].wrapping_sub(i_face_min_max.min[0])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[1].wrapping_sub(min_max_vert.min[1])) & REJECT_MASK) != 0
                || ((min_max_vert.max[1].wrapping_sub(i_face_min_max.min[1])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[2].wrapping_sub(min_max_vert.min[2])) & REJECT_MASK) != 0
                || ((min_max_vert.max[2].wrapping_sub(i_face_min_max.min[2])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[3].wrapping_sub(min_max_vert.min[3])) & REJECT_MASK) != 0
                || ((min_max_vert.max[3].wrapping_sub(i_face_min_max.min[3])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[4].wrapping_sub(min_max_vert.min[4])) & REJECT_MASK) != 0
                || ((min_max_vert.max[4].wrapping_sub(i_face_min_max.min[4])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[5].wrapping_sub(min_max_vert.min[5])) & REJECT_MASK) != 0
                || ((min_max_vert.max[5].wrapping_sub(i_face_min_max.min[5])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[6].wrapping_sub(min_max_vert.min[6])) & REJECT_MASK) != 0
                || ((min_max_vert.max[6].wrapping_sub(i_face_min_max.min[6])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[7].wrapping_sub(min_max_vert.min[7])) & REJECT_MASK) != 0
                || ((min_max_vert.max[7].wrapping_sub(i_face_min_max.min[7])) & REJECT_MASK) != 0
            {
                //-- Rejection en Z
                return state;
            }
        }

        counters::NB_CAL3_INTERSECTION.with(|c| c.set(c.get() + 1));
        let mut plim: rcad_kernel::geom::Point3;
        let psta: DVec2;
        psta = ec.value(sta);
        plim = ec.value_3d(sta);

        // OCCT L2185-2190: static int aff = 0; ... — the dead debug print
        // block (aff is always 0).
        let aff = 0;
        if aff != 0 {
            let mut nump1 = 0;
            nump1 += 1;
            let _ = nump1; // printf("\npoint PNR%d  %g %g %g", ...)
        }

        let lin_l = self.my_proj.shoot(psta.x, psta.y); // OCCT L2192: gp_Lin L
        // OCCT gp_Lin → the rcad kernel line (the Lin/Line3 field mapping).
        let l = Line3 {
            origin: lin_l.pos,
            direction: lin_l.dir,
        };
        let mut w_lim = elclib::line_parameter(&l, plim);
        self.my_intersector.perform_line(&l, w_lim);
        let trace_cl = Self::trace_enabled();
        let trace_hits: &mut Vec<(f64, f64, f64, bool)> = &mut Vec::new();
        let _ = trace_hits;
        if trace_cl {
            eprintln!(
                "[TRACE] classify: E={} iFace={} param={:.6} plim=({:.4},{:.4},{:.4}) w_lim0={:.6} done={} nbpts={}",
                e,
                self.i_face,
                param,
                plim.x,
                plim.y,
                plim.z,
                w_lim,
                self.my_intersector.is_done(),
                if self.my_intersector.is_done() { self.my_intersector.nb_points() } else { 0 },
            );
        }
        if self.my_intersector.is_done() {
            let nb_points = self.my_intersector.nb_points();
            if nb_points > 0 {
                let mut tol_z = self.my_big_size * 0.000001;
                if self.i_face_test {
                    if !self.my_le_out_line && !self.my_le_internal {
                        tol_z = self.my_big_size * 0.001;
                    } else {
                        tol_z = self.my_big_size * 0.01;
                    }
                }
                w_lim -= tol_z;
                let mut period_u: f64;
                let mut period_v: f64;
                let mut u_min: f64 = 0.0;
                let mut u_max: f64 = 0.0;
                let mut v_min: f64 = 0.0;
                let mut v_max: f64 = 0.0;
                let i_face_geom = unsafe { &*self.i_face_geom };
                if i_face_geom.is_u_periodic() {
                    period_u = i_face_geom.u_period();
                    u_min = i_face_geom.first_u_parameter();
                    u_max = i_face_geom.last_u_parameter();
                } else {
                    period_u = 0.0;
                }
                if i_face_geom.is_v_periodic() {
                    period_v = i_face_geom.v_period();
                    v_min = i_face_geom.first_v_parameter();
                    v_max = i_face_geom.last_v_parameter();
                } else {
                    period_v = 0.0;
                }

                i = 1;
                while i <= nb_points as i32 {
                    // OCCT L2240: myIntersector.CSPoint(i).Values(PInter, u, v, w, Tr);
                    // (the PInter/Tr out values are not read in OCCT either).
                    let (_p_inter, mut u, mut v, w, _tr) =
                        self.my_intersector.cs_point(i as usize).values();
                    let mut in_domain = false;
                    if w < w_lim {
                        let mut a_dummy_shift: f64 = 0.0; // OCCT uninitialized; neutral default 0.
                        if period_u > 0.0 {
                            // OCCT L2246: GeomInt::AdjustPeriodic(u, UMin, UMax, PeriodU, u, aDummyShift);
                            let _ = geom_int_adjust_periodic(
                                u, u_min, u_max, period_u, &mut u, &mut a_dummy_shift, 0.0,
                            );
                        }
                        if period_v > 0.0 {
                            // OCCT L2250: GeomInt::AdjustPeriodic(v, VMin, VMax, PeriodV, v, aDummyShift);
                            let _ = geom_int_adjust_periodic(
                                v, v_min, v_max, period_v, &mut v, &mut a_dummy_shift, 0.0,
                            );
                        }

                        let pnt2d = DVec2::new(u, v);
                        // OCCT L2254: myClassifier->Classify(pnt2d,
                        // Precision::PConfusion()) — the default
                        // RecadreOnPeriodic = Standard_False.
                        let cl = self.my_classifier.classify(pnt2d, PRECISION_P_CONFUSION, false);
                        in_domain = cl != State::Out;
                        if trace_cl {
                            trace_hits.push((u, v, w, in_domain));
                        }
                        if in_domain {
                            state = State::In;
                            *level += 1;
                            if !level_flag {
                                return state;
                            }
                        }
                    }
                    i += 1;
                }
                if trace_cl {
                    eprintln!("[TRACE] classify hits (u,v,w,in): {:?}", trace_hits);
                }
            }
        }
        state
    }

    /// OCCT SimplClassify (cxx L2272-2323).
    pub fn simpl_classify(
        &mut self,
        _e: i32,
        ed: &EdgeData<'a>,
        nbp: i32,
        p1: f64,
        p2: f64,
    ) -> State {
        counters::NB_CLASSIFICATION.with(|c| c.set(c.get() + 1));
        let mut vert_min = MinMaxIndices::default();
        let mut vert_max = MinMaxIndices::default();
        let mut min_max_vert = MinMaxIndices::default();
        let mut tot_min = [0.0f64; 16]; // OCCT uninitialized; InitMinMax fills before use.
        let mut tot_max = [0.0f64; 16];

        let mut i: i32;
        let mut state = State::In;
        //  bool rej = false;
        let ec: &Curve<'a> = ed.geometry();
        let mut sta: f64;
        let mut xsta: f64 = 0.0; // OCCT uninitialized; neutral default 0.
        let mut ysta: f64 = 0.0;
        let mut zsta: f64 = 0.0;
        let mut dp: f64;
        let tol = ed.tolerance() as f64;

        dp = (p2 - p1) / (nbp as f64 + 1.0);

        // OCCT L2291: for (sta = p1 + dp, i = 1; i <= Nbp; ++i, sta += dp)
        sta = p1 + dp;
        i = 1;
        while i <= nbp {
            self.my_proj
                .project_xyz(ec.value_3d(sta), &mut xsta, &mut ysta, &mut zsta);

            //-- les rejections sont faites dans l intersecteur a moindre frais
            //-- puisque la surface sera chargee
            HLRAlgo::init_min_max(PRECISION_INFINITE, &mut tot_min, &mut tot_max);
            HLRAlgo::update_min_max(xsta, ysta, zsta, &mut tot_min, &mut tot_max);
            HLRAlgo::enlarge_min_max(tol, &mut tot_min, &mut tot_max);
            reject1(
                &self.my_deca,
                &tot_min,
                &tot_max,
                &self.my_sur_d,
                &mut vert_min,
                &mut vert_max,
            );

            HLRAlgo::encode_min_max(&vert_min, &vert_max, &mut min_max_vert);
            let i_face_min_max = unsafe { &*self.i_face_min_max };
            if ((i_face_min_max.max[0].wrapping_sub(min_max_vert.min[0])) & REJECT_MASK) != 0
                || ((min_max_vert.max[0].wrapping_sub(i_face_min_max.min[0])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[1].wrapping_sub(min_max_vert.min[1])) & REJECT_MASK) != 0
                || ((min_max_vert.max[1].wrapping_sub(i_face_min_max.min[1])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[2].wrapping_sub(min_max_vert.min[2])) & REJECT_MASK) != 0
                || ((min_max_vert.max[2].wrapping_sub(i_face_min_max.min[2])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[3].wrapping_sub(min_max_vert.min[3])) & REJECT_MASK) != 0
                || ((min_max_vert.max[3].wrapping_sub(i_face_min_max.min[3])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[4].wrapping_sub(min_max_vert.min[4])) & REJECT_MASK) != 0
                || ((min_max_vert.max[4].wrapping_sub(i_face_min_max.min[4])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[5].wrapping_sub(min_max_vert.min[5])) & REJECT_MASK) != 0
                || ((min_max_vert.max[5].wrapping_sub(i_face_min_max.min[5])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[6].wrapping_sub(min_max_vert.min[6])) & REJECT_MASK) != 0
                || ((min_max_vert.max[6].wrapping_sub(i_face_min_max.min[6])) & REJECT_MASK) != 0
                || ((i_face_min_max.max[7].wrapping_sub(min_max_vert.min[7])) & REJECT_MASK) != 0
                || ((min_max_vert.max[7].wrapping_sub(i_face_min_max.min[7])) & REJECT_MASK) != 0
            {
                //-- Rejection en Z
                return State::Out;
            }
            i += 1;
            sta += dp;
        }
        state
    }

    /// OCCT RejectedPoint (cxx L2330-2593) — build an interference if non
    /// Rejected intersection point.
    pub fn rejected_point(
        &mut self,
        p_inter: &IntersectionPoint,
        bound_ori: Orientation,
        num_seg: i32,
    ) -> bool {
        let mut ind: i32 = 0; // OCCT L2334: int Ind = 0;
        let mut decal: i32;
        let mut p1: f64;
        let mut p2: f64;
        let mut dz: f64;
        let mut t1: f32 = 0.0; // OCCT uninitialized float; neutral default 0.
        let mut t2: f32 = 0.0;
        let st: State;
        let mut orie = Orientation::Forward; // TopAbs_Orientation Orie = TopAbs_FORWARD;
        let mut or2 = Orientation::Internal; // TopAbs_Orientation Or2 = TopAbs_INTERNAL;
        let mut inverted = false;
        let tr1: crate::geomalgo::int_res2d::Transition;
        let tr2: crate::geomalgo::int_res2d::Transition;
        let tol_z = self.my_big_size * 0.00001;

        let my_le_geom = unsafe { &*self.my_le_geom };
        let my_fe_geom = unsafe { &*self.my_fe_geom };
        p1 = my_le_geom.parameter_3d(p_inter.param_on_first());
        p2 = my_fe_geom.parameter_3d(p_inter.param_on_second());
        dz = my_le_geom.z(p1) - my_fe_geom.z(p2);

        if self.my_le == self.my_fe {
            // auto intersection can be inverted
            if dz >= tol_z {
                inverted = true;
                let p = p1;
                p1 = p2;
                p2 = p;
                dz = -dz;
            }
        }

        if dz >= tol_z {
            self.my_above_intf = true;
            return true;
        }
        self.my_above_intf = false;
        st = if dz <= -tol_z { State::In } else { State::On };

        if inverted {
            tr1 = *p_inter.transition_of_second();
            tr2 = *p_inter.transition_of_first();
        } else {
            tr1 = *p_inter.transition_of_first();
            tr2 = *p_inter.transition_of_second();
        }

        if self.i_face_test {
            if self.my_le == self.my_fe {
                if st == State::In {
                    unsafe { &mut *self.my_le_data }.set_simple(false);
                }
            } else {
                if self.my_same_vertex {
                    if (st == State::On)
                        || (tr1.position_on_curve() != Position::Middle)
                        || (tr2.position_on_curve() != Position::Middle)
                    {
                        return true;
                    }
                }
            }
            if st == State::In {
                self.i_face_smpl = false;
            }
        }

        match tr1.transition_type() {
            // compute the transition
            TypeTrans::In => {
                orie = if self.my_fe_ori == Orientation::Reversed {
                    Orientation::Reversed
                } else {
                    Orientation::Forward
                };
            }
            TypeTrans::Out => {
                orie = if self.my_fe_ori == Orientation::Reversed {
                    Orientation::Forward
                } else {
                    Orientation::Reversed
                };
            }
            TypeTrans::Touch => {
                match tr1.situation() {
                    Situation::Inside => {
                        orie = if self.my_fe_ori == Orientation::Reversed {
                            Orientation::External
                        } else {
                            Orientation::Internal
                        };
                    }
                    Situation::Outside => {
                        orie = if self.my_fe_ori == Orientation::Reversed {
                            Orientation::Internal
                        } else {
                            Orientation::External
                        };
                    }
                    Situation::Unknown => {
                        return true;
                    }
                }
            }
            TypeTrans::Undecided => {
                return true;
            }
        }

        if self.i_face_back {
            orie = top_abs_complement(orie); // change the transition
        }
        let mut ori = Orientation::Forward;
        match tr1.position_on_curve() {
            Position::Head => {
                ori = Orientation::Forward;
            }
            Position::Middle => {
                ori = Orientation::Internal;
            }
            Position::End => {
                ori = Orientation::Reversed;
            }
        }

        if st != State::Out {
            if tr2.position_on_curve() != Position::Middle {
                // correction de la transition  sur myFE
                if self.my_same_vertex {
                    return true; // si intersection a une extremite verticale !
                }

                let mut douteux = false;
                let psav = p2;
                let mut ptsav = DVec2::ZERO; // OCCT uninitialized; neutral default 0.
                let mut tgsav = DVec2::ZERO;
                let mut nmsav = DVec2::ZERO;
                if tr2.position_on_curve() == Position::Head {
                    let my_fe_data = unsafe { &mut *self.my_fe_data };
                    ind = my_fe_data.v_sta();
                    or2 = Orientation::Forward;
                    adjust_parameter(my_fe_data, true, &mut p2, &mut t2);
                    if unsafe { &*self.my_fe_data }.ver_at_sta() {
                        douteux = true;
                        my_fe_geom.d2_2d(psav, &mut ptsav, &mut tgsav, &mut nmsav);
                        if tgsav.length_squared() <= DERIVEE_PREMIERE_NULLE {
                            tgsav = nmsav;
                        }
                    }
                } else {
                    let my_fe_data = unsafe { &mut *self.my_fe_data };
                    ind = my_fe_data.v_end();
                    or2 = Orientation::Reversed;
                    adjust_parameter(my_fe_data, false, &mut p2, &mut t2);
                    if unsafe { &*self.my_fe_data }.ver_at_end() {
                        douteux = true;
                        my_fe_geom.d2_2d(psav, &mut ptsav, &mut tgsav, &mut nmsav);
                        if tgsav.length_squared() <= DERIVEE_PREMIERE_NULLE {
                            tgsav = nmsav;
                        }
                    }
                }
                let mut tgfe = DVec2::ZERO; // gp_Vec2d TgFE;
                my_fe_geom.d1_2d(p2, &mut ptsav, &mut tgfe);
                if douteux {
                    // OCCT L2497: TgFE.XY().Dot(Tgsav.XY())
                    if tgfe.dot(tgsav) < 0.0 {
                        if orie == Orientation::Forward {
                            orie = Orientation::Reversed;
                        } else if orie == Orientation::Reversed {
                            orie = Orientation::Forward;
                        }
                    }
                }
                self.my_intf.change_boundary().set_2d(self.my_fe as i32, p2);
            }
            if ori != Orientation::Internal {
                // correction de la transition  sur myLE
                let mut douteux = false; // si intersection a une extremite verticale !
                let psav = p1;
                let mut ptsav = DVec2::ZERO; // OCCT uninitialized; neutral default 0.
                let mut tgsav = DVec2::ZERO;
                let mut nmsav = DVec2::ZERO;
                if ori == Orientation::Forward {
                    adjust_parameter(unsafe { &mut *self.my_le_data }, true, &mut p1, &mut t1);
                    if unsafe { &*self.my_le_data }.ver_at_sta() {
                        douteux = true;
                        my_le_geom.d2_2d(psav, &mut ptsav, &mut tgsav, &mut nmsav);
                        if tgsav.length_squared() <= DERIVEE_PREMIERE_NULLE {
                            tgsav = nmsav;
                        }
                    }
                } else {
                    adjust_parameter(unsafe { &mut *self.my_le_data }, false, &mut p1, &mut t1);
                    if unsafe { &*self.my_le_data }.ver_at_end() {
                        douteux = true;
                        my_le_geom.d2_2d(psav, &mut ptsav, &mut tgsav, &mut nmsav);
                        if tgsav.length_squared() <= DERIVEE_PREMIERE_NULLE {
                            tgsav = nmsav;
                        }
                    }
                }
                if douteux {
                    let mut tgle = DVec2::ZERO; // gp_Vec2d TgLE;
                    my_le_geom.d1_2d(p1, &mut ptsav, &mut tgle);
                    // OCCT L2547: TgLE.XY().Dot(Tgsav.XY())
                    if tgle.dot(tgsav) < 0.0 {
                        if orie == Orientation::Forward {
                            orie = Orientation::Reversed;
                        } else if orie == Orientation::Reversed {
                            orie = Orientation::Forward;
                        }
                    }
                }
            }
            if st == State::On {
                let mut stbef = State::In; // OCCT uninitialized; neutral default (TopAbs 0).
                let mut staft = State::In;
                self.edge_state(p1, p2, &mut stbef, &mut staft);
                self.my_intf.change_boundary().set_state_3d(stbef, staft);
            }
        }

        if self.my_fe_internal {
            decal = 2;
        } else {
            decal = 1;
            if st == State::In && ori == Orientation::Forward && orie == Orientation::Forward {
                decal = 0;
            }
        }
        // OCCT L2580: HLRAlgo_Intersection& inter = myIntf.ChangeIntersection();
        let inter = self.my_intf.change_intersection();
        inter.set_orientation(ori);
        inter.set_level(decal);
        inter.set_seg_index(num_seg);
        inter.set_index(ind);
        inter.set_parameter(p1);
        inter.set_tolerance(self.my_le_tol);
        inter.set_state(st);
        self.my_intf.set_orientation(or2);
        self.my_intf.set_transition(orie);
        self.my_intf.set_boundary_transition(bound_ori);
        self.my_intf.change_boundary().set_2d(self.my_fe as i32, p2);
        false
    }

    /// OCCT SameVertex (cxx L2597-2651).
    pub fn same_vertex(&mut self, h1: bool, h2: bool) -> bool {
        let v1: i32;
        let v2: i32;
        if h1 {
            v1 = unsafe { &*self.my_le_data }.v_sta();
        } else {
            v1 = unsafe { &*self.my_le_data }.v_end();
        }
        if h2 {
            v2 = unsafe { &*self.my_fe_data }.v_sta();
        } else {
            v2 = unsafe { &*self.my_fe_data }.v_end();
        }
        let same_v = v1 == v2;
        if same_v {
            self.my_intersected = true; // compute the intersections
            if (self.my_le_type == GeomAbsCurveType::Line
                || self.my_le_type == GeomAbsCurveType::Circle
                || self.my_le_type == GeomAbsCurveType::Ellipse)
                && (self.my_fe_type == GeomAbsCurveType::Line
                    || self.my_fe_type == GeomAbsCurveType::Circle
                    || self.my_fe_type == GeomAbsCurveType::Ellipse)
            {
                self.my_intersected = false; // no other intersection
            }

            let mut other_case = true;

            if (h1 && unsafe { &*self.my_le_data }.out_lv_sta())
                || (!h1 && unsafe { &*self.my_le_data }.out_lv_end())
            {
                if self.i_face_test || self.my_le_internal {
                    other_case = false;
                }
            } else if self.i_face_test {
                other_case = false;
            }

            if other_case {
                if (h1 && unsafe { &*self.my_le_data }.cut_at_sta())
                    || (!h1 && unsafe { &*self.my_le_data }.cut_at_end())
                {
                    self.my_intersected = false; // two connected OutLines do not
                } // intersect themselves.
            }
        }
        same_v
    }

    /// OCCT IsBadFace (cxx L2655-2683).
    pub fn is_bad_face(&self) -> bool {
        if !self.i_face_geom.is_null() {
            // check for garbage data - if periodic then bounds must not exceed period
            let p_geom = unsafe { &*self.i_face_geom };
            if p_geom.is_u_periodic() {
                let a_period = p_geom.u_period();
                let a_min = p_geom.first_u_parameter();
                let a_max = p_geom.last_u_parameter();
                if a_period * 2.0 < a_max - a_min {
                    return true;
                }
            }
            if p_geom.is_v_periodic() {
                let a_period = p_geom.v_period();
                let a_min = p_geom.first_v_parameter();
                let a_max = p_geom.last_v_parameter();
                if a_period * 2.0 < a_max - a_min {
                    return true;
                }
            }
        }
        false
    }
}

