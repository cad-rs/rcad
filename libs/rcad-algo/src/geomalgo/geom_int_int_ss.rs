// OCCT GeomInt_IntSS (TKGeomAlgo/GeomInt) — 1:1 Rust translation.
//
// OCCT sources:
//   GeomInt_IntSS.hxx L38-206  (class, members, method inventory)
//   GeomInt_IntSS.lxx L22-132  (ctor x2, Perform(HS1,HS2,...) x2, IsDone,
//                               TolReached2d/3d, NbLines, NbBoundaries,
//                               NbPoints, Point)
//   GeomInt_IntSS.cxx L24-227  (Perform x2, InternalPerform, Line, Boundary,
//                               Pnt2d, HasLineOnS1/S2, LineOnS1/S2)
//   GeomInt_IntSS_1.cxx        (MakeCurve / TreatRLine / BuildPCurves /
//                               TrimILineOnSurfBoundaries / MakeBSpline x2 —
//                               see geom_int_int_ss_1.rs)
//
// Architecture differences (Rust vs C++):
//   - `occ::handle<GeomAdaptor_Surface> myHS1/myHS2` map to
//     [`GeomSurfaceAdapter`] (the landed 1:1 `GeomAdaptor_Surface`), whose
//     `SurfaceAdapter` trait is `Adaptor3d_HSurfaceTool`.  `Option<...>` is the
//     rcad encoding of the OCCT null handle.
//   - `occ::handle<Adaptor3d_TopolTool> dom1/dom2` are rebuilt on demand from
//     the adaptor surface they are built over (see the same note in
//     `geom_int_line_constructor.rs`): rcad's `TopolTool` borrows its adaptor.
//   - `NCollection_Sequence<handle<Geom_Curve>> sline` and the two pcurve
//     sequences map to `Vec<Option<Curve3>>` / `Vec<Option<Curve2d>>` — OCCT
//     appends a NULL handle `H1` on the "no approximation" arms, which is the
//     `None` entry here.
//   - OCCT sequences are 1-based; `sline(Index + myNbrestr)` becomes
//     `sline[Index + my_nbrestr - 1]`.

use glam::DVec2;

use rcad_kernel::geom::{Curve2d, Curve3};

use super::geom_int_int_ss_1::{
    adapter_identity, define_uv_max_step, int_patch_intersection_prepare_surfaces,
    surface3_identity,
};
use super::geom_int_line_constructor::GeomIntLineConstructor;
use crate::geomalgo::int_patch::IntPatchIntersection;
use crate::hlr::contap::geom_tool::GeomTool;
use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};
use crate::topalgo::adaptor3d::topol_tool::TopolTool;

/// OCCT GeomInt_IntSS (hxx L38-206).
pub struct GeomIntIntSS {
    // OCCT hxx L194: IntPatch_Intersection myIntersector
    my_intersector: IntPatchIntersection,
    // OCCT hxx L195: GeomInt_LineConstructor myLConstruct
    my_l_construct: GeomIntLineConstructor,
    // OCCT hxx L196-197: myHS1 / myHS2
    my_hs1: Option<GeomSurfaceAdapter>,
    my_hs2: Option<GeomSurfaceAdapter>,
    // OCCT hxx L198: int myNbrestr
    my_nbrestr: usize,
    // OCCT hxx L199-201: sline / slineS1 / slineS2
    sline: Vec<Option<Curve3>>,
    sline_s1: Vec<Option<Curve2d>>,
    sline_s2: Vec<Option<Curve2d>>,
    // OCCT hxx L202-205
    my_tol_reached_2d: f64,
    my_tol_reached_3d: f64,
    my_tol_check: f64,
    my_tol_ang_check: f64,
}

impl Default for GeomIntIntSS {
    fn default() -> Self {
        Self::new()
    }
}

impl GeomIntIntSS {
    // ========================================================================
    // OCCT GeomInt_IntSS.lxx L22-29: GeomInt_IntSS()
    // ========================================================================

    /// OCCT GeomInt_IntSS() (lxx L22-29): myNbrestr(0), myTolReached2d(0.0),
    /// myTolReached3d(0.0), myTolCheck(1.e-7), myTolAngCheck(1.e-6).
    pub fn new() -> Self {
        GeomIntIntSS {
            my_intersector: IntPatchIntersection::new(),
            my_l_construct: GeomIntLineConstructor::new(),
            my_hs1: None,
            my_hs2: None,
            my_nbrestr: 0,
            sline: Vec::new(),
            sline_s1: Vec::new(),
            sline_s2: Vec::new(),
            my_tol_reached_2d: 0.0,
            my_tol_reached_3d: 0.0,
            my_tol_check: 1.0e-7,
            my_tol_ang_check: 1.0e-6,
        }
    }

    /// OCCT GeomInt_IntSS(S1, S2, Tol, Approx, ApproxS1, ApproxS2) (lxx
    /// L33-46): the member initializer list + Perform.
    pub fn new_with_surfaces(
        s1: &Surface3Like,
        s2: &Surface3Like,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) -> Self {
        let mut this = GeomIntIntSS::new();
        this.perform_surfaces(s1, s2, tol, approx, approx_s1, approx_s2);
        this
    }

    /// OCCT GeomInt_IntSS::Perform(S1, S2, Tol, Approx, ApproxS1, ApproxS2)
    /// (cxx L24-41).
    pub fn perform_surfaces(
        &mut self,
        s1: &Surface3Like,
        s2: &Surface3Like,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        self.my_hs1 = Some(GeomSurfaceAdapter::new(s1.clone()));
        // OCCT L32-39: if (S1 == S2) myHS2 = myHS1; else myHS2 = new
        // GeomAdaptor_Surface(S2).  The OCCT test is a handle identity test;
        // rcad compares the carried Geom_Surface values (see
        // `surface3_identity`).
        if surface3_identity(s1, s2) {
            self.my_hs2 = self.my_hs1.clone();
        } else {
            self.my_hs2 = Some(GeomSurfaceAdapter::new(s2.clone()));
        }
        self.internal_perform(tol, approx, approx_s1, approx_s2, false, 0., 0., 0., 0.);
    }

    /// OCCT GeomInt_IntSS::Perform(S1, S2, Tol, U1, V1, U2, V2, Approx,
    /// ApproxS1, ApproxS2) (cxx L47-68) — general intersection with a starting
    /// point.
    pub fn perform_surfaces_with_start(
        &mut self,
        s1: &Surface3Like,
        s2: &Surface3Like,
        tol: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        self.my_hs1 = Some(GeomSurfaceAdapter::new(s1.clone()));
        if surface3_identity(s1, s2) {
            self.my_hs2 = self.my_hs1.clone();
        } else {
            self.my_hs2 = Some(GeomSurfaceAdapter::new(s2.clone()));
        }
        self.internal_perform(tol, approx, approx_s1, approx_s2, true, u1, v1, u2, v2);
    }

    /// OCCT GeomInt_IntSS::Perform(HS1, HS2, Tol, Approx, ApproxS1, ApproxS2)
    /// (lxx L52-62) — intersection of adapted surfaces.
    pub fn perform_adaptors(
        &mut self,
        hs1: &GeomSurfaceAdapter,
        hs2: &GeomSurfaceAdapter,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        self.my_hs1 = Some(hs1.clone());
        self.my_hs2 = Some(hs2.clone());
        self.internal_perform(tol, approx, approx_s1, approx_s2, false, 0., 0., 0., 0.);
    }

    /// OCCT GeomInt_IntSS::Perform(HS1, HS2, Tol, U1, V1, U2, V2, Approx,
    /// ApproxS1, ApproxS2) (lxx L68-82) — intersection of adapted surfaces with
    /// a starting point.
    pub fn perform_adaptors_with_start(
        &mut self,
        hs1: &GeomSurfaceAdapter,
        hs2: &GeomSurfaceAdapter,
        tol: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        self.my_hs1 = Some(hs1.clone());
        self.my_hs2 = Some(hs2.clone());
        self.internal_perform(tol, approx, approx_s1, approx_s2, true, u1, v1, u2, v2);
    }

    // ========================================================================
    // OCCT GeomInt_IntSS::InternalPerform (cxx L72-161)
    // ========================================================================

    /// OCCT GeomInt_IntSS::InternalPerform (cxx L72-161).
    #[allow(clippy::too_many_arguments)]
    fn internal_perform(
        &mut self,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
        use_start: bool,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        // OCCT L82-84.
        self.my_tol_reached_2d = 0.0;
        self.my_tol_reached_3d = 0.0;
        self.my_nbrestr = 0;
        self.sline.clear();
        // OCCT L84-85: slineS1/slineS2 are NOT cleared by InternalPerform in
        // OCCT either (only sline is; slineS1/S2 grow with sline).
        self.sline_s1.clear();
        self.sline_s2.clear();

        // OCCT L86-92.
        let tol_arc = tol;
        let tol_tang = tol;
        let mut deflection = 0.1;
        let my_hs1 = self.my_hs1.as_ref().expect("GeomInt_IntSS::myHS1");
        let my_hs2 = self.my_hs2.as_ref().expect("GeomInt_IntSS::myHS2");
        if my_hs1.get_type() == GeomAbsSurfaceType::BSplineSurface
            && my_hs2.get_type() == GeomAbsSurfaceType::BSplineSurface
        {
            deflection /= 10.;
        }

        // OCCT L94-95.
        let dom1 = TopolTool::<GeomSurfaceAdapter, GeomTool>::new(my_hs1);
        let dom2 = TopolTool::<GeomSurfaceAdapter, GeomTool>::new(my_hs2);

        // OCCT L97-98: NCollection_DynamicArray<handle<Adaptor3d_Surface>>.
        let mut a_vec_hs1: Vec<GeomSurfaceAdapter> = Vec::new();
        let mut a_vec_hs2: Vec<GeomSurfaceAdapter> = Vec::new();

        // OCCT L100-108.
        if adapter_identity(my_hs1, my_hs2) {
            a_vec_hs1.push(my_hs1.clone());
            a_vec_hs2.push(my_hs2.clone());
        } else {
            int_patch_intersection_prepare_surfaces(
                my_hs1,
                &dom1,
                my_hs2,
                &dom2,
                tol,
                &mut a_vec_hs1,
                &mut a_vec_hs2,
            );
        }

        // OCCT L110-159.
        let mut a_num_of_hs1 = 0usize;
        while a_num_of_hs1 < a_vec_hs1.len() {
            let a_hs1 = a_vec_hs1[a_num_of_hs1].clone();
            let mut a_num_of_hs2 = 0usize;
            while a_num_of_hs2 < a_vec_hs2.len() {
                let a_hs2 = a_vec_hs2[a_num_of_hs2].clone();

                // OCCT L118-119.
                let mut a_dom1 = TopolTool::<GeomSurfaceAdapter, GeomTool>::new(&a_hs1);
                let mut a_dom2 = TopolTool::<GeomSurfaceAdapter, GeomTool>::new(&a_hs2);

                // OCCT L121-124: myLConstruct.Load(aDom1, aDom2, aHS1, aHS2).
                // rcad stores the adaptors the tools are built over (see the
                // file-header note in geom_int_line_constructor.rs).
                self.my_l_construct.load(&a_hs1, &a_hs2, &a_hs1, &a_hs2);

                // OCCT L126: UVMaxStep = DefineUVMaxStep(aHS1, aDom1, aHS2, aDom2).
                let uv_max_step =
                    define_uv_max_step(&a_hs1, &mut a_dom1, &a_hs2, &mut a_dom2);

                // OCCT L128.
                self.my_intersector
                    .set_tolerances(tol_arc, tol_tang, uv_max_step, deflection);

                // OCCT L130-148.
                if adapter_identity(&a_hs1, &a_hs2) {
                    // OCCT L132: myIntersector.Perform(aHS1, aDom1, TolArc,
                    // TolTang) — the SELF-intersection overload.
                    // rcad: IntPatch_Intersection's self-intersection overload
                    // (IntPatch_Intersection.cxx L177-...) is not translated;
                    // a fresh intersector reproduces its reset (done = false,
                    // spnt/slin cleared), which is OCCT's failure path here —
                    // the IsDone() gate below then skips the MakeCurve loop.
                    self.my_intersector = IntPatchIntersection::new();
                } else if !use_start {
                    // OCCT L136.
                    self.my_intersector.perform(
                        a_hs1.surface3(),
                        a_hs2.surface3(),
                        uv_rect(&a_hs1),
                        uv_rect(&a_hs2),
                        tol_arc,
                        tol_tang,
                    );
                } else {
                    // OCCT L140-147.
                    let a_state1 = a_dom1.classify(DVec2::new(u1, v1), tol, true);
                    let a_state2 = a_dom2.classify(DVec2::new(u2, v2), tol, true);

                    if (a_state1 == State::In || a_state1 == State::On)
                        && (a_state2 == State::In || a_state2 == State::On)
                    {
                        // OCCT L146: myIntersector.Perform(aHS1, aDom1, aHS2,
                        // aDom2, U1, V1, U2, V2, TolArc, TolTang).
                        self.my_intersector.perform_with_start(
                            a_hs1.surface3(),
                            a_hs2.surface3(),
                            uv_rect(&a_hs1),
                            uv_rect(&a_hs2),
                            u1,
                            v1,
                            u2,
                            v2,
                            tol_arc,
                            tol_tang,
                        );
                    }
                }

                // OCCT L151-158.
                if self.my_intersector.is_done() {
                    let nblin = self.my_intersector.nb_lines();
                    let mut i = 1usize;
                    while i <= nblin {
                        self.make_curve(i, &a_dom1, &a_dom2, tol, approx, approx_s1, approx_s2);
                        i += 1;
                    }
                }
                a_num_of_hs2 += 1;
            }
            a_num_of_hs1 += 1;
        }
    }

    // ========================================================================
    // OCCT GeomInt_IntSS.cxx L165-227: accessors
    // ========================================================================

    /// OCCT GeomInt_IntSS::IsDone (lxx L86-89).
    pub fn is_done(&self) -> bool {
        self.my_intersector.is_done()
    }

    /// OCCT GeomInt_IntSS::TolReached2d (lxx L93-96).
    pub fn tol_reached_2d(&self) -> f64 {
        self.my_tol_reached_2d
    }

    /// OCCT GeomInt_IntSS::TolReached3d (lxx L100-103).
    pub fn tol_reached_3d(&self) -> f64 {
        self.my_tol_reached_3d
    }

    /// OCCT GeomInt_IntSS::NbLines (lxx L107-110):
    /// sline.Length() - myNbrestr.
    pub fn nb_lines(&self) -> usize {
        self.sline.len() - self.my_nbrestr
    }

    /// OCCT GeomInt_IntSS::Line(Index) (cxx L165-169) — 1-based Index, the
    /// returned curve is the OCCT handle (a null handle is `None`).
    pub fn line(&self, index: usize) -> Option<Curve3> {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::Line");
        }
        self.sline[index + self.my_nbrestr - 1].clone()
    }

    /// OCCT GeomInt_IntSS::HasLineOnS1(index) (cxx L199-203):
    /// !slineS1(index).IsNull().
    pub fn has_line_on_s1(&self, index: usize) -> bool {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::HasLineOnS1");
        }
        self.sline_s1[index - 1].is_some()
    }

    /// OCCT GeomInt_IntSS::LineOnS1(Index) (cxx L215-219).
    pub fn line_on_s1(&self, index: usize) -> Option<Curve2d> {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::LineOnS1");
        }
        self.sline_s1[index - 1].clone()
    }

    /// OCCT GeomInt_IntSS::HasLineOnS2(index) (cxx L207-211).
    pub fn has_line_on_s2(&self, index: usize) -> bool {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::HasLineOnS2");
        }
        self.sline_s2[index - 1].is_some()
    }

    /// OCCT GeomInt_IntSS::LineOnS2(Index) (cxx L223-227).
    pub fn line_on_s2(&self, index: usize) -> Option<Curve2d> {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::LineOnS2");
        }
        self.sline_s2[index - 1].clone()
    }

    /// OCCT GeomInt_IntSS::NbBoundaries (lxx L114-118).
    pub fn nb_boundaries(&self) -> usize {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::NbBoundaries() - no result");
        }
        self.my_nbrestr
    }

    /// OCCT GeomInt_IntSS::Boundary(Index) (cxx L173-178).
    pub fn boundary(&self, index: usize) -> Option<Curve3> {
        if !self.my_intersector.is_done() {
            panic!("StdFail_NotDone: GeomInt_IntSS::Line");
        }
        if index == 0 || index > self.my_nbrestr {
            panic!("Standard_OutOfRange: GeomInt_IntSS::Boundary");
        }
        self.sline[index - 1].clone()
    }

    /// OCCT GeomInt_IntSS::NbPoints (lxx L122-125).
    pub fn nb_points(&self) -> usize {
        self.my_intersector.nb_points()
    }

    /// OCCT GeomInt_IntSS::Point(Index) (lxx L129-132):
    /// myIntersector.Point(Index).Value().
    pub fn point(&self, index: usize) -> glam::DVec3 {
        // rcad note: IntPatchIntersection::point is 0-based while OCCT's
        // IntPatch_Intersection::Point(Index) is 1-based.
        self.my_intersector.point(index - 1).p1
    }

    /// OCCT GeomInt_IntSS::Pnt2d(Index, OnFirst) (cxx L182-195).
    pub fn pnt2d(&self, index: usize, on_first: bool) -> DVec2 {
        let thept = self.my_intersector.point(index - 1);
        if on_first {
            // thept.ParametersOnS1(U, V)
            DVec2::new(thept.u1, thept.v1)
        } else {
            DVec2::new(thept.u2, thept.v2)
        }
    }

    /// OCCT GeomInt_IntSS::SetTolFixTangents (hxx L121).
    pub fn set_tol_fix_tangents(&mut self, a_tol_check: f64, a_tol_ang_check: f64) {
        self.my_tol_check = a_tol_check;
        self.my_tol_ang_check = a_tol_ang_check;
    }

    /// OCCT GeomInt_IntSS::TolFixTangents (hxx L123).
    pub fn tol_fix_tangents(&self) -> (f64, f64) {
        (self.my_tol_check, self.my_tol_ang_check)
    }

    // ---- crate-internal accessors used by geom_int_int_ss_1.rs ------------

    pub(crate) fn my_intersector(&self) -> &IntPatchIntersection {
        &self.my_intersector
    }
    pub(crate) fn my_l_construct(&self) -> &GeomIntLineConstructor {
        &self.my_l_construct
    }
    pub(crate) fn my_l_construct_mut(&mut self) -> &mut GeomIntLineConstructor {
        &mut self.my_l_construct
    }
    pub(crate) fn my_hs1(&self) -> &GeomSurfaceAdapter {
        self.my_hs1.as_ref().expect("GeomInt_IntSS::myHS1")
    }
    pub(crate) fn my_hs2(&self) -> &GeomSurfaceAdapter {
        self.my_hs2.as_ref().expect("GeomInt_IntSS::myHS2")
    }
    pub(crate) fn my_tol_reached_2d_mut(&mut self) -> &mut f64 {
        &mut self.my_tol_reached_2d
    }
    pub(crate) fn my_tol_reached_3d_mut(&mut self) -> &mut f64 {
        &mut self.my_tol_reached_3d
    }
    pub(crate) fn my_tol_check(&self) -> f64 {
        self.my_tol_check
    }
    pub(crate) fn my_tol_ang_check(&self) -> f64 {
        self.my_tol_ang_check
    }
    pub(crate) fn sline_mut(&mut self) -> &mut Vec<Option<Curve3>> {
        &mut self.sline
    }
    pub(crate) fn sline_s1_mut(&mut self) -> &mut Vec<Option<Curve2d>> {
        &mut self.sline_s1
    }
    pub(crate) fn sline_s2_mut(&mut self) -> &mut Vec<Option<Curve2d>> {
        &mut self.sline_s2
    }
}

/// The rcad carrier of `occ::handle<Geom_Surface>` at the GeomInt_IntSS API
/// boundary: the kernel `Geom_Surface` value.
pub type Surface3Like = rcad_kernel::geom::Surface3;

/// OCCT the adaptor's UV domain used as the `Adaptor3d_TopolTool` domain
/// rectangle (rcad's IntPatchIntersection takes the corrected UV rectangles in
/// place of the tool handles — the same encoding the FF pipeline uses).
pub(crate) fn uv_rect(hs: &GeomSurfaceAdapter) -> [f64; 4] {
    use crate::hlr::contap::surface_adaptor::SurfaceAdapter;
    [
        hs.first_u_parameter(),
        hs.last_u_parameter(),
        hs.first_v_parameter(),
        hs.last_v_parameter(),
    ]
}

/// `GeomInt_LineConstructor::Load(aDom1, aDom2, aHS1, aHS2)` needs the tools;
/// rcad stores the adaptors they are built over (see the file-header note in
/// `geom_int_line_constructor.rs`).

use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use rcad_kernel::topods::State;
