// OCCT IntWalk_IWalking.gxx / IntPatch_TheIWalking 1:1 Rust translation —
// the implicit-surface walking algorithm.
//
// Walks along the intersection curve F(u,v)=0 on the parametric surface,
// starting from boundary points (IntSurf_PathPoint) and interior points
// (IntSurf_InteriorPoint), producing polylines (IntPatch_TheIWLineOfTheIWalking).
//
// OCCT IntWalk_IWalking.gxx L1-3152, IntWalk_IWalking.lxx L18-60.

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};

use crate::geomalgo::int_surf::{LineOn2S, PntOn2S};

use super::function_set_root::FunctionSetRoot;
use super::path_point::{InteriorPoint, PathPoint};
use super::surf_function::SurfFunction;

// OCCT constants (IntWalk_IWalking.gxx L36-40).
const COS_REF_3D: f64 = 0.98; // correspond to 11.478 deg
const COS_REF_2D: f64 = 0.88; // correspond to 25 deg
const MAX_DIVISION: i32 = 60; // max number of step division

/// OCCT IntWalk_StatusDeflection.hxx.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatusDeflection {
    PasTropGrand,
    StepTooSmall,
    PointConfondu,
    ArretSurPointPrecedent,
    ArretSurPoint,
    OK,
}

impl StatusDeflection {
    fn index(&self) -> i32 {
        match self {
            StatusDeflection::PasTropGrand => 0,
            StatusDeflection::StepTooSmall => 1,
            StatusDeflection::PointConfondu => 2,
            StatusDeflection::ArretSurPointPrecedent => 3,
            StatusDeflection::ArretSurPoint => 4,
            StatusDeflection::OK => 5,
        }
    }
}

/// OCCT IntWalk_WalkingData — state + UV of a start point (wd1/wd2).
#[derive(Clone, Debug)]
struct WalkingData {
    etat: i32,
    ustart: f64,
    vstart: f64,
}

impl WalkingData {
    fn dummy() -> Self {
        WalkingData {
            etat: -10,
            ustart: 0.0,
            vstart: 0.0,
        }
    }
}

/// OCCT Bnd_Range — an interval used for the estimated U/V range of the
/// section curve (mySRangeU / mySRangeV).
#[derive(Clone, Debug)]
struct BndRange {
    is_void: bool,
    a: f64,
    b: f64,
}

impl BndRange {
    fn new() -> Self {
        BndRange {
            is_void: true,
            a: 0.0,
            b: 0.0,
        }
    }
    fn from_bounds(a: f64, b: f64) -> Self {
        BndRange {
            is_void: false,
            a,
            b,
        }
    }
    fn add(&mut self, x: f64) {
        if self.is_void {
            self.is_void = false;
            self.a = x;
            self.b = x;
        } else {
            self.a = self.a.min(x);
            self.b = self.b.max(x);
        }
    }
    fn delta(&self) -> f64 {
        if self.is_void {
            0.0
        } else {
            self.b - self.a
        }
    }
    fn is_void(&self) -> bool {
        self.is_void
    }
    fn enlarge(&mut self, d: f64) {
        if !self.is_void {
            self.a -= d;
            self.b += d;
        }
    }
    fn common(&mut self, other: &BndRange) {
        if self.is_void || other.is_void {
            self.is_void = true;
            return;
        }
        self.a = self.a.max(other.a);
        self.b = self.b.min(other.b);
        if self.a > self.b {
            self.is_void = true;
        }
    }
}

/// OCCT IntSurf_Couple — (line point index, passing point index).
#[derive(Debug, Clone, Copy)]
struct Couple {
    first: usize,
    second: i32,
}

/// OCCT IntPatch_TheIWLineOfTheIWalking — a polyline produced by the walking.
#[derive(Clone, Debug)]
pub struct IWLine {
    line: LineOn2S,
    couple: Vec<Couple>,
    closed: bool,
    has_first: bool,
    has_last: bool,
    first_index: i32,
    last_index: i32,
    the_first_point: PathPoint,
    the_last_point: PathPoint,
    indextg: i32,
    vcttg: DVec3,
    istgtbeg: bool,
    istgtend: bool,
}

impl IWLine {
    /// OCCT constructor (IntPatch_TheIWLineOfTheIWalking_0.cxx L25-40).
    pub fn new() -> Self {
        IWLine {
            line: LineOn2S::new(),
            couple: Vec::new(),
            closed: false,
            has_first: false,
            has_last: false,
            first_index: -1,
            last_index: -1,
            the_first_point: PathPoint::new(),
            the_last_point: PathPoint::new(),
            indextg: -1,
            vcttg: DVec3::ZERO,
            istgtbeg: false,
            istgtend: false,
        }
    }

    /// OCCT Reverse() — reverse the points; HasFirst/HasLast are kept.
    pub fn reverse(&mut self) {
        self.line.reverse();
        let n = self.line.nb_points();
        let nb_couple = self.couple.len();
        for i in 0..nb_couple {
            let c = self.couple[i];
            self.couple[i] = Couple {
                first: n - c.first + 1,
                second: c.second,
            };
        }
    }

    /// OCCT Cut(Index) — split the line at the 1-based point of rank Index.
    pub fn cut(&mut self, index: usize) {
        let _lost = self.line.split(index);
    }

    /// OCCT AddPoint(P).
    pub fn add_point(&mut self, p: &PntOn2S) {
        self.line.add(p);
    }

    /// OCCT AddStatusFirst(Closed, HasFirst).
    pub fn add_status_first(&mut self, closed: bool, has_first: bool) {
        self.closed = closed;
        self.has_first = has_first;
    }

    /// OCCT AddStatusFirst(Closed, HasFirst, Index, P).
    pub fn add_status_first_full(&mut self, closed: bool, has_first: bool, index: i32, p: &PathPoint) {
        self.closed = closed;
        self.has_first = has_first;
        self.first_index = index;
        self.the_first_point = p.clone();
    }

    /// OCCT AddStatusFirstLast(Closed, HasFirst, HasLast).
    pub fn add_status_first_last(&mut self, closed: bool, has_first: bool, has_last: bool) {
        self.closed = closed;
        self.has_first = has_first;
        self.has_last = has_last;
    }

    /// OCCT AddStatusLast(HasLast).
    pub fn add_status_last(&mut self, has_last: bool) {
        self.has_last = has_last;
    }

    /// OCCT AddStatusLast(HasLast, Index, P).
    pub fn add_status_last_full(&mut self, has_last: bool, index: i32, p: &PathPoint) {
        self.has_last = has_last;
        self.last_index = index;
        self.the_last_point = p.clone();
    }

    /// OCCT AddIndexPassing(Index).
    pub fn add_index_passing(&mut self, index: i32) {
        self.couple.push(Couple {
            first: self.line.nb_points() + 1,
            second: index,
        });
    }

    /// OCCT SetTangentVector(V, Index).
    pub fn set_tangent_vector(&mut self, v: DVec3, index: i32) {
        self.indextg = index;
        self.vcttg = v;
    }

    /// OCCT SetTangencyAtBegining.
    pub fn set_tangency_at_begining(&mut self, is_tangent: bool) {
        self.istgtbeg = is_tangent;
    }

    /// OCCT SetTangencyAtEnd.
    pub fn set_tangency_at_end(&mut self, is_tangent: bool) {
        self.istgtend = is_tangent;
    }

    /// OCCT NbPoints() — 1-based count.
    pub fn nb_points(&self) -> usize {
        self.line.nb_points()
    }

    /// OCCT Value(Index) — 1-based.
    pub fn value(&self, index: usize) -> &PntOn2S {
        self.line.value(index - 1)
    }

    /// OCCT Line().
    pub fn line(&self) -> &LineOn2S {
        &self.line
    }

    /// OCCT IsClosed().
    pub fn is_closed(&self) -> bool {
        self.closed
    }
    /// OCCT HasFirstPoint().
    pub fn has_first_point(&self) -> bool {
        self.has_first
    }
    /// OCCT HasLastPoint().
    pub fn has_last_point(&self) -> bool {
        self.has_last
    }
    /// OCCT FirstPoint().
    pub fn first_point(&self) -> &PathPoint {
        &self.the_first_point
    }
    /// OCCT FirstPointIndex().
    pub fn first_point_index(&self) -> i32 {
        self.first_index
    }
    /// OCCT LastPoint().
    pub fn last_point(&self) -> &PathPoint {
        &self.the_last_point
    }
    /// OCCT LastPointIndex().
    pub fn last_point_index(&self) -> i32 {
        self.last_index
    }
    /// OCCT TangentVector(Index) — sets Index = indextg, returns vcttg.
    pub fn tangent_vector(&self) -> (DVec3, i32) {
        (self.vcttg, self.indextg)
    }
    /// OCCT IsTangentAtBegining().
    pub fn is_tangent_at_begining(&self) -> bool {
        self.istgtbeg
    }
    /// OCCT IsTangentAtEnd().
    pub fn is_tangent_at_end(&self) -> bool {
        self.istgtend
    }
    /// OCCT NbPassingPoint().
    pub fn nb_passing_point(&self) -> usize {
        self.couple.len()
    }
    /// OCCT PassingPoint(Index, IndexLine, IndexPnts).
    pub fn passing_point(&self, index: usize, index_line: &mut usize, index_pnts: &mut i32) {
        let c = self.couple[index - 1];
        *index_line = c.first;
        *index_pnts = c.second;
    }
}

impl Default for IWLine {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT IntPatch_TheIWalking.
pub struct IWalking {
    done: bool,
    seq_single: Vec<PathPoint>,
    fleche: f64,
    pas: f64,
    tolerance: [f64; 2],
    epsilon: f64,
    reversed: bool,
    wd1: Vec<WalkingData>,
    wd2: Vec<WalkingData>,
    nb_multiplicities: Vec<i32>,
    my_s_range_u: BndRange,
    my_s_range_v: BndRange,
    um: f64,
    um_max: f64,
    vm: f64,
    vm_max: f64,
    previous_point: PntOn2S,
    previous_d3d: DVec3,
    previous_d2d: DVec2,
    seq_ajout: Vec<i32>,
    seq_alone: Vec<i32>,
    point_line_line: HashMap<i32, Vec<i32>>,
    lines: Vec<IWLine>,
    to_fill_holes: bool,
}

impl IWalking {
    /// OCCT IntWalk_IWalking(Epsilon, Deflection, Step, ToFillHoles = false)
    /// (gxx L90-109).
    pub fn new(epsilon: f64, deflection: f64, step: f64, to_fill_holes: bool) -> Self {
        IWalking {
            done: false,
            seq_single: Vec::new(),
            fleche: deflection,
            pas: step,
            tolerance: [0.0; 2],
            epsilon: epsilon * epsilon,
            reversed: false,
            wd1: Vec::new(),
            wd2: Vec::new(),
            nb_multiplicities: Vec::new(),
            my_s_range_u: BndRange::new(),
            my_s_range_v: BndRange::new(),
            um: 0.0,
            um_max: 0.0,
            vm: 0.0,
            vm_max: 0.0,
            previous_point: PntOn2S::new(),
            previous_d3d: DVec3::ZERO,
            previous_d2d: DVec2::ZERO,
            seq_ajout: Vec::new(),
            seq_alone: Vec::new(),
            point_line_line: HashMap::new(),
            lines: Vec::new(),
            to_fill_holes,
        }
    }

    /// OCCT SetTolerance (lxx L21-29).
    pub fn set_tolerance(&mut self, epsilon: f64, deflection: f64, step: f64) {
        self.fleche = deflection;
        self.pas = step;
        self.epsilon = epsilon * epsilon;
    }

    /// OCCT Clear (gxx L119-134).
    fn clear(&mut self) {
        self.wd1.clear();
        self.wd2.clear();
        self.wd1.push(WalkingData::dummy());
        self.wd2.push(WalkingData::dummy());
        self.nb_multiplicities.clear();
        self.nb_multiplicities.push(-1);
        self.done = false;
        self.seq_ajout.clear();
        self.lines.clear();
    }

    /// OCCT IsDone (lxx L33-35).
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT NbLines (lxx L37-41).
    pub fn nb_lines(&self) -> usize {
        self.lines.len()
    }
    /// OCCT Value(Index) — 1-based.
    pub fn value(&self, index: usize) -> &IWLine {
        &self.lines[index - 1]
    }
    /// OCCT NbSinglePnts (lxx L43-47).
    pub fn nb_single_points(&self) -> usize {
        self.seq_single.len()
    }
    /// OCCT SinglePnt(Index) — 1-based.
    pub fn single_pnt(&self, index: usize) -> &PathPoint {
        &self.seq_single[index - 1]
    }

    /// OCCT IsTangentExtCheck (gxx L52-88).
    fn is_tangent_ext_check(
        func: &mut SurfFunction,
        u: f64,
        v: f64,
        step_u: f64,
        step_v: f64,
        u_inf: f64,
        u_sup: f64,
        v_inf: f64,
        v_sup: f64,
    ) -> bool {
        let a_tol = func.tolerance();
        let a_par_u = [
            (u + step_u).min(u_sup),
            (u - step_u).max(u_inf),
            u,
            u,
        ];
        let a_par_v = [
            v,
            v,
            (v + step_v).min(v_sup),
            (v - step_v).max(v_inf),
        ];
        for i in 0..4 {
            let x = [a_par_u[i], a_par_v[i]];
            if func.value(&x).is_some() {
                if func.root().abs() > a_tol {
                    return false;
                }
            }
        }
        true
    }

    /// OCCT Perform(Pnts1, Pnts2, Func, Caro, Reversed) (gxx L153-311).
    /// `domain` is the corrected face UV rectangle ([u_min, u_max, v_min,
    /// v_max]) — the OCCT adaptor surface carries the restricted face domain,
    /// while the rcad Surface3 exposes the natural (possibly infinite) domain.
    pub fn perform(
        &mut self,
        pnts1: &[PathPoint],
        pnts2: &[InteriorPoint],
        func: &mut SurfFunction,
        caro: &Surface3,
        domain: [f64; 4],
        reversed: bool,
    ) {
        let nb_pnts1 = pnts1.len();
        let nb_pnts2 = pnts2.len();
        let mut rajout = false;

        self.clear();
        self.reversed = reversed;

        self.um = domain[0];
        self.vm = domain[2];
        self.um_max = domain[1];
        self.vm_max = domain[3];

        if self.um_max < self.um {
            let utemp = self.um_max;
            self.um_max = self.um;
            self.um = utemp;
        }
        if self.vm_max < self.vm {
            let vtemp = self.vm_max;
            self.vm_max = self.vm;
            self.vm = vtemp;
        }

        let a_step_u = self.pas * (self.um_max - self.um);
        let a_step_v = self.pas * (self.vm_max - self.vm);

        let mut u_mult: Vec<f64> = Vec::new();
        let mut v_mult: Vec<f64> = Vec::new();

        let decal = 0;
        for i in 1..=(nb_pnts1 + decal) {
            let path_pnt = &pnts1[i - 1 - decal];
            let mut a_wd1 = WalkingData {
                etat: 1,
                ustart: 0.0,
                vstart: 0.0,
            };
            if !path_pnt.is_passing_pnt() {
                a_wd1.etat = 11;
            }
            if !path_pnt.is_tangent() {
                a_wd1.etat += 1;
            }
            if a_wd1.etat == 2 {
                a_wd1.etat = 11;
            }
            let uv2 = path_pnt.value_2d();
            a_wd1.ustart = uv2.x;
            a_wd1.vstart = uv2.y;
            self.my_s_range_u.add(a_wd1.ustart);
            self.my_s_range_v.add(a_wd1.vstart);

            self.wd1.push(a_wd1);
            let a_nb_mult = path_pnt.multiplicity();
            self.nb_multiplicities.push(a_nb_mult);

            for j in 1..=a_nb_mult {
                let mut u = 0.0;
                let mut v = 0.0;
                path_pnt.parameters(j, &mut u, &mut v);
                u_mult.push(u);
                v_mult.push(v);
            }
        }

        for i in 1..=nb_pnts2 {
            let an_ip = &pnts2[i - 1];
            let mut a_wd2 = WalkingData {
                etat: 1,
                ustart: 0.0,
                vstart: 0.0,
            };
            a_wd2.ustart = an_ip.u_parameter();
            a_wd2.vstart = an_ip.v_parameter();
            self.my_s_range_u.add(a_wd2.ustart);
            self.my_s_range_v.add(a_wd2.vstart);

            if !Self::is_tangent_ext_check(
                func,
                a_wd2.ustart,
                a_wd2.vstart,
                a_step_u,
                a_step_v,
                self.um,
                self.um_max,
                self.vm,
                self.vm_max,
            ) {
                a_wd2.etat = 13;
            }
            self.wd2.push(a_wd2);
        }

        self.tolerance = [
            u_resolution(domain, rcad_kernel::precision::CONFUSION),
            v_resolution(domain, rcad_kernel::precision::CONFUSION),
        ];

        func.set_surface(caro.clone());

        if self.my_s_range_u.delta() > self.tolerance[0].max(rcad_kernel::precision::PCONFUSION) {
            self.my_s_range_u.enlarge(self.my_s_range_u.delta());
            let full = BndRange::from_bounds(self.um, self.um_max);
            self.my_s_range_u.common(&full);
        } else {
            self.my_s_range_u = BndRange::from_bounds(self.um, self.um_max);
        }

        if self.my_s_range_v.delta() > self.tolerance[1].max(rcad_kernel::precision::PCONFUSION) {
            self.my_s_range_v.enlarge(self.my_s_range_v.delta());
            let full = BndRange::from_bounds(self.vm, self.vm_max);
            self.my_s_range_v.common(&full);
        } else {
            self.my_s_range_v = BndRange::from_bounds(self.vm, self.vm_max);
        }

        // Calculation of all open lines.
        if nb_pnts1 != 0 {
            self.compute_open_line(&u_mult, &v_mult, pnts1, func, &mut rajout);
        }

        // Calculation of all closed lines.
        if nb_pnts2 != 0 {
            self.compute_close_line(&u_mult, &v_mult, pnts1, pnts2, func, &mut rajout);
        }

        if self.to_fill_holes {
            let max_nb_iter = 10;
            let mut nb_iter = 0;
            while self.seq_alone.len() > 1 && nb_iter < max_nb_iter {
                nb_iter += 1;
                let copy_seq_alone = self.seq_alone.clone();
                let mut pnts_in_holes: Vec<InteriorPoint> = Vec::new();
                self.fill_pnts_in_holes(func, copy_seq_alone, &mut pnts_in_holes);
                self.wd2.clear();
                self.wd2.push(WalkingData::dummy());
                let nb_holes = pnts_in_holes.len();
                for i in 1..=nb_holes {
                    let an_ip = &pnts_in_holes[i - 1];
                    let mut a_wd2 = WalkingData {
                        etat: 13,
                        ustart: 0.0,
                        vstart: 0.0,
                    };
                    a_wd2.ustart = an_ip.u_parameter();
                    a_wd2.vstart = an_ip.v_parameter();
                    self.wd2.push(a_wd2);
                }
                self.compute_close_line(&u_mult, &v_mult, pnts1, &pnts_in_holes, func, &mut rajout);
            }
        }

        for i in 1..=nb_pnts1 {
            if self.wd1[i].etat > 0 {
                self.seq_single.push(pnts1[i - 1].clone());
            }
        }
        self.done = true;
    }

    /// OCCT Perform(Pnts1, Func, Caro, Reversed) (gxx L317-391) — without
    /// interior points.
    #[allow(dead_code)]
    pub fn perform_no_interior(
        &mut self,
        pnts1: &[PathPoint],
        func: &mut SurfFunction,
        caro: &Surface3,
        domain: [f64; 4],
        reversed: bool,
    ) {
        let mut rajout = false;
        let nb_pnts1 = pnts1.len();

        self.reversed = reversed;

        let mut u_mult: Vec<f64> = Vec::new();
        let mut v_mult: Vec<f64> = Vec::new();

        for i in 1..=nb_pnts1 {
            let path_pnt = &pnts1[i - 1];
            let mut a_wd1 = WalkingData {
                etat: 1,
                ustart: 0.0,
                vstart: 0.0,
            };
            if !path_pnt.is_passing_pnt() {
                a_wd1.etat = 11;
            }
            if !path_pnt.is_tangent() {
                a_wd1.etat += 1;
            }
            let uv2 = path_pnt.value_2d();
            a_wd1.ustart = uv2.x;
            a_wd1.vstart = uv2.y;
            self.wd1.push(a_wd1);
            let a_nb_mult = path_pnt.multiplicity();
            self.nb_multiplicities.push(a_nb_mult);

            for j in 1..=a_nb_mult {
                let mut u = 0.0;
                let mut v = 0.0;
                path_pnt.parameters(j, &mut u, &mut v);
                u_mult.push(u);
                v_mult.push(v);
            }
        }

        self.tolerance = [
            u_resolution(domain, rcad_kernel::precision::CONFUSION),
            v_resolution(domain, rcad_kernel::precision::CONFUSION),
        ];

        self.um = domain[0];
        self.vm = domain[2];
        self.um_max = domain[1];
        self.vm_max = domain[3];

        if self.um_max < self.um {
            let utemp = self.um_max;
            self.um_max = self.um;
            self.um = utemp;
        }
        if self.vm_max < self.vm {
            let vtemp = self.vm_max;
            self.vm_max = self.vm;
            self.vm = vtemp;
        }

        func.set_surface(caro.clone());

        if nb_pnts1 != 0 {
            self.compute_open_line(&u_mult, &v_mult, pnts1, func, &mut rajout);
        }

        for i in 1..=nb_pnts1 {
            if self.wd1[i].etat > 0 {
                self.seq_single.push(pnts1[i - 1].clone());
            }
        }
        self.done = true;
    }

    // =====================================================================
    // Cadrage (gxx L413-591)
    // =====================================================================
    fn cadrage(
        &mut self,
        born_inf: &mut [f64; 2],
        born_sup: &mut [f64; 2],
        uvap: &mut [f64; 2],
        step: &mut f64,
        step_sign: i32,
    ) -> bool {
        let duvx = self.previous_d2d.x;
        let duvy = self.previous_d2d.y;

        if !self.reversed {
            let (u, v) = self.previous_point.parameters_on_surface(false);
            uvap[0] = u;
            uvap[1] = v;
        } else {
            let (u, v) = self.previous_point.parameters_on_surface(true);
            uvap[0] = u;
            uvap[1] = v;
        }

        let u1 = uvap[0] + *step * duvx * step_sign as f64;
        let v1 = uvap[1] + *step * duvy * step_sign as f64;

        let infu = u1 <= born_inf[0] + rcad_kernel::precision::PCONFUSION;
        let supu = u1 >= born_sup[0] - rcad_kernel::precision::PCONFUSION;
        let infv = v1 <= born_inf[1] + rcad_kernel::precision::PCONFUSION;
        let supv = v1 >= born_sup[1] - rcad_kernel::precision::PCONFUSION;

        let the_step_u;
        let the_step_v;

        if !infu && !supu && !infv && !supv {
            uvap[0] = u1;
            uvap[1] = v1;
            return false;
        }

        if (infu || supu) && (infv || supv) {
            if infu {
                the_step_u = if duvx != 0.0 {
                    ((born_inf[0] - uvap[0]) / duvx).abs()
                } else {
                    *step
                };
            } else {
                the_step_u = if duvx != 0.0 {
                    ((born_sup[0] - uvap[0]) / duvx).abs()
                } else {
                    *step
                };
            }
            if infv {
                the_step_v = if duvy != 0.0 {
                    ((born_inf[1] - uvap[1]) / duvy).abs()
                } else {
                    *step
                };
            } else {
                the_step_v = if duvy != 0.0 {
                    ((born_sup[1] - uvap[1]) / duvy).abs()
                } else {
                    *step
                };
            }

            if the_step_u <= the_step_v {
                *step = the_step_u;
                if infu {
                    uvap[0] = born_inf[0];
                    born_sup[0] = born_inf[0];
                } else {
                    uvap[0] = born_sup[0];
                    born_inf[0] = born_sup[0];
                }
                uvap[1] += *step * duvy * step_sign as f64;
            } else {
                *step = the_step_v;
                if infv {
                    uvap[1] = born_inf[1];
                    born_sup[1] = born_inf[1];
                } else {
                    uvap[1] = born_sup[1];
                    born_inf[1] = born_sup[1];
                }
                uvap[0] += *step * duvx * step_sign as f64;
            }
            return true;
        } else if infu {
            if duvx != 0.0 {
                let a_step = ((born_inf[0] - uvap[0]) / duvx).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_sup[0] = born_inf[0];
            uvap[0] = born_inf[0];
            uvap[1] += *step * duvy * step_sign as f64;
            return true;
        } else if supu {
            if duvx != 0.0 {
                let a_step = ((born_sup[0] - uvap[0]) / duvx).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_inf[0] = born_sup[0];
            uvap[0] = born_sup[0];
            uvap[1] += *step * duvy * step_sign as f64;
            return true;
        } else if infv {
            if duvy != 0.0 {
                let a_step = ((born_inf[1] - uvap[1]) / duvy).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_sup[1] = born_inf[1];
            uvap[0] += *step * duvx * step_sign as f64;
            uvap[1] = born_inf[1];
            return true;
        } else if supv {
            if duvy != 0.0 {
                let a_step = ((born_sup[1] - uvap[1]) / duvy).abs();
                if a_step < *step {
                    *step = a_step;
                }
            }
            born_inf[1] = born_sup[1];
            uvap[0] += *step * duvx * step_sign as f64;
            uvap[1] = born_sup[1];
            return true;
        }
        true
    }

    // =====================================================================
    // TestArretPassage — open lines (gxx L593-775)
    // =====================================================================
    fn test_arret_passage_open(
        &mut self,
        u_mult: &[f64],
        v_mult: &[f64],
        func: &mut SurfFunction,
        uv: &mut [f64; 2],
        irang: &mut i32,
    ) -> bool {
        let tolu = self.tolerance[0];
        let tolv = self.tolerance[1];
        let tolu2 = 10.0 * self.tolerance[0];
        let tolv2 = 10.0 * self.tolerance[1];

        let mut arrive = false;

        let (up, vp) = if !self.reversed {
            self.previous_point.parameters_on_surface(false)
        } else {
            self.previous_point.parameters_on_surface(true)
        };

        // Crossing test on interior points (wd2).
        for i in 1..self.wd2.len() {
            if self.wd2[i].etat > 0 {
                let utest = self.wd2[i].ustart;
                let vtest = self.wd2[i].vstart;

                let du = uv[0] - utest;
                let dv = uv[1] - vtest;
                let dup = up - utest;
                let dvp = vp - vtest;

                if (du.abs() < tolu2 && dv.abs() < tolv2)
                    || (dup.abs() < tolu2 && dvp.abs() < tolv2)
                {
                    self.wd2[i].etat = -self.wd2[i].etat;
                } else {
                    let d_du = uv[0] - up;
                    let d_dv = uv[1] - vp;
                    let d_dd = d_du * d_du + d_dv * d_dv;
                    let dd1 = du * du + dv * dv;
                    if dd1 <= d_dd {
                        let dd2 = dup * dup + dvp * dvp;
                        if dd2 <= d_dd && ((du * dup) + (dv * dvp * tolu / tolv) <= 0.0) {
                            self.wd2[i].etat = -self.wd2[i].etat;
                        }
                    }
                }
            }
        }

        // Stop test on points given at input and not yet processed.
        let mut i_candidates: Vec<usize> = Vec::new();
        let mut sq_dist_candidates: Vec<f64> = Vec::new();

        let mut l = 1;
        while l <= 2 && !arrive {
            for i in 1..self.wd1.len() {
                let is_to_check = if l == 1 {
                    self.wd1[i].etat > 0
                } else {
                    self.wd1[i].etat < 0
                };
                if is_to_check {
                    let utest = self.wd1[i].ustart;
                    let vtest = self.wd1[i].vstart;
                    let dup = up - utest;
                    let dvp = vp - vtest;
                    if dup.abs() >= tolu || dvp.abs() >= tolv {
                        let uv1m_utest = uv[0] - utest;
                        let uv2m_vtest = uv[1] - vtest;
                        if ((dup * uv1m_utest + dvp * uv2m_vtest) < 0.0)
                            || (uv1m_utest.abs() < tolu && uv2m_vtest.abs() < tolv)
                        {
                            i_candidates.push(i);
                            sq_dist_candidates.push(dup * dup + dvp * dvp);
                        } else if i < self.nb_multiplicities.len()
                            && self.nb_multiplicities[i] > 0
                            && i_candidates.is_empty()
                        {
                            let mut n: usize = 0;
                            for k in 1..i {
                                n += self.nb_multiplicities[k] as usize;
                            }
                            let mut j = n;
                            while j < n + self.nb_multiplicities[i] as usize {
                                if j < u_mult.len()
                                    && ((up - u_mult[j]) * (uv[0] - u_mult[j])
                                        + (vp - v_mult[j]) * (uv[1] - v_mult[j])
                                        < 0.0
                                        || (uv[0] - u_mult[j]).abs() < tolu
                                            && (uv[1] - v_mult[j]).abs() < tolv)
                                {
                                    *irang = i as i32;
                                    arrive = true;
                                    uv[0] = utest;
                                    uv[1] = vtest;
                                    break;
                                }
                                j += 1;
                            }
                        }
                        if arrive {
                            let _ = func.values(&[uv[0], uv[1]]);
                            break;
                        }
                    }
                }
            }
            if !i_candidates.is_empty() {
                let mut min_sq_dist = f64::MAX;
                for ind in 0..i_candidates.len() {
                    if sq_dist_candidates[ind] < min_sq_dist {
                        min_sq_dist = sq_dist_candidates[ind];
                        *irang = i_candidates[ind] as i32;
                    }
                }
                arrive = true;
                uv[0] = self.wd1[*irang as usize].ustart;
                uv[1] = self.wd1[*irang as usize].vstart;
            }
            l += 1;
        }
        arrive
    }

    // =====================================================================
    // TestArretPassage — closed lines (gxx L777-963)
    // =====================================================================
    fn test_arret_passage_close(
        &mut self,
        u_mult: &[f64],
        v_mult: &[f64],
        uv: &[f64; 2],
        index: i32,
        irang: &mut i32,
    ) -> bool {
        let mut tolu = self.tolerance[0];
        let mut tolv = self.tolerance[1];

        let (up0, vp0) = if !self.reversed {
            self.previous_point.parameters_on_surface(false)
        } else {
            self.previous_point.parameters_on_surface(true)
        };
        let mut up = up0;
        let mut vp = vp0;
        let mut uv1 = uv[0];
        let mut uv2 = uv[1];

        // Normalizing factor.
        let is_highly_anisotropic = tolu.max(tolv) > 1000.0 * tolu.min(tolv);
        let deltau = if self.my_s_range_u.is_void() {
            self.um_max - self.um
        } else if is_highly_anisotropic {
            self.my_s_range_u.delta()
        } else {
            self.my_s_range_u.delta().max(1.0)
        };
        let deltav = if self.my_s_range_v.is_void() {
            self.vm_max - self.vm
        } else if is_highly_anisotropic {
            self.my_s_range_v.delta()
        } else {
            self.my_s_range_v.delta().max(1.0)
        };

        up /= deltau;
        uv1 /= deltau;
        vp /= deltav;
        uv2 /= deltav;

        tolu /= deltau;
        tolv /= deltav;

        let tolu2 = tolu + tolu;
        let tolv2 = tolv + tolv;

        let d_previous_current = (up - uv1) * (up - uv1) + (vp - uv2) * (vp - uv2);
        let mut arrive = false;

        for k in 1..self.wd2.len() {
            if self.wd2[k].etat > 0 {
                let utest = self.wd2[k].ustart / deltau;
                let vtest = self.wd2[k].vstart / deltav;

                let uv1m_utest = uv1 - utest;
                let uv2m_vtest = uv2 - vtest;
                if (uv1m_utest < tolu2 && uv1m_utest > -tolu2)
                    && (uv2m_vtest < tolv2 && uv2m_vtest > -tolv2)
                {
                    if index != k as i32 {
                        self.wd2[k].etat = -self.wd2[k].etat;
                    } else {
                        arrive = true;
                    }
                } else {
                    let upm_utest = up - utest;
                    let vpm_vtest = vp - vtest;
                    let d_previous_start = upm_utest * upm_utest + vpm_vtest * vpm_vtest;
                    let d_current_start = uv1m_utest * uv1m_utest + uv2m_vtest * uv2m_vtest;

                    let scal = upm_utest * uv1m_utest + vpm_vtest * uv2m_vtest;
                    if upm_utest.abs() < tolu && vpm_vtest.abs() < tolv {
                        if index != k as i32 {
                            self.wd2[k].etat = -self.wd2[k].etat;
                        }
                    } else if scal < 0.0 && (d_previous_start + d_current_start < d_previous_current) {
                        if index == k as i32 {
                            arrive = true;
                        } else {
                            self.wd2[k].etat = -self.wd2[k].etat;
                        }
                    } else if k != index as usize {
                        if d_previous_start < d_previous_current * 0.25 {
                            self.wd2[k].etat = -self.wd2[k].etat;
                        } else if d_current_start < d_previous_current * 0.25 {
                            self.wd2[k].etat = -self.wd2[k].etat;
                        } else {
                            let u_mid_utest = 0.5 * (uv1 + up) - utest;
                            let v_mid_vtest = 0.5 * (uv2 + vp) - vtest;
                            let d_middle_start = u_mid_utest * u_mid_utest + v_mid_vtest * v_mid_vtest;
                            if d_middle_start < d_previous_current * 0.5 {
                                self.wd2[k].etat = -self.wd2[k].etat;
                            }
                        }
                    }
                }
            }
        }

        // Crossing test on crossing points.
        *irang = 0;
        for i in 1..self.wd1.len() {
            if self.wd1[i].etat > 0 && self.wd1[i].etat < 11 {
                let utest = self.wd1[i].ustart / deltau;
                let vtest = self.wd1[i].vstart / deltav;

                if ((up - utest) * (uv1 - utest) + (vp - vtest) * (uv2 - vtest) < 0.0)
                    || ((uv1 - utest).abs() < tolu && (uv2 - vtest).abs() < tolv)
                {
                    *irang = i as i32;
                } else if i < self.nb_multiplicities.len() && self.nb_multiplicities[i] > 0 {
                    let mut n: usize = 0;
                    for k in 1..i {
                        n += self.nb_multiplicities[k] as usize;
                    }
                    let mut j = n;
                    while j < n + self.nb_multiplicities[i] as usize {
                        if j < u_mult.len() {
                            let u_multj = u_mult[j] / deltau;
                            let v_multj = v_mult[j] / deltav;
                            if ((up - u_multj) * (uv1 - u_multj) + (vp - v_multj) * (uv2 - v_multj) < 0.0)
                                || ((uv1 - u_multj).abs() < tolu && (uv2 - v_multj).abs() < tolv)
                            {
                                *irang = i as i32;
                                break;
                            }
                        }
                        j += 1;
                    }
                }
            }
        }
        arrive
    }

    // =====================================================================
    // TestArretAjout (gxx L965-1031)
    // =====================================================================
    fn test_arret_ajout(
        &mut self,
        func: &mut SurfFunction,
        uv: &mut [f64; 2],
        irang: &mut i32,
        psol: &mut PntOn2S,
    ) -> bool {
        let mut arrive = false;
        let (up, vp) = if !self.reversed {
            self.previous_point.parameters_on_surface(false)
        } else {
            self.previous_point.parameters_on_surface(true)
        };

        let nb_ajout = self.seq_ajout.len();
        for i in 1..=nb_ajout {
            *irang = self.seq_ajout[i - 1];
            if (*irang as isize).unsigned_abs() as usize <= self.lines.len() {
                let line = &self.lines[(*irang as isize).unsigned_abs() as usize - 1];
                let p;
                if *irang > 0 {
                    p = line.value(line.nb_points()).clone();
                } else {
                    p = line.value(1).clone();
                }
                let (u1, v1) = if !self.reversed {
                    p.parameters_on_surface(false)
                } else {
                    p.parameters_on_surface(true)
                };
                if ((up - u1) * (uv[0] - u1) + (vp - v1) * (uv[1] - v1)) < 0.0
                    || ((uv[0] - u1).abs() < self.tolerance[0]
                        && (uv[1] - v1).abs() < self.tolerance[1])
                {
                    arrive = true;
                    uv[0] = u1;
                    uv[1] = v1;
                    let _ = func.values(&[uv[0], uv[1]]);
                    *psol = p;
                    break;
                }
            }
        }
        arrive
    }

    // =====================================================================
    // FillPntsInHoles (gxx L1033-1181)
    // =====================================================================
    fn fill_pnts_in_holes(
        &mut self,
        func: &mut SurfFunction,
        mut copy_seq_alone: Vec<i32>,
        pnts_in_holes: &mut Vec<InteriorPoint>,
    ) {
        let born_inf = [self.um, self.vm];
        let born_sup = [self.um_max, self.vm_max];
        self.point_line_line.clear();
        let mut seq_to_remove: Vec<i32> = Vec::new();
        let mut bad_solutions: Vec<i32> = Vec::new();

        let mut i: usize = 1;
        while i < copy_seq_alone.len() {
            let irang1 = copy_seq_alone[i - 1];
            if irang1 == 0 {
                i += 1;
                continue;
            }
            let mut to_remove = false;
            let mut point_alone1 = PntOn2S::new();
            let mut point_alone2 = PntOn2S::new();
            let line1 = &self.lines[(irang1.abs() as usize) - 1];
            if irang1 > 0 {
                point_alone1 = line1.value(line1.nb_points()).clone();
            } else {
                point_alone1 = line1.value(1).clone();
            }
            let p2d1 = point_alone1.value_on_surface(self.reversed);
            let mut min_sq_dist = f64::MAX;
            let mut min_rang = 0i32;
            let mut min_index = 0usize;
            for j in (i + 1)..=copy_seq_alone.len() {
                let irang2 = copy_seq_alone[j - 1];
                if irang2 == 0 || bad_solutions.contains(&irang2) {
                    continue;
                }
                let line2 = &self.lines[(irang2.abs() as usize) - 1];
                if irang2 > 0 {
                    point_alone2 = line2.value(line2.nb_points()).clone();
                } else {
                    point_alone2 = line2.value(1).clone();
                }
                let p2d2 = point_alone2.value_on_surface(self.reversed);
                let a_sq_dist = p2d1.distance_squared(p2d2);
                if a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                    min_rang = irang2;
                    min_index = j;
                }
            }
            if min_rang == 0 {
                seq_to_remove.push(irang1);
                bad_solutions.clear();
                i += 1;
                continue;
            }
            // Ends of same line.
            if irang1.abs() == min_rang.abs()
                && self.lines[(irang1.abs() as usize) - 1].nb_points() == 2
            {
                seq_to_remove.push(irang1);
                seq_to_remove.push(min_rang);
                copy_seq_alone[i - 1] = 0;
                copy_seq_alone[min_index - 1] = 0;
                bad_solutions.clear();
                i += 1;
                continue;
            }

            let line2 = &self.lines[(min_rang.abs() as usize) - 1];
            if min_rang > 0 {
                point_alone2 = line2.value(line2.nb_points()).clone();
            } else {
                point_alone2 = line2.value(1).clone();
            }
            let pnt1 = point_alone1.value();
            let pnt2 = point_alone2.value();
            let p2d2 = point_alone2.value_on_surface(self.reversed);
            let min_sq_dist_3d = pnt1.distance_squared(pnt2);
            if min_sq_dist_3d <= self.epsilon
                || ((p2d1.x - p2d2.x).abs() <= self.tolerance[0]
                    && (p2d1.y - p2d2.y).abs() <= self.tolerance[1])
            {
                to_remove = true;
            } else {
                // Real curve.
                let uvap = [
                    (p2d1.x + p2d2.x) / 2.0,
                    (p2d1.y + p2d2.y) / 2.0,
                ];
                let mut rs_nld = FunctionSetRoot::new(func, self.tolerance);
                rs_nld.perform(func, uvap, born_inf, born_sup);
                if rs_nld.is_done() && func.root().abs() <= func.tolerance() && !func.is_tangent() {
                    let uv = rs_nld.root();
                    let pmid = DVec2::new(uv[0], uv[1]);
                    let p1p2 = p2d2 - p2d1;
                    let p1pmid = pmid - p2d1;
                    let p2pmid = pmid - p2d2;
                    let scal_prod1 = p1p2.dot(p1pmid);
                    let scal_prod2 = p1p2.dot(p2pmid);
                    let mut is_pmid_valid = scal_prod1 > 0.0 && scal_prod2 < 0.0;
                    if is_pmid_valid {
                        for iline in 1..=self.lines.len() {
                            if self.is_point_on_line_2d(pmid, iline as i32) {
                                is_pmid_valid = false;
                                break;
                            }
                        }
                    }
                    if is_pmid_valid {
                        let a_point = InteriorPoint::new_full(
                            func.point(),
                            uv[0],
                            uv[1],
                            func.direction_3d(),
                            func.direction_2d(),
                        );
                        pnts_in_holes.push(a_point);
                        self.point_line_line
                            .entry(pnts_in_holes.len() as i32)
                            .or_default()
                            .extend_from_slice(&[irang1, min_rang]);
                    } else {
                        if !bad_solutions.contains(&min_rang) {
                            bad_solutions.push(min_rang);
                        }
                        continue;
                    }
                } else {
                    if !bad_solutions.contains(&min_rang) {
                        bad_solutions.push(min_rang);
                    }
                    continue;
                }
            }
            copy_seq_alone[i - 1] = 0;
            copy_seq_alone[min_index - 1] = 0;
            if to_remove {
                seq_to_remove.push(irang1);
                seq_to_remove.push(min_rang);
            }
            bad_solutions.clear();
            i += 1;
        }

        for s in &seq_to_remove {
            let mut j = 1;
            while j <= self.seq_alone.len() {
                if self.seq_alone[j - 1] == *s {
                    self.seq_alone.remove(j - 1);
                    break;
                }
                j += 1;
            }
        }
    }

    // =====================================================================
    // TestArretCadre (gxx L1183-1397)
    // =====================================================================
    fn test_arret_cadre(
        &mut self,
        u_mult: &[f64],
        v_mult: &[f64],
        line: &mut IWLine,
        func: &mut SurfFunction,
        uv: &mut [f64; 2],
        irang: &mut i32,
    ) {
        let mut found = false;
        *irang = 0;
        for i in 1..self.wd1.len() {
            if self.wd1[i].etat < 0 {
                let mut n: usize = 0;
                if self.nb_multiplicities[i] > 0 {
                    for k in 1..i {
                        n += self.nb_multiplicities[k] as usize;
                    }
                }
                let (up0, vp0) = if !self.reversed {
                    line.value(1).parameters_on_surface(false)
                } else {
                    line.value(1).parameters_on_surface(true)
                };
                let mut up = up0;
                let mut vp = vp0;
                let nbp = line.nb_points();
                let mut j = 2;
                while j <= nbp {
                    let (uc, vc) = if !self.reversed {
                        line.value(j).parameters_on_surface(false)
                    } else {
                        line.value(j).parameters_on_surface(true)
                    };

                    let mut a_vec1 = DVec2::new(up - self.wd1[i].ustart, vp - self.wd1[i].vstart);
                    let mut a_vec2 = DVec2::new(uc - self.wd1[i].ustart, vc - self.wd1[i].vstart);
                    cut_vector_by_tolerances(&mut a_vec1, &self.tolerance);
                    cut_vector_by_tolerances(&mut a_vec2, &self.tolerance);

                    let scal = a_vec1.dot(a_vec2);

                    if scal < 0.0 {
                        line.cut(j);
                        *irang = i as i32;
                        uv[0] = self.wd1[*irang as usize].ustart;
                        uv[1] = self.wd1[*irang as usize].vstart;
                        found = true;
                    } else if (uc - self.wd1[i].ustart).abs() < self.tolerance[0]
                        && (vc - self.wd1[i].vstart).abs() < self.tolerance[1]
                    {
                        line.cut(j);
                        *irang = i as i32;
                        uv[0] = self.wd1[*irang as usize].ustart;
                        uv[1] = self.wd1[*irang as usize].vstart;
                        found = true;
                    } else if self.nb_multiplicities[i] > 0 {
                        let mut k = n;
                        while k < n + self.nb_multiplicities[i] as usize {
                            if k >= u_mult.len() {
                                break;
                            }
                            a_vec1 = DVec2::new(up - u_mult[k], vp - v_mult[k]);
                            a_vec2 = DVec2::new(uc - u_mult[k], vc - v_mult[k]);
                            cut_vector_by_tolerances(&mut a_vec1, &self.tolerance);
                            cut_vector_by_tolerances(&mut a_vec2, &self.tolerance);

                            let scal = a_vec1.dot(a_vec2);
                            if scal < 0.0 {
                                line.cut(j);
                                *irang = i as i32;
                                uv[0] = self.wd1[*irang as usize].ustart;
                                uv[1] = self.wd1[*irang as usize].vstart;
                                found = true;
                                break;
                            } else if (uc - u_mult[k]).abs() < self.tolerance[0]
                                && (vc - v_mult[k]).abs() < self.tolerance[1]
                            {
                                line.cut(j);
                                *irang = i as i32;
                                uv[0] = self.wd1[*irang as usize].ustart;
                                uv[1] = self.wd1[*irang as usize].vstart;
                                found = true;
                                break;
                            }
                            k += 1;
                        }
                    }
                    if found {
                        let _ = func.values(&[uv[0], uv[1]]);
                        let nbp = line.nb_points();
                        let (vcttg, indextg) = line.tangent_vector();
                        if indextg > nbp as i32 {
                            if j > 3 && j <= nbp + 1 {
                                let dir3d = func.direction_3d();
                                let mut dir3d1 = line.value(j - 1).value() - line.value(j - 2).value();
                                let dot = dir3d.dot(dir3d1);
                                if dot < 0.0 {
                                    // Normally this Function should not be used often!
                                    dir3d1 = -dir3d1;
                                }
                                let _ = (dir3d1, vcttg);
                                line.set_tangent_vector(func.direction_3d(), j as i32 - 1);
                            }
                        }
                        return;
                    }
                    up = uc;
                    vp = vc;
                    j += 1;
                }

                // Now the last point of the line and the last calculated point
                // are compared; there will be no need to "Cut".
                let mut a_vec1 = DVec2::new(up - self.wd1[i].ustart, vp - self.wd1[i].vstart);
                let mut a_vec2 = DVec2::new(uv[0] - self.wd1[i].ustart, uv[1] - self.wd1[i].vstart);
                cut_vector_by_tolerances(&mut a_vec1, &self.tolerance);
                cut_vector_by_tolerances(&mut a_vec2, &self.tolerance);

                let scal = a_vec1.dot(a_vec2);

                if scal < 0.0 {
                    *irang = i as i32;
                    uv[0] = self.wd1[*irang as usize].ustart;
                    uv[1] = self.wd1[*irang as usize].vstart;
                    found = true;
                } else if (uv[0] - self.wd1[i].ustart).abs() < self.tolerance[0]
                    && (uv[1] - self.wd1[i].vstart).abs() < self.tolerance[1]
                {
                    *irang = i as i32;
                    uv[0] = self.wd1[*irang as usize].ustart;
                    uv[1] = self.wd1[*irang as usize].vstart;
                    found = true;
                } else if self.nb_multiplicities[i] > 0 {
                    let mut j = n;
                    while j < n + self.nb_multiplicities[i] as usize {
                        if j >= u_mult.len() {
                            break;
                        }
                        a_vec1 = DVec2::new(up - u_mult[j], vp - v_mult[j]);
                        a_vec2 = DVec2::new(uv[0] - u_mult[j], uv[1] - v_mult[j]);
                        cut_vector_by_tolerances(&mut a_vec1, &self.tolerance);
                        cut_vector_by_tolerances(&mut a_vec2, &self.tolerance);

                        let scal = a_vec1.dot(a_vec2);
                        if scal < 0.0 {
                            *irang = i as i32;
                            uv[0] = self.wd1[*irang as usize].ustart;
                            uv[1] = self.wd1[*irang as usize].vstart;
                            found = true;
                            break;
                        } else if (uv[0] - u_mult[j]).abs() < self.tolerance[0]
                            && (uv[1] - v_mult[j]).abs() < self.tolerance[1]
                        {
                            *irang = i as i32;
                            uv[0] = self.wd1[*irang as usize].ustart;
                            uv[1] = self.wd1[*irang as usize].vstart;
                            found = true;
                            break;
                        }
                        j += 1;
                    }
                }
                if found {
                    *irang = -*irang;
                    let _ = func.values(&[uv[0], uv[1]]);
                    return;
                }
            }
        }
    }

    // =====================================================================
    // TestDeflection (gxx L2628-2899)
    // =====================================================================
    fn test_deflection(
        &mut self,
        func: &mut SurfFunction,
        finished: bool,
        uv: &[f64; 2],
        status_precedent: StatusDeflection,
        nb_division: &mut i32,
        step: &mut f64,
        step_sign: i32,
    ) -> StatusDeflection {
        let mut a_status = StatusDeflection::OK;

        let (paramu, paramv) = if !self.reversed {
            self.previous_point.parameters_on_surface(false)
        } else {
            self.previous_point.parameters_on_surface(true)
        };

        let du = uv[0] - paramu;
        let dv = uv[1] - paramv;
        let duv = du * du + dv * dv;

        // OCCT: Corde = gp_Vec(previousPoint.Value(), sp.Point()) = sp.Point() -
        // previousPoint.Value() = current - previous.
        let corde = func.point() - self.previous_point.value();
        let norme = corde.length_squared();

        if (norme <= 4.0 * rcad_kernel::precision::SQUARE_CONFUSION)
            && ((duv <= rcad_kernel::precision::square_p_confusion())
                || (status_precedent != StatusDeflection::OK))
        {
            a_status = StatusDeflection::PointConfondu;
            if status_precedent == StatusDeflection::PasTropGrand {
                return StatusDeflection::ArretSurPointPrecedent;
            }
        } else {
            let mut cosi = corde.dot(self.previous_d3d);
            let mut cosi2 = 0.0;

            if cosi * step_sign as f64 >= 0.0 {
                // angle 3d <= pi/2.
                let a_div = self.previous_d3d.length_squared() * norme;
                if a_div == 0.0 {
                    return a_status;
                }
                cosi2 = cosi * cosi / a_div;
            }
            if cosi2 < COS_REF_3D {
                // angle 3d too great.
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x).abs();
                let step_v = (*step * self.previous_d2d.y).abs();
                if step_u < self.tolerance[0] && step_v < self.tolerance[1] {
                    a_status = StatusDeflection::ArretSurPointPrecedent;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                return a_status;
            }
        }

        let a_min_tol_u = 0.1 * (*step * self.previous_d2d.x).abs();
        let a_min_tol_v = 0.1 * (*step * self.previous_d2d.y).abs();
        let a_tol_u = if a_min_tol_u > 0.0 {
            self.tolerance[0].min(a_min_tol_u)
        } else {
            self.tolerance[0]
        };
        let a_tol_v = if a_min_tol_v > 0.0 {
            self.tolerance[1].min(a_min_tol_v)
        } else {
            self.tolerance[1]
        };

        if du.abs() < a_tol_u && dv.abs() < a_tol_v {
            return StatusDeflection::ArretSurPointPrecedent; // confused point 2d
        }

        let mut cosi = step_sign as f64 * (du * self.previous_d2d.x + dv * self.previous_d2d.y);

        if cosi < 0.0 && a_status == StatusDeflection::PointConfondu {
            return StatusDeflection::ArretSurPointPrecedent; // leave as step back with confused point
        }

        if func.is_tangent() {
            return StatusDeflection::ArretSurPoint;
        }

        if (*nb_division < MAX_DIVISION)
            && (a_status != StatusDeflection::PointConfondu)
            && (status_precedent != StatusDeflection::PointConfondu)
        {
            let mut cosi2 = cosi * cosi / duv;
            if cosi2 < COS_REF_2D || cosi < 0.0 {
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x).abs();
                let step_v = (*step * self.previous_d2d.y).abs();
                if step_u < self.tolerance[0] && step_v < self.tolerance[1] {
                    a_status = StatusDeflection::ArretSurPointPrecedent;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                *nb_division += 1;
                return a_status;
            }

            cosi = corde.dot(func.direction_3d());
            let dir3d_sq = func.direction_3d().length_squared();
            cosi2 = if dir3d_sq > 0.0 && norme > 0.0 {
                cosi * cosi / dir3d_sq / norme
            } else {
                0.0
            };
            if cosi2 < COS_REF_3D {
                // angle 3d too great.
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x).abs();
                let step_v = (*step * self.previous_d2d.y).abs();
                if step_u < self.tolerance[0] && step_v < self.tolerance[1] {
                    a_status = StatusDeflection::ArretSurPoint;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                return a_status;
            }
            let d2d = func.direction_2d();
            cosi = du * d2d.x + dv * d2d.y;
            cosi2 = cosi * cosi / duv;
            if cosi2 < COS_REF_2D || d2d.dot(self.previous_d2d) < 0.0 {
                // angle 2d too great or change the side.
                *step /= 2.0;
                let step_u = (*step * self.previous_d2d.x).abs();
                let step_v = (*step * self.previous_d2d.y).abs();
                if step_u < self.tolerance[0] && step_v < self.tolerance[1] {
                    a_status = StatusDeflection::ArretSurPointPrecedent;
                } else {
                    a_status = StatusDeflection::PasTropGrand;
                }
                return a_status;
            }
        }

        if !finished {
            if a_status == StatusDeflection::PointConfondu {
                let step_u = (1.5 * du).abs().min(self.pas * (self.um_max - self.um));
                let step_v = (1.5 * dv).abs().min(self.pas * (self.vm_max - self.vm));

                let d2dx = self.previous_d2d.x.abs();
                let d2dy = self.previous_d2d.y.abs();

                if d2dx < self.tolerance[0] {
                    *step = if d2dy != 0.0 { step_v / d2dy } else { *step };
                } else if d2dy < self.tolerance[1] {
                    *step = if d2dx != 0.0 { step_u / d2dx } else { *step };
                } else {
                    *step = if d2dx != 0.0 && d2dy != 0.0 {
                        (step_u / d2dx).min(step_v / d2dy)
                    } else {
                        *step
                    };
                }
            } else {
                let fleche_courante = (self.previous_d3d.normalize_or_zero()
                    - func.direction_3d().normalize_or_zero())
                .length_squared()
                    * norme
                    / 64.0;

                if fleche_courante <= 0.25 * self.fleche * self.fleche {
                    let d2dx = func.direction_2d().x.abs();
                    let d2dy = func.direction_2d().y.abs();

                    let step_u = (1.5 * du).abs().min(self.pas * (self.um_max - self.um));
                    let step_v = (1.5 * dv).abs().min(self.pas * (self.vm_max - self.vm));

                    if d2dx < self.tolerance[0] {
                        *step = if d2dy != 0.0 { step_v / d2dy } else { *step };
                    } else if d2dy < self.tolerance[1] {
                        *step = if d2dx != 0.0 { step_u / d2dx } else { *step };
                    } else {
                        *step = if d2dx != 0.0 && d2dy != 0.0 {
                            (step_u / d2dx).min(step_v / d2dy)
                        } else {
                            *step
                        };
                    }
                } else if fleche_courante > self.fleche * self.fleche {
                    // step too great.
                    *step /= 2.0;
                    let step_u = (*step * self.previous_d2d.x).abs();
                    let step_v = (*step * self.previous_d2d.y).abs();
                    if step_u < self.tolerance[0] && step_v < self.tolerance[1] {
                        a_status = StatusDeflection::ArretSurPointPrecedent;
                    } else {
                        a_status = StatusDeflection::PasTropGrand;
                    }
                } else {
                    let d2dx = func.direction_2d().x.abs();
                    let d2dy = func.direction_2d().y.abs();

                    let step_u = (1.5 * du).abs().min(self.pas * (self.um_max - self.um));
                    let step_v = (1.5 * dv).abs().min(self.pas * (self.vm_max - self.vm));

                    if d2dx < self.tolerance[0] {
                        if d2dy != 0.0 {
                            *step = (*step).min(step_v / d2dy);
                        }
                    } else if d2dy < self.tolerance[1] {
                        if d2dx != 0.0 {
                            *step = (*step).min(step_u / d2dx);
                        }
                    } else if d2dx != 0.0 && d2dy != 0.0 {
                        *step = (*step).min((step_u / d2dx).min(step_v / d2dy));
                    }
                }
            }
        }
        a_status
    }

    // =====================================================================
    // ComputeOpenLine (gxx L1414-1928)
    // =====================================================================
    fn compute_open_line(
        &mut self,
        u_mult: &[f64],
        v_mult: &[f64],
        pnts1: &[PathPoint],
        func: &mut SurfFunction,
        rajout: &mut bool,
    ) {
        let mut i: usize = 1;
        let mut n = 0i32;
        let mut born_inf = [self.um, self.vm];
        let mut born_sup = [self.um_max, self.vm_max];
        let mut uvap = [0.0f64; 2];
        let mut pas_c = 0.0f64;
        let mut pas_cu = 0.0f64;
        let mut pas_cv = 0.0f64;
        let mut arrive = false;
        let mut cadre = false;
        let mut arret_ajout = false;
        let mut psol = PntOn2S::new();
        let mut current_line: Option<IWLine> = None;
        let mut tgtend = false;
        let mut a_status = StatusDeflection::OK;
        let mut status_precedent = StatusDeflection::OK;
        let mut nb_division = 0i32;
        let mut step_sign = 0i32;
        let mut path_pnt = PathPoint::new();

        born_inf[0] = self.um;
        born_sup[0] = self.um_max;
        born_inf[1] = self.vm;
        born_sup[1] = self.vm_max;

        let mut rs_nld = FunctionSetRoot::new(func, self.tolerance);
        let nb_path = pnts1.len();

        let mut movementdirectioninfo = vec![0i32; nb_path + 1];

        let mut a_func_for_duplicate = func.clone();

        while i <= nb_path {
            if (self.wd1[i].etat > 11)
                || ((self.wd1[i].etat < -11) && (movementdirectioninfo[i] != 0))
            {
                path_pnt = pnts1[i - 1].clone();
                uvap[0] = self.wd1[i].ustart;
                uvap[1] = self.wd1[i].vstart;
                make_walking_point(
                    self.reversed,
                    11,
                    uvap[0],
                    uvap[1],
                    func,
                    &mut self.previous_point,
                );

                let prev = self.previous_point.clone();
                if self.is_point_on_line(
                    &prev,
                    &born_inf,
                    &born_sup,
                    &mut rs_nld,
                    &mut a_func_for_duplicate,
                ) {
                    self.wd1[i].etat = -self.wd1[i].etat.abs();
                    i += 1;
                    continue;
                }

                let mut cl = IWLine::new();
                cl.set_tangency_at_begining(false);
                tgtend = false;
                cl.add_status_first_full(false, true, i as i32, &path_pnt);
                self.previous_d3d = func.direction_3d();
                self.previous_d2d = func.direction_2d();
                cl.add_point(&self.previous_point);

                if movementdirectioninfo[i] != 0 {
                    if movementdirectioninfo[i] < 0 {
                        step_sign = -1;
                        cl.set_tangent_vector(-self.previous_d3d, 1);
                    } else {
                        step_sign = 1;
                        cl.set_tangent_vector(self.previous_d3d, 1);
                    }
                } else {
                    let tyutuyt = path_pnt.direction_3d().dot(self.previous_d3d);
                    if tyutuyt < 0.0 {
                        step_sign = -1;
                        cl.set_tangent_vector(-self.previous_d3d, 1);
                    } else {
                        step_sign = 1;
                        cl.set_tangent_vector(self.previous_d3d, 1);
                    }
                }

                self.wd1[i].etat = -self.wd1[i].etat.abs();
                movementdirectioninfo[i] = if movementdirectioninfo[i] == 0 {
                    step_sign
                } else {
                    0
                };

                // First step of advancement.
                let d2dx = self.previous_d2d.x.abs();
                let d2dy = self.previous_d2d.y.abs();
                if d2dx < self.tolerance[0] {
                    pas_c = if d2dy != 0.0 {
                        self.pas * (self.vm_max - self.vm) / d2dy
                    } else {
                        self.pas * (self.vm_max - self.vm)
                    };
                } else if d2dy < self.tolerance[1] {
                    pas_c = if d2dx != 0.0 {
                        self.pas * (self.um_max - self.um) / d2dx
                    } else {
                        self.pas * (self.um_max - self.um)
                    };
                } else {
                    // OCCT: pas * min((UM-Um)/d2dx, (VM-Vm)/d2dy).
                    pas_c = self.pas
                        * ((self.um_max - self.um) / d2dx).min((self.vm_max - self.vm) / d2dy);
                }

                arrive = false;
                arret_ajout = false;
                nb_division = 0;
                status_precedent = StatusDeflection::OK;
                let mut index_of_path_point_do_not_check = 0usize;
                let mut a_nb_iter = 10;
                let mut a_nb_bad_root_iter = 0i32;

                while !arrive {
                    cadre = self.cadrage(&mut born_inf, &mut born_sup, &mut uvap, &mut pas_c, step_sign);
                    rs_nld.perform(func, uvap, born_inf, born_sup);

                    if cadre {
                        born_inf[0] = self.um;
                        born_sup[0] = self.um_max;
                        born_inf[1] = self.vm;
                        born_sup[1] = self.vm_max;
                    }
                    if rs_nld.is_done() {
                        if func.root().abs() > func.tolerance() {
                            pas_c /= 2.0;
                            pas_cu = (pas_c * self.previous_d2d.x).abs();
                            pas_cv = (pas_c * self.previous_d2d.y).abs();
                            a_nb_bad_root_iter += 1;
                            if (pas_cu <= self.tolerance[0] && pas_cv <= self.tolerance[1])
                                || a_nb_bad_root_iter > MAX_DIVISION
                            {
                                if cl.nb_points() == 1 {
                                    break;
                                }
                                arrive = true;
                                cl.add_status_last(false);
                                tgtend = true;
                                *rajout = true;
                                self.seq_alone.push(self.lines.len() as i32 + 1);
                                self.seq_ajout.push(self.lines.len() as i32 + 1);
                            }
                        } else {
                            // test stop.
                            a_nb_bad_root_iter = 0;
                            let r = rs_nld.root();
                            uvap[0] = r[0];
                            uvap[1] = r[1];
                            arrive =
                                self.test_arret_passage_open(u_mult, v_mult, func, &mut uvap, &mut n);
                            if arrive {
                                cadre = false;
                            } else {
                                if *rajout {
                                    arret_ajout =
                                        self.test_arret_ajout(func, &mut uvap, &mut n, &mut psol);
                                    if arret_ajout {
                                        tgtend = self.lines[n as usize].is_tangent_at_end();
                                        n = -n;
                                    }
                                }
                                if !(*rajout && arret_ajout) {
                                    let (prev_up, prev_vp) = if !self.reversed {
                                        self.previous_point.parameters_on_surface(false)
                                    } else {
                                        self.previous_point.parameters_on_surface(true)
                                    };
                                    arrive = test_passed_solution_with_negative_state(
                                        &self.wd1,
                                        u_mult,
                                        v_mult,
                                        prev_up,
                                        prev_vp,
                                        &self.nb_multiplicities,
                                        &self.tolerance,
                                        func,
                                        &mut uvap,
                                        &mut n,
                                    );
                                    if arrive {
                                        cadre = false;
                                    }
                                }
                                if !arret_ajout && cadre {
                                    if cl.nb_points() == 1 {
                                        break; // cancel the line
                                    }
                                    self.test_arret_cadre(
                                        u_mult, v_mult, &mut cl, func, &mut uvap, &mut n,
                                    );
                                    if n <= 0 {
                                        self.make_walking_point(2, uvap[0], uvap[1], func, &mut psol);
                                        tgtend = func.is_tangent();
                                        n = -n;
                                    }
                                }
                            }
                            a_status = self.test_deflection(
                                func,
                                arrive,
                                &uvap,
                                status_precedent,
                                &mut nb_division,
                                &mut pas_c,
                                step_sign,
                            );
                            status_precedent = a_status;
                            if a_status == StatusDeflection::PasTropGrand {
                                arrive = false;
                                arret_ajout = false;
                                tgtend = false;
                                if !self.reversed {
                                    let (_, _) = self.previous_point.parameters_on_surface(false);
                                    self.previous_point.set_value_uv(false, uvap[0], uvap[1]);
                                } else {
                                    self.previous_point.set_value_uv(true, uvap[0], uvap[1]);
                                }
                            } else if arret_ajout || cadre {
                                arrive = true;
                                cl.add_status_last(false);
                                cl.add_point(&psol);
                                // Remove <n> from <seq_alone>.
                                let mut iseq = 1;
                                while iseq <= self.seq_alone.len() {
                                    if self.seq_alone[iseq - 1] == n {
                                        self.seq_alone.remove(iseq - 1);
                                        break;
                                    }
                                    iseq += 1;
                                }

                                if cadre && n == 0 {
                                    *rajout = true;
                                    self.seq_ajout.push(self.lines.len() as i32 + 1);
                                }
                            } else if a_status == StatusDeflection::ArretSurPointPrecedent {
                                if cl.nb_points() == 1 {
                                    // cancel the line.
                                    arrive = false;
                                    break;
                                }
                                arrive = true;
                                *rajout = true;
                                self.seq_alone.push(self.lines.len() as i32 + 1);
                                self.seq_ajout.push(self.lines.len() as i32 + 1);
                                cl.add_status_last(false);
                                tgtend = true;
                            } else if arrive {
                                if cl.nb_points() == 1
                                    && (n == i as i32 || a_status == StatusDeflection::PointConfondu)
                                {
                                    arrive = false;
                                    break;
                                }
                                // Point of stop given at input.
                                path_pnt = pnts1[n as usize - 1].clone();
                                let etat1_n = self.wd1[n as usize].etat;
                                if etat1_n.abs() < 11 {
                                    // Passing point that is a stop.
                                    if a_status == StatusDeflection::ArretSurPoint {
                                        cl.add_status_last(false);
                                        tgtend = true;
                                    } else {
                                        arrive = false;
                                    }
                                    cl.add_index_passing(n);
                                } else {
                                    // Point of stop given at input.
                                    if etat1_n == 11 {
                                        tgtend = true;
                                    }
                                    cl.add_status_last_full(true, n, &path_pnt);
                                }
                                self.add_point_in_current_line(n, &path_pnt, &mut cl);
                                if etat1_n != 1 && etat1_n != 11 {
                                    self.wd1[n as usize].etat = -etat1_n.abs();
                                    movementdirectioninfo[n as usize] =
                                        if movementdirectioninfo[n as usize] == 0 {
                                            step_sign
                                        } else {
                                            0
                                        };
                                    if arrive && movementdirectioninfo[n as usize] != 0 {
                                        index_of_path_point_do_not_check = n as usize;
                                    }
                                    if arrive {
                                        *rajout = true;
                                        self.seq_ajout.push(self.lines.len() as i32 + 1);
                                    }
                                }
                            } else if a_status == StatusDeflection::ArretSurPoint {
                                arrive = true;
                                cl.add_status_last(false);
                                tgtend = true;
                                self.make_walking_point(1, uvap[0], uvap[1], func, &mut psol);
                                cl.add_point(&psol);
                                *rajout = true;
                                self.seq_alone.push(self.lines.len() as i32 + 1);
                                self.seq_ajout.push(self.lines.len() as i32 + 1);
                            } else if a_status == StatusDeflection::OK {
                                make_walking_point(
                                    self.reversed,
                                    2,
                                    uvap[0],
                                    uvap[1],
                                    func,
                                    &mut self.previous_point,
                                );
                                self.previous_d3d = func.direction_3d();
                                self.previous_d2d = func.direction_2d();
                                cl.add_point(&self.previous_point);
                            } else if a_status == StatusDeflection::PointConfondu {
                                a_nb_iter -= 1;
                            }
                        }
                    } else {
                        // No numerical solution.
                        pas_c /= 2.0;
                        pas_cu = (pas_c * self.previous_d2d.x).abs();
                        pas_cv = (pas_c * self.previous_d2d.y).abs();
                        if pas_cu <= self.tolerance[0] && pas_cv <= self.tolerance[1] {
                            if cl.nb_points() == 1 {
                                break;
                            }
                            arrive = true;
                            cl.add_status_last(false);
                            tgtend = true;
                            *rajout = true;
                            self.seq_alone.push(self.lines.len() as i32 + 1);
                            self.seq_ajout.push(self.lines.len() as i32 + 1);
                        }
                    }

                    if a_nb_iter < 0 {
                        break;
                    }
                }

                if arrive {
                    cl.set_tangency_at_end(tgtend);
                    self.lines.push(cl);
                    movementdirectioninfo[i] = 0;
                    if self.wd1[i].etat > 0 {
                        self.wd1[i].etat = -self.wd1[i].etat;
                    }

                    // lbr le 5 juin 97 (Pb ds Contap).
                    for av in 1..=nb_path {
                        if (self.wd1[av].etat > 11)
                            || ((av != i)
                                && (av != index_of_path_point_do_not_check)
                                && (self.wd1[av].etat < -11)
                                && (movementdirectioninfo[av] != 0))
                        {
                            let mut uav = self.wd1[av].ustart;
                            let mut vav = self.wd1[av].vstart;
                            let av_p = self.lines.last().unwrap().value(self.lines.last().unwrap().nb_points()).clone();
                            let (uavp, vavp) = if !self.reversed {
                                av_p.parameters_on_surface(false)
                            } else {
                                av_p.parameters_on_surface(true)
                            };
                            uav -= uavp;
                            vav -= vavp;
                            uav *= 0.001;
                            vav *= 0.001;
                            if uav.abs() < self.tolerance[0] && vav.abs() < self.tolerance[1] {
                                if self.wd1[av].etat < 0 {
                                    movementdirectioninfo[av] = 0;
                                } else {
                                    self.wd1[av].etat = -self.wd1[av].etat;
                                    movementdirectioninfo[av] = step_sign;
                                }
                                let last_line = self.lines.last_mut().unwrap();
                                last_line.add_status_last_full(true, av as i32, &pnts1[av - 1].clone());
                            }

                            let av_pp = self.lines.last().unwrap().value(1).clone();
                            let (uavp, vavp) = if !self.reversed {
                                av_pp.parameters_on_surface(false)
                            } else {
                                av_pp.parameters_on_surface(true)
                            };
                            uav = self.wd1[av].ustart;
                            vav = self.wd1[av].vstart;
                            uav -= uavp;
                            vav -= vavp;
                            uav *= 0.001;
                            vav *= 0.001;
                            if uav.abs() < self.tolerance[0] && vav.abs() < self.tolerance[1] {
                                if self.wd1[av].etat < 0 {
                                    movementdirectioninfo[av] = 0;
                                } else {
                                    self.wd1[av].etat = -self.wd1[av].etat;
                                    movementdirectioninfo[av] = -step_sign;
                                }
                                let last_line = self.lines.last_mut().unwrap();
                                last_line.add_status_first_full(false, true, av as i32, &pnts1[av - 1].clone());
                            }
                        }
                    }
                }
            }
            i += 1;
        }
        let _ = (pas_cu, pas_cv, current_line);
    }

    // =====================================================================
    // ComputeCloseLine (gxx L2007-2624)
    // =====================================================================
    fn compute_close_line(
        &mut self,
        u_mult: &[f64],
        v_mult: &[f64],
        pnts1: &[PathPoint],
        pnts2: &[InteriorPoint],
        func: &mut SurfFunction,
        rajout: &mut bool,
    ) {
        let mut i: usize = 1;
        let mut n = 0i32;
        let mut born_inf = [self.um, self.vm];
        let mut born_sup = [self.um_max, self.vm_max];
        let mut uvap = [0.0f64; 2];
        let mut pas_c = 0.0f64;
        let mut pas_cu = 0.0f64;
        let mut pas_cv = 0.0f64;
        let mut pas_sav = 0.0f64;
        let mut arrive = false;
        let mut cadre = false;
        let mut arret_ajout = false;
        let mut psol = PntOn2S::new();
        let mut current_line: Option<IWLine> = None;
        let mut path_pnt = PathPoint::new();
        let mut loop_pnt = InteriorPoint::new();

        let mut tgtbeg = false;
        let mut tgtend = false;
        let mut step_sign = 0i32;

        let mut a_status = StatusDeflection::OK;
        let mut status_precedent;
        let mut nb_division = 0i32;

        let mut ipass = 0i32;

        born_inf[0] = self.um;
        born_sup[0] = self.um_max;
        born_inf[1] = self.vm;
        born_sup[1] = self.vm_max;

        let mut rs_nld = FunctionSetRoot::new(func, self.tolerance);
        let nb_loop = pnts2.len();

        // Check borders for degeneracy.
        let a_nb_sample_pnt = 10;
        let mut is_left_degenerated_border = [true, true];
        let mut is_right_degenerated_border = [true, true];
        let a_step = [
            (born_sup[0] - born_inf[0]) / (a_nb_sample_pnt - 1) as f64,
            (born_sup[1] - born_inf[1]) / (a_nb_sample_pnt - 1) as f64,
        ];
        for a_border_idx in 1..=2 {
            let a_change_idx = if a_border_idx == 2 { 1 } else { 2 };
            let mut uv = [0.0f64; 2];

            // Left border.
            uv[a_border_idx - 1] = born_inf[a_border_idx - 1];
            for a_param_idx in 0..a_nb_sample_pnt {
                let a_param = born_inf[a_change_idx - 1] + a_param_idx as f64 * a_step[a_change_idx - 1];
                uv[a_change_idx - 1] = a_param;
                let d = func.derivatives(&uv);
                if d.is_some() {
                    let d = d.unwrap();
                    if d[a_change_idx - 1].abs() > rcad_kernel::precision::CONFUSION {
                        is_left_degenerated_border[a_border_idx - 1] = false;
                        break;
                    }
                }
            }

            // Right border.
            uv[a_border_idx - 1] = born_sup[a_border_idx - 1];
            for a_param_idx in 0..a_nb_sample_pnt {
                let a_param = born_inf[a_change_idx - 1] + a_param_idx as f64 * a_step[a_change_idx - 1];
                uv[a_change_idx - 1] = a_param;
                let d = func.derivatives(&uv);
                if d.is_some() {
                    let d = d.unwrap();
                    if d[a_change_idx - 1].abs() > rcad_kernel::precision::CONFUSION {
                        is_right_degenerated_border[a_border_idx - 1] = false;
                        break;
                    }
                }
            }
        }

        let mut a_func_for_duplicate = func.clone();

        while i <= nb_loop {
            if self.wd2[i].etat > 12 {
                // Start point of closed line.
                loop_pnt = pnts2[i - 1].clone();
                self.previous_point.set_value(
                    loop_pnt.value(),
                    self.reversed,
                    self.wd2[i].ustart,
                    self.wd2[i].vstart,
                );

                let prev = self.previous_point.clone();
                if self.is_point_on_line(
                    &prev,
                    &born_inf,
                    &born_sup,
                    &mut rs_nld,
                    &mut a_func_for_duplicate,
                ) {
                    self.wd2[i].etat = -self.wd2[i].etat;
                    i += 1;
                    continue;
                }

                self.previous_d3d = loop_pnt.direction();
                self.previous_d2d = loop_pnt.direction_2d();

                let mut cl = IWLine::new();
                cl.add_point(&self.previous_point);
                cl.set_tangent_vector(self.previous_d3d, 1);
                tgtbeg = false;
                tgtend = false;
                uvap[0] = self.wd2[i].ustart;
                uvap[1] = self.wd2[i].vstart;

                step_sign = 1;

                // First step of advancement.
                let d2dx = self.previous_d2d.x.abs();
                let d2dy = self.previous_d2d.y.abs();
                if d2dx < self.tolerance[0] {
                    pas_c = if d2dy != 0.0 {
                        self.pas * (self.vm_max - self.vm) / d2dy
                    } else {
                        self.pas * (self.vm_max - self.vm)
                    };
                } else if d2dy < self.tolerance[1] {
                    pas_c = if d2dx != 0.0 {
                        self.pas * (self.um_max - self.um) / d2dx
                    } else {
                        self.pas * (self.um_max - self.um)
                    };
                } else {
                    pas_c = self.pas
                        * (self.um_max - self.um)
                            .min(self.vm_max - self.vm)
                        / (d2dx.max(d2dy));
                }

                pas_sav = pas_c;

                arrive = false;
                arret_ajout = false;
                nb_division = 0;
                status_precedent = StatusDeflection::OK;
                let mut a_nb_iter = 10;
                let mut a_nb_bad_root_iter = 0i32;
                while !arrive {
                    cadre = self.cadrage(&mut born_inf, &mut born_sup, &mut uvap, &mut pas_c, step_sign);
                    rs_nld.perform(func, uvap, born_inf, born_sup);
                    let mut is_on_degenerated_border = false;

                    if cadre {
                        born_inf[0] = self.um;
                        born_sup[0] = self.um_max;
                        born_inf[1] = self.vm;
                        born_sup[1] = self.vm_max;
                    }
                    if rs_nld.is_done() {
                        if func.root().abs() > func.tolerance() {
                            // No solution for the tolerance.
                            pas_c /= 2.0;
                            pas_cu = (pas_c * self.previous_d2d.x).abs();
                            pas_cv = (pas_c * self.previous_d2d.y).abs();
                            a_nb_bad_root_iter += 1;

                            if (pas_cu <= self.tolerance[0] && pas_cv <= self.tolerance[1])
                                || a_nb_bad_root_iter > MAX_DIVISION
                            {
                                if cl.nb_points() == 1 {
                                    self.remove_two_end_points(i as i32);
                                    break; // cancel the line
                                }
                                if self.wd2[i].etat > 12 {
                                    // The line should become open.
                                    self.wd2[i].etat = 12;
                                    arret_ajout = false;
                                    self.open_line(0, &psol, pnts1, func, &mut cl);
                                    step_sign = -1;
                                    status_precedent = StatusDeflection::OK;
                                    arrive = false;
                                    pas_c = pas_sav;
                                    a_nb_bad_root_iter = 0;
                                    *rajout = true;
                                    self.seq_alone.push(-(self.lines.len() as i32) - 1);
                                    self.seq_ajout.push(-(self.lines.len() as i32) - 1);
                                } else {
                                    // Line s is open.
                                    arrive = true;
                                    cl.add_status_last(false);
                                    *rajout = true;
                                    self.seq_alone.push(self.lines.len() as i32 + 1);
                                    self.seq_ajout.push(self.lines.len() as i32 + 1);
                                    tgtend = true;
                                }
                            }
                        } else {
                            // There is a solution.
                            a_nb_bad_root_iter = 0;
                            let r = rs_nld.root();
                            uvap[0] = r[0];
                            uvap[1] = r[1];

                            // Avoid uninitialized memory access.
                            if cl.nb_points() > 2 {
                                for a_coord_idx in 1..=2 {
                                    if (is_left_degenerated_border[a_coord_idx - 1]
                                        && (uvap[a_coord_idx - 1] - born_inf[a_coord_idx - 1]).abs()
                                            < rcad_kernel::precision::PCONFUSION)
                                        || (is_right_degenerated_border[a_coord_idx - 1]
                                            && (uvap[a_coord_idx - 1] - born_sup[a_coord_idx - 1]).abs()
                                                < rcad_kernel::precision::PCONFUSION)
                                    {
                                        let (uvprev0, uvprev1) = if !self.reversed {
                                            cl.value(cl.nb_points() - 1).parameters_on_surface(false)
                                        } else {
                                            cl.value(cl.nb_points() - 1).parameters_on_surface(true)
                                        };
                                        let (uv0, uv1) = if !self.reversed {
                                            cl.value(cl.nb_points()).parameters_on_surface(false)
                                        } else {
                                            cl.value(cl.nb_points()).parameters_on_surface(true)
                                        };
                                        let uvprev = [uvprev0, uvprev1];
                                        let uv_l = [uv0, uv1];

                                        let mut a_scale_coeff = 0.0;

                                        // Avoid finite cycle which leads to stop
                                        // computing iline.
                                        if a_status != StatusDeflection::PasTropGrand {
                                            if (uv_l[a_coord_idx - 1] - uvprev[a_coord_idx - 1]).abs()
                                                > f64::EPSILON
                                            {
                                                a_scale_coeff = ((uvap[a_coord_idx - 1]
                                                    - uv_l[a_coord_idx - 1])
                                                    / (uv_l[a_coord_idx - 1] - uvprev[a_coord_idx - 1]))
                                                .abs();
                                            }
                                            let a_fix_idx = if a_coord_idx == 1 { 2 } else { 1 };
                                            uvap[a_fix_idx - 1] = uv_l[a_fix_idx - 1]
                                                + (uv_l[a_fix_idx - 1] - uvprev[a_fix_idx - 1])
                                                    * a_scale_coeff;
                                            is_on_degenerated_border = true;
                                        }
                                    }
                                }
                            }

                            arrive = self.test_arret_passage_close(u_mult, v_mult, &uvap, i as i32, &mut ipass);
                            if arrive {
                                // Reset proper parameter to test the arrow.
                                psol = cl.value(1).clone();
                                if !self.reversed {
                                    psol.set_value_uv(false, uvap[0], uvap[1]);
                                } else {
                                    psol.set_value_uv(true, uvap[0], uvap[1]);
                                }
                                cadre = false;
                            } else {
                                if *rajout {
                                    // Test on added points.
                                    arret_ajout = self.test_arret_ajout(func, &mut uvap, &mut n, &mut psol);
                                    if arret_ajout {
                                        if n > 0 {
                                            tgtend = self.lines[n as usize].is_tangent_at_end();
                                            n = -n;
                                        } else {
                                            tgtend = self.lines[(-n) as usize].is_tangent_at_begining();
                                        }
                                        arrive = self.wd2[i].etat == 12;
                                    }
                                }

                                if !arret_ajout && cadre {
                                    if cl.nb_points() == 1 {
                                        self.remove_two_end_points(i as i32);
                                        break; // cancel the line
                                    }
                                    self.test_arret_cadre(u_mult, v_mult, &mut cl, func, &mut uvap, &mut n);
                                    if n <= 0 {
                                        self.make_walking_point(2, uvap[0], uvap[1], func, &mut psol);
                                        tgtend = func.is_tangent();
                                        if is_on_degenerated_border {
                                            tgtend = true;
                                        }
                                        n = -n;
                                    }
                                    arrive = self.wd2[i].etat == 12; // the line is open
                                }
                            }
                            a_status = self.test_deflection(
                                func,
                                arrive,
                                &uvap,
                                status_precedent,
                                &mut nb_division,
                                &mut pas_c,
                                step_sign,
                            );

                            if is_on_degenerated_border && tgtend {
                                a_status = StatusDeflection::ArretSurPoint;
                            }

                            status_precedent = a_status;
                            if a_status == StatusDeflection::PasTropGrand {
                                // Division of the step.
                                arrive = false;
                                arret_ajout = false;
                                tgtend = false;
                                if !self.reversed {
                                    self.previous_point.set_value_uv(false, uvap[0], uvap[1]);
                                } else {
                                    self.previous_point.set_value_uv(true, uvap[0], uvap[1]);
                                }
                            } else if arret_ajout || cadre {
                                if arrive {
                                    // Line s is open.
                                    cl.add_status_last(false);
                                    cl.add_point(&psol);

                                    // Remove <SaveN> from <seq_alone> and, if it is
                                    // first found point, from <seq_ajout> too.
                                    if self.is_valid_end_point(i as i32, n) {
                                        let mut iseq = 1;
                                        while iseq <= self.seq_alone.len() {
                                            if self.seq_alone[iseq - 1] == n {
                                                self.seq_alone.remove(iseq - 1);
                                                break;
                                            }
                                            iseq += 1;
                                        }
                                        if cl.nb_points() <= 3 {
                                            let mut iseq = 1;
                                            while iseq <= self.seq_ajout.len() {
                                                if self.seq_ajout[iseq - 1] == n {
                                                    self.seq_ajout.remove(iseq - 1);
                                                    break;
                                                }
                                                iseq += 1;
                                            }
                                        }
                                    } else {
                                        if self.seq_alone.last() == Some(&(-(self.lines.len() as i32) - 1)) {
                                            self.seq_alone.pop();
                                            self.seq_ajout.pop();
                                        }
                                        self.remove_two_end_points(i as i32);
                                        arrive = false;
                                        break; // cancel the line
                                    }

                                    if cadre && n == 0 {
                                        *rajout = true;
                                        self.seq_ajout.push(self.lines.len() as i32 + 1);
                                    }
                                } else {
                                    // Open.
                                    self.wd2[i].etat = 12; // declare it open
                                    tgtbeg = tgtend;
                                    tgtend = false;
                                    arret_ajout = false;
                                    step_sign = -1;
                                    status_precedent = StatusDeflection::OK;
                                    pas_c = pas_sav;
                                    // Check if <Psol> has been really updated.
                                    if arrive || *rajout || (!arret_ajout && cadre && n <= 0) {
                                        if a_status == StatusDeflection::ArretSurPointPrecedent {
                                            cl.add_point(&psol);
                                            self.open_line(0, &psol, pnts1, func, &mut cl);
                                        } else {
                                            self.open_line(-(self.lines.len() as i32) - 1, &psol, pnts1, func, &mut cl);
                                        }
                                    }
                                    // Remove <SaveN> from <seq_alone>.
                                    if self.is_valid_end_point(i as i32, n) {
                                        let mut iseq = 1;
                                        while iseq <= self.seq_alone.len() {
                                            if self.seq_alone[iseq - 1] == n {
                                                self.seq_alone.remove(iseq - 1);
                                                break;
                                            }
                                            iseq += 1;
                                        }
                                        if cl.nb_points() <= 2 {
                                            let mut iseq = 1;
                                            while iseq <= self.seq_ajout.len() {
                                                if self.seq_ajout[iseq - 1] == n {
                                                    self.seq_ajout.remove(iseq - 1);
                                                    break;
                                                }
                                                iseq += 1;
                                            }
                                        }
                                    } else {
                                        self.remove_two_end_points(i as i32);
                                        break; // cancel the line
                                    }

                                    if cadre && n == 0 {
                                        *rajout = true;
                                        self.seq_ajout.push(-(self.lines.len() as i32) - 1);
                                    }
                                }
                            } else if a_status == StatusDeflection::ArretSurPointPrecedent {
                                if cl.nb_points() == 1 {
                                    // Cancel the line.
                                    arrive = false;
                                    self.remove_two_end_points(i as i32);
                                    break;
                                }
                                if self.wd2[i].etat > 12 {
                                    // The line should become open.
                                    self.wd2[i].etat = 12;
                                    arret_ajout = false;
                                    self.open_line(0, &psol, pnts1, func, &mut cl);
                                    step_sign = -1;
                                    status_precedent = StatusDeflection::OK;
                                    arrive = false;
                                    pas_c = pas_sav;
                                    a_nb_bad_root_iter = 0;
                                    *rajout = true;
                                    self.seq_alone.push(-(self.lines.len() as i32) - 1);
                                    self.seq_ajout.push(-(self.lines.len() as i32) - 1);
                                } else {
                                    // Line s is open.
                                    arrive = true;
                                    cl.add_status_last(false);
                                    *rajout = true;
                                    self.seq_alone.push(self.lines.len() as i32 + 1);
                                    self.seq_ajout.push(self.lines.len() as i32 + 1);
                                }
                            } else if arrive {
                                if self.wd2[i].etat > 12 {
                                    // Line closed good case.
                                    cl.add_status_first_last(true, false, false);
                                    let first = cl.value(1).clone();
                                    cl.add_point(&first);
                                } else if (n > 0) && (pnts1.len() as i32 >= n) {
                                    // Point of stop given at input.
                                    path_pnt = pnts1[n as usize - 1].clone();
                                    cl.add_status_last_full(true, n, &path_pnt);
                                    self.add_point_in_current_line(n, &path_pnt, &mut cl);
                                }
                            } else if a_status == StatusDeflection::ArretSurPoint {
                                if self.wd2[i].etat > 12 {
                                    // Line should become open.
                                    self.wd2[i].etat = 12;
                                    tgtbeg = true;
                                    tgtend = false;
                                    n = -(self.lines.len() as i32) - 1;
                                    psol.set_value(func.point(), self.reversed, uvap[0], uvap[1]);
                                    self.open_line(n, &psol, pnts1, func, &mut cl);
                                    step_sign = -1;
                                    *rajout = true;
                                    self.seq_alone.push(n);
                                    self.seq_ajout.push(n);
                                    status_precedent = StatusDeflection::OK;
                                    arrive = false;
                                    pas_c = pas_sav;
                                    a_nb_bad_root_iter = 0;
                                } else {
                                    arrive = true;
                                    if ipass != 0 {
                                        // Point of passage, point of stop.
                                        path_pnt = pnts1[ipass as usize - 1].clone();
                                        cl.add_status_last_full(true, ipass, &path_pnt);
                                        self.add_point_in_current_line(ipass, &path_pnt, &mut cl);
                                    } else {
                                        cl.add_status_last(false);
                                        let mut new_p = PntOn2S::new();
                                        new_p.set_value(func.point(), self.reversed, uvap[0], uvap[1]);
                                        cl.add_point(&new_p);
                                        *rajout = true;
                                        self.seq_alone.push(self.lines.len() as i32 + 1);
                                        self.seq_ajout.push(self.lines.len() as i32 + 1);
                                    }
                                }
                            } else if a_status == StatusDeflection::OK {
                                if ipass != 0 {
                                    cl.add_index_passing(ipass);
                                }
                                self.previous_point
                                    .set_value(func.point(), self.reversed, uvap[0], uvap[1]);
                                self.previous_d3d = func.direction_3d();
                                self.previous_d2d = func.direction_2d();
                                cl.add_point(&self.previous_point);
                            } else if a_status == StatusDeflection::PointConfondu {
                                a_nb_iter -= 1;
                            }
                        }
                    } else {
                        // No numerical solution NotDone.
                        pas_c /= 2.0;
                        pas_cu = (pas_c * self.previous_d2d.x).abs();
                        pas_cv = (pas_c * self.previous_d2d.y).abs();

                        if pas_cu <= self.tolerance[0] && pas_cv <= self.tolerance[1] {
                            if cl.nb_points() == 1 {
                                self.remove_two_end_points(i as i32);
                                break; // cancel the line
                            }
                            if self.wd2[i].etat > 12 {
                                // The line should become open.
                                self.wd2[i].etat = 12;
                                arret_ajout = false;
                                self.open_line(0, &psol, pnts1, func, &mut cl);
                                step_sign = -1;
                                status_precedent = StatusDeflection::OK;
                                arrive = false;
                                pas_c = pas_sav;
                                a_nb_bad_root_iter = 0;
                                *rajout = true;
                                self.seq_alone.push(-(self.lines.len() as i32) - 1);
                                self.seq_ajout.push(-(self.lines.len() as i32) - 1);
                            } else {
                                // Line s is open.
                                arrive = true;
                                cl.add_status_last(false);
                                tgtend = true;
                                *rajout = true;
                                self.seq_alone.push(self.lines.len() as i32 + 1);
                                self.seq_ajout.push(self.lines.len() as i32 + 1);
                            }
                        }
                    }

                    if a_nb_iter < 0 {
                        break;
                    }
                }
                if arrive {
                    cl.set_tangency_at_begining(tgtbeg);
                    cl.set_tangency_at_end(tgtend);

                    self.lines.push(cl);
                    self.wd2[i].etat = -self.wd2[i].etat;
                }
            }
            i += 1;
        }
        let _ = (pas_cu, pas_cv, current_line, loop_pnt);
    }

    // =====================================================================
    // AddPointInCurrentLine (gxx L2908-2916)
    // =====================================================================
    fn add_point_in_current_line(&mut self, n: i32, path_pnt: &PathPoint, current_line: &mut IWLine) {
        let mut psol = PntOn2S::new();
        psol.set_value(
            path_pnt.value(),
            self.reversed,
            self.wd1[n as usize].ustart,
            self.wd1[n as usize].vstart,
        );
        current_line.add_point(&psol);
    }

    // =====================================================================
    // MakeWalkingPoint (gxx L2918-2951)
    // =====================================================================
    fn make_walking_point(&mut self, case: i32, u: f64, v: f64, func: &mut SurfFunction, psol: &mut PntOn2S) {
        make_walking_point(self.reversed, case, u, v, func, psol);
    }

    // =====================================================================
    // OpenLine (gxx L2953-2998)
    // =====================================================================
    fn open_line(
        &mut self,
        n: i32,
        psol: &PntOn2S,
        pnts1: &[PathPoint],
        func: &mut SurfFunction,
        line: &mut IWLine,
    ) {
        let mut uv = [0.0f64; 2];
        let first = line.value(1).clone();
        self.previous_point = first;
        if !self.reversed {
            let (u, v) = self.previous_point.parameters_on_surface(false);
            uv[0] = u;
            uv[1] = v;
        } else {
            let (u, v) = self.previous_point.parameters_on_surface(true);
            uv[0] = u;
            uv[1] = v;
        }
        let _ = func.values(&uv);
        self.previous_d3d = func.direction_3d();
        self.previous_d2d = func.direction_2d();

        if n > 0 {
            // Departure point given at input.
            let path_pnt = pnts1[n as usize - 1].clone();
            line.add_status_first_full(false, true, n, &path_pnt);
            self.add_point_in_current_line(n, &path_pnt, line);
        } else {
            if n < 0 {
                line.add_point(psol);
            }
            line.add_status_first(false, false);
            // Mark the line as open without a given stop point.
        }
        line.reverse();
        let (_, indextg) = line.tangent_vector();
        line.set_tangent_vector(-self.previous_d3d, line.nb_points() as i32);
        let _ = indextg;
    }

    // =====================================================================
    // IsValidEndPoint (gxx L3000-3013)
    // =====================================================================
    fn is_valid_end_point(&mut self, ind_of_point: i32, ind_of_line: i32) -> bool {
        if self.point_line_line.is_empty() {
            return true;
        }
        if let Some(list) = self.point_line_line.get_mut(&ind_of_point) {
            if let Some(pos) = list.iter().position(|v| *v == ind_of_line) {
                list.remove(pos);
                return true;
            }
        }
        false
    }

    // =====================================================================
    // RemoveTwoEndPoints (gxx L3015-3027)
    // =====================================================================
    fn remove_two_end_points(&mut self, ind_of_point: i32) {
        if let Some(list) = self.point_line_line.get(&ind_of_point) {
            let line1 = list[0];
            let line2 = list[list.len() - 1];
            let mut iseq = 1;
            while iseq <= self.seq_alone.len() {
                if self.seq_alone[iseq - 1] == line1 || self.seq_alone[iseq - 1] == line2 {
                    self.seq_alone.remove(iseq - 1);
                    iseq = iseq.saturating_sub(1);
                }
                iseq += 1;
            }
        }
    }

    // =====================================================================
    // IsPointOnLine(gp_Pnt2d, Irang) (gxx L3029-3048)
    // =====================================================================
    fn is_point_on_line_2d(&self, p2d: DVec2, irang: i32) -> bool {
        let a_line = &self.lines[irang.abs() as usize - 1];
        for i in 1..=a_line.nb_points() {
            let p2d1 = a_line.value(i).value_on_surface(self.reversed);
            if (p2d1.x - p2d.x).abs() <= self.tolerance[0]
                && (p2d1.y - p2d.y).abs() <= self.tolerance[1]
            {
                return true;
            }
            if i < a_line.nb_points() {
                let p2d2 = a_line.value(i + 1).value_on_surface(self.reversed);
                let pp1 = p2d - p2d1;
                let pp2 = p2d - p2d2;
                if pp1.dot(pp2) < 0.0 {
                    return true;
                }
            }
        }
        false
    }

    // =====================================================================
    // IsPointOnLine(IntSurf_PntOn2S, Binf, Bsup, Solver, Func) (gxx L3059-3152)
    // =====================================================================
    fn is_point_on_line(
        &mut self,
        p_on_2s: &PntOn2S,
        inf_bounds: &[f64; 2],
        sup_bounds: &[f64; 2],
        solver: &mut FunctionSetRoot,
        func: &mut SurfFunction,
    ) -> bool {
        let eps = 1.0;
        let a_p3d = p_on_2s.value();

        for a_l_idx in 1..=self.lines.len() {
            let a_l = &self.lines[a_l_idx - 1].line;

            if a_l.is_out_box(a_p3d) {
                continue;
            }

            // Look for the nearest segment.
            let mut a_umin = 0.0;
            let mut a_vmin = 0.0;
            let mut a_min_sq_dist = f64::MAX;
            for a_pt_idx in 1..a_l.nb_points() {
                let a_p1 = a_l.value(a_pt_idx - 1).value();
                let a_p2 = a_l.value(a_pt_idx).value();

                let a_p1p = a_p3d - a_p1;
                let a_p1p2 = a_p2 - a_p1;

                let a_sq12 = a_p1p2.length_squared();

                if a_sq12 < f64::EPSILON {
                    continue;
                }

                let a_dp = a_p1p.dot(a_p1p2);

                let mut a_sq_d = f64::MAX;
                if a_dp < 0.0 {
                    continue;
                } else if a_dp > a_sq12 {
                    continue;
                } else {
                    a_sq_d = a_p1p.cross(a_p1p2).length_squared() / a_sq12;
                }

                if a_sq_d < a_min_sq_dist {
                    a_min_sq_dist = a_sq_d;

                    let a_l1 = a_dp / a_sq12;
                    let a_l2 = 1.0 - a_l1;

                    if a_l1 < eps || a_l2 < eps {
                        return true;
                    }

                    let (a_u1, a_v1) = a_l.value(a_pt_idx - 1).parameters_on_surface(self.reversed);
                    let (a_u2, a_v2) = a_l.value(a_pt_idx).parameters_on_surface(self.reversed);

                    a_umin = a_l1 * a_u2 + a_l2 * a_u1;
                    a_vmin = a_l1 * a_v2 + a_l2 * a_v1;
                }
            }

            if a_min_sq_dist > rcad_kernel::precision::INFINITE_VALUE {
                continue;
            }

            let a_vec_prms = [a_umin, a_vmin];
            solver.perform(func, a_vec_prms, *inf_bounds, *sup_bounds);
            if !solver.is_done() {
                continue;
            }

            let a_vec_prms2 = solver.root();

            let pa = func.p_surface().point_at(a_umin, a_vmin);
            let pb = func.p_surface().point_at(a_vec_prms2[0], a_vec_prms2[1]);
            let a_sq_d1 = pb.distance_squared(a_p3d);
            let a_sq_d2 = pa.distance_squared(pb);

            if a_sq_d1 < 4.0 * a_sq_d2 {
                return true;
            }
        }
        false
    }
}

/// OCCT IntWalk_IWalking::MakeWalkingPoint (gxx L2918-2951) — free function.
fn make_walking_point(
    reversed: bool,
    case: i32,
    u: f64,
    v: f64,
    func: &mut SurfFunction,
    psol: &mut PntOn2S,
) {
    if case == 1 || case == 2 {
        psol.set_value(func.point(), reversed, u, v);
    } else if case == 11 || case == 12 {
        let _ = func.values(&[u, v]);
        make_walking_point(reversed, case - 10, u, v, func, psol);
    }
    // OCCT throws Standard_ConstructionError otherwise.
}

/// OCCT TestPassedSolutionWithNegativeState (gxx L1931-2002).
fn test_passed_solution_with_negative_state(
    wd: &[WalkingData],
    u_mult: &[f64],
    v_mult: &[f64],
    prev_up: f64,
    prev_vp: f64,
    nb_multiplicities: &[i32],
    tolerance: &[f64; 2],
    func: &mut SurfFunction,
    uv: &mut [f64; 2],
    irang: &mut i32,
) -> bool {
    let mut arrive = false;
    let tolu = tolerance[0];
    let tolv = tolerance[1];
    for i in 1..wd.len() {
        if wd[i].etat < -11 {
            let utest = wd[i].ustart;
            let vtest = wd[i].vstart;
            let dup = prev_up - utest;
            let dvp = prev_vp - vtest;
            if dup.abs() >= tolu || dvp.abs() >= tolv {
                let uv1m_utest = uv[0] - utest;
                let uv2m_vtest = uv[1] - vtest;
                if ((dup * uv1m_utest + dvp * uv2m_vtest) < 0.0)
                    || (uv1m_utest.abs() < tolu && uv2m_vtest.abs() < tolv)
                {
                    *irang = i as i32;
                    arrive = true;
                    uv[0] = utest;
                    uv[1] = vtest;
                } else if i < nb_multiplicities.len() && nb_multiplicities[i] > 0 {
                    let mut n: usize = 0;
                    for k in 1..i {
                        n += nb_multiplicities[k] as usize;
                    }
                    let mut j = n;
                    while j < n + nb_multiplicities[i] as usize {
                        if j < u_mult.len()
                            && ((prev_up - u_mult[j]) * (uv[0] - u_mult[j])
                                + (prev_vp - v_mult[j]) * (uv[1] - v_mult[j])
                                < 0.0
                                || (uv[0] - u_mult[j]).abs() < tolu && (uv[1] - v_mult[j]).abs() < tolv)
                        {
                            *irang = i as i32;
                            arrive = true;
                            uv[0] = utest;
                            uv[1] = vtest;
                            break;
                        }
                        j += 1;
                    }
                }
                if arrive {
                    let _ = func.values(&[uv[0], uv[1]]);
                    break;
                }
            }
        }
    }
    arrive
}

/// OCCT CutVectorByTolerances (gxx L400-406).
fn cut_vector_by_tolerances(v: &mut DVec2, tolerance: &[f64; 2]) {
    if v.x.abs() < tolerance[0] {
        v.x = 0.0;
    }
    if v.y.abs() < tolerance[1] {
        v.y = 0.0;
    }
}

/// rcad adaptation of Adaptor3d_HSurfaceTool::UResolution / VResolution
/// (the corrected face domain).
fn u_resolution(domain: [f64; 4], tol3d: f64) -> f64 {
    let u_extent = (domain[1] - domain[0]).abs();
    if u_extent.is_finite() && u_extent > 1e-12 {
        tol3d.max(1e-9) / u_extent
    } else {
        rcad_kernel::precision::PCONFUSION
    }
}

fn v_resolution(domain: [f64; 4], tol3d: f64) -> f64 {
    let v_extent = (domain[3] - domain[2]).abs();
    if v_extent.is_finite() && v_extent > 1e-12 {
        tol3d.max(1e-9) / v_extent
    } else {
        rcad_kernel::precision::PCONFUSION
    }
}
