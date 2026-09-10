// OCCT LocOpe_CSIntersector.hxx L36-151 + LocOpe_CSIntersector.cxx L31-588 +
// LocOpe_CSIntersector.lxx L21-43 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_CSIntersector.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_CSIntersector.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_CSIntersector.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. IntCurvesFace_Intersector (TKTopAlgo) — the per-face line/curve
//    intersector with UV classification is not translated yet; rcad carries
//    the reduced curve-surface intersection vehicle
//    topalgo::brep_int_curve_surface::Inter (UV-domain containment, In
//    transition only). It is re-hosted below as IntCurvesFaceIntersector
//    with the OCCT accessor surface (NbPnt/Pnt/WParameter/UParameter/
//    VParameter/Transition); the reduced semantics are a GAP to close when
//    the full IntCurvesFace translation lands.
// 2. OCCT myPoints is a void* to an array of NCollection_Sequence; rcad is
//    Vec<Vec<LocOpePntFace>> ( Destroy() clears it; the OCCT delete[] runs
//    in the destructor — the Rust drop is the destructor equivalent).
// 3. OCCT LocalizeBefore(I, FromInd, Tol, ...) calls the static
//    LocBefore with the int FromInd implicitly converted to double (only
//    the double overload exists — cxx L294-315 against L38-43/L381-438);
//    the conversion is kept explicit in the translation.
// 4. TopExp_Explorer is feat::brep_feat_builder::explorer (same crate).
//
// first consumer: BRepFeat_MakeCylindricalHole (3a) — Perform(BRepFeat) uses
// LocOpe_CSIntersector over the cylinder axis; BRepFeat_Form family (3b).

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_pnt_face::LocOpePntFace;
use rcad_kernel::geom::{ Circle3, Curve3, Line3, Surface3 };
use rcad_kernel::math::gp::Lin;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ Orientation, ShapeType, TShape };
use rcad_kernel::SurfaceEval;

// OCCT IntCurveSurface_TransitionOnCurve — carried by the rcad intersector.
use crate::topalgo::brep_int_curve_surface::inter::TransitionOnCurve;

/// OCCT IntCurvesFace_Intersector re-host (architecture difference #1) —
/// intersects one face with a line or a curve over a parameter interval.
pub(crate) struct IntCurvesFaceIntersector {
    my_face: Shape,                       // OCCT: myFace
    my_uvw_bounds: [f64; 4],              // face UV domain (rcad vehicle input)
    my_points: Vec<InterPointRecord>,     // OCCT: myPntPoints sequence
}

/// One intersection record of the re-hosted intersector.
pub(crate) struct InterPointRecord {
    pub pnt: glam::DVec3,
    pub w: f64,
    pub u: f64,
    pub v: f64,
    pub transition: TransitionOnCurve,
}

impl IntCurvesFaceIntersector {
    /// OCCT IntCurvesFace_Intersector(F, Tol) — loads the face.
    pub(crate) fn new(the_face: &Shape, _the_tol: f64) -> Self {
        let uv = face_uv_domain(the_face);
        IntCurvesFaceIntersector {
            my_face: the_face.clone(),
            my_uvw_bounds: uv.unwrap_or([0.0, 0.0, 0.0, 0.0]),
            my_points: Vec::new(),
        }
    }

    /// OCCT IntCurvesFace_Intersector::Perform(L, PInf, PSup) — line overload.
    pub(crate) fn perform_lin(&mut self, the_lin: &Lin, p_inf: f64, p_sup: f64) {
        let curve = Curve3::Line(Line3 {
            origin: the_lin.pos,
            direction: the_lin.dir,
        });
        self.perform_curve(&curve, p_inf, p_sup);
    }

    /// OCCT IntCurvesFace_Intersector::Perform(HC, PInf, PSup) — curve
    /// overload over the [PInf, PSup] parameter interval.
    pub(crate) fn perform_curve(&mut self, the_curve: &Curve3, p_inf: f64, p_sup: f64) {
        self.my_points.clear();
        // rcad vehicle: topalgo::brep_int_curve_surface::Inter (arch. diff. #1).
        let Some(surface) = face_surface_world(&self.my_face) else {
            return;
        };
        let mut the_int = crate::topalgo::brep_int_curve_surface::inter::Inter::new();
        the_int.load(&self.my_face, rcad_kernel::precision::CONFUSION);
        the_int.init_curve(
            the_curve,
            &surface,
            self.my_uvw_bounds[0],
            self.my_uvw_bounds[1],
            self.my_uvw_bounds[2],
            self.my_uvw_bounds[3],
        );
        while the_int.more() {
            the_int.next();
            let w = the_int.current_w();
            // OCCT Perform restricts the curve parameter to [PInf, PSup].
            if w < p_inf || w > p_sup {
                continue;
            }
            self.my_points.push(InterPointRecord {
                pnt: the_int.current_point(),
                w,
                u: the_int.current_u(),
                v: the_int.current_v(),
                transition: the_int.current_transition(),
            });
        }
    }

    /// OCCT IntCurvesFace_Intersector::IsDone().
    pub(crate) fn is_done(&self) -> bool {
        true
    }

    /// OCCT IntCurvesFace_Intersector::NbPnt().
    pub(crate) fn nb_pnt(&self) -> i32 {
        self.my_points.len() as i32
    }

    /// OCCT IntCurvesFace_Intersector::Pnt(j) — 1-based.
    pub(crate) fn pnt(&self, j: i32) -> glam::DVec3 {
        self.my_points[(j - 1) as usize].pnt
    }

    /// OCCT IntCurvesFace_Intersector::WParameter(j) — 1-based.
    pub(crate) fn w_parameter(&self, j: i32) -> f64 {
        self.my_points[(j - 1) as usize].w
    }

    /// OCCT IntCurvesFace_Intersector::UParameter(j) — 1-based.
    pub(crate) fn u_parameter(&self, j: i32) -> f64 {
        self.my_points[(j - 1) as usize].u
    }

    /// OCCT IntCurvesFace_Intersector::VParameter(j) — 1-based.
    pub(crate) fn v_parameter(&self, j: i32) -> f64 {
        self.my_points[(j - 1) as usize].v
    }

    /// OCCT IntCurvesFace_Intersector::Transition(j) — 1-based.
    pub(crate) fn transition(&self, j: i32) -> TransitionOnCurve {
        self.my_points[(j - 1) as usize].transition
    }
}

/// OCCT BRep_Tool::Surface(F) with location applied (identity on standalone
/// feat shapes — loc_ope_find_edges.rs architecture difference #1).
fn face_surface_world(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// Face UV domain (BRepAdaptor_Surface::FirstUParameter/... vehicle).
fn face_uv_domain(face: &Shape) -> Option<[f64; 4]> {
    match face.data.as_ref() {
        TShape::Face(fd) => match (&fd.uv_domain, &fd.surface) {
            (Some(d), _) => Some(*d),
            (None, Some(s)) => Some(s.default_domain()),
            (None, None) => None,
        },
        _ => None,
    }
}

/// OCCT LocOpe_CSIntersector (LocOpe_CSIntersector.hxx L36-151).
pub struct LocOpeCSIntersector {
    my_done: bool,                        // OCCT: myDone
    my_shape: Shape,                      // OCCT: myShape
    my_points: Vec<Vec<LocOpePntFace>>,   // OCCT: myPoints (arch. diff. #2)
    my_nbelem: i32,                       // OCCT: myNbelem
}

impl Default for LocOpeCSIntersector {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeCSIntersector {
    /// OCCT LocOpe_CSIntersector::LocOpe_CSIntersector() (lxx L21-26).
    pub fn new() -> Self {
        LocOpeCSIntersector {
            my_done: false,
            my_shape: Shape::null(),
            my_points: Vec::new(),
            my_nbelem: 0,
        }
    }

    /// OCCT LocOpe_CSIntersector::LocOpe_CSIntersector(S) (lxx L30-36).
    pub fn with_shape(the_s: &Shape) -> Self {
        LocOpeCSIntersector {
            my_done: false,
            my_shape: the_s.clone(),
            my_points: Vec::new(),
            my_nbelem: 0,
        }
    }

    /// OCCT LocOpe_CSIntersector::Init(S) (cxx L58-65).
    pub fn init(&mut self, the_s: &Shape) {
        self.my_done = false;
        self.my_shape = the_s.clone();
        self.my_points = Vec::new();
        self.my_nbelem = 0;
    }

    /// OCCT LocOpe_CSIntersector::Perform(Slin) (cxx L69-99) — line sequence
    /// overload.
    pub fn perform_lin(&mut self, the_slin: &[Lin]) {
        if self.my_shape.is_null() || the_slin.is_empty() {
            // OCCT cxx L71-74: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        }
        self.my_done = false;

        // OCCT cxx L77-80: myNbelem = Slin.Length(); new sequences array.
        self.my_nbelem = the_slin.len() as i32;
        self.my_points = vec![Vec::new(); self.my_nbelem as usize];

        // OCCT cxx L82-83: binf = RealFirst(); bsup = RealLast().
        let binf = f64::MIN; // OCCT RealFirst() = -DBL_MAX
        let bsup = f64::MAX; // OCCT RealLast() = DBL_MAX
        // OCCT cxx L84-97: per face, per line.
        for theface in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            let mut the_int = IntCurvesFaceIntersector::new(&theface, rcad_kernel::precision::PCONFUSION);
            for (i, lin) in the_slin.iter().enumerate() {
                the_int.perform_lin(lin, binf, bsup);
                if the_int.is_done() {
                    add_points(&the_int, &mut self.my_points[i], &theface);
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT LocOpe_CSIntersector::Perform(Scir) (cxx L103-137) — circle
    /// sequence overload.
    pub fn perform_cir(&mut self, the_scir: &[Circle3]) {
        if self.my_shape.is_null() || the_scir.is_empty() {
            // OCCT cxx L105-108: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        }
        self.my_done = false;

        // OCCT cxx L111-114.
        self.my_nbelem = the_scir.len() as i32;
        self.my_points = vec![Vec::new(); self.my_nbelem as usize];

        // OCCT cxx L118-119: binf = 0.; bsup = 2.*M_PI.
        let binf = 0.0;
        let bsup = 2.0 * std::f64::consts::PI;

        // OCCT cxx L116-135: per face, per circle (HC->Load(Geom_Circle)).
        for theface in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            let mut the_int = IntCurvesFaceIntersector::new(&theface, 0.0);
            for (i, cir) in the_scir.iter().enumerate() {
                let curve = Curve3::Circle(*cir);
                the_int.perform_curve(&curve, binf, bsup);
                if the_int.is_done() {
                    add_points(&the_int, &mut self.my_points[i], &theface);
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT LocOpe_CSIntersector::Perform(Scur) (cxx L141-177) — generic
    /// curve sequence overload (None carries the OCCT null handle).
    pub fn perform_cur(&mut self, the_scur: &[Option<Curve3>]) {
        if self.my_shape.is_null() || the_scur.is_empty() {
            // OCCT cxx L143-146: throw Standard_ConstructionError().
            panic!("Standard_ConstructionError");
        }
        self.my_done = false;

        // OCCT cxx L149-152.
        self.my_nbelem = the_scur.len() as i32;
        self.my_points = vec![Vec::new(); self.my_nbelem as usize];

        // OCCT cxx L154-175: per face, per curve (null handle -> continue).
        for theface in explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape) {
            let mut the_int = IntCurvesFaceIntersector::new(&theface, 0.0);
            for (i, cur) in the_scur.iter().enumerate() {
                let Some(cur) = cur else {
                    continue;
                };
                // OCCT cxx L167-168: binf/bsup from the curve itself.
                use rcad_kernel::geom::CurveEval;
                let binf = cur.default_domain()[0];
                let bsup = cur.default_domain()[1];
                the_int.perform_curve(cur, binf, bsup);
                if the_int.is_done() {
                    add_points(&the_int, &mut self.my_points[i], &theface);
                }
            }
        }
        self.my_done = true;
    }

    /// OCCT LocOpe_CSIntersector::IsDone() (lxx L40-43).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT LocOpe_CSIntersector::NbPoints(I) (cxx L181-192).
    pub fn nb_points(&self, the_i: i32) -> i32 {
        if !self.my_done {
            // OCCT cxx L183-186: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if the_i <= 0 || the_i > self.my_nbelem {
            // OCCT cxx L187-190: throw Standard_OutOfRange().
            panic!("Standard_OutOfRange");
        }
        self.my_points[(the_i - 1) as usize].len() as i32
    }

    /// OCCT LocOpe_CSIntersector::Point(I, Index) (cxx L196-207).
    pub fn point(&self, the_i: i32, the_index: i32) -> LocOpePntFace {
        if !self.my_done {
            // OCCT cxx L198-201: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if the_i <= 0 || the_i > self.my_nbelem {
            // OCCT cxx L202-205: throw Standard_OutOfRange().
            panic!("Standard_OutOfRange");
        }
        self.my_points[(the_i - 1) as usize][(the_index - 1) as usize].clone()
    }

    /// OCCT LocOpe_CSIntersector::Destroy() (cxx L211-215) — the OCCT
    /// destructor calls it; the Rust drop runs the same release.
    pub fn destroy(&mut self) {
        self.my_points = Vec::new();
    }

    /// OCCT LocOpe_CSIntersector::LocalizeAfter(I, From, Tol, Or, IndFrom,
    /// IndTo) (cxx L219-240) — the double From overload.
    pub fn localize_after(
        &self,
        the_i: i32,
        the_from: f64,
        the_tol: f64,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L226-229: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if the_i <= 0 || the_i > self.my_nbelem {
            // OCCT cxx L230-233: throw Standard_OutOfRange().
            panic!("Standard_OutOfRange");
        }
        loc_after_param(
            &self.my_points[(the_i - 1) as usize],
            the_from,
            the_tol,
            the_or,
            the_ind_from,
            the_ind_to,
        )
    }

    /// OCCT LocOpe_CSIntersector::LocalizeBefore(I, From, Tol, Or, IndFrom,
    /// IndTo) (cxx L244-265) — the double From overload.
    pub fn localize_before(
        &self,
        the_i: i32,
        the_from: f64,
        the_tol: f64,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L251-254: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if the_i <= 0 || the_i > self.my_nbelem {
            // OCCT cxx L255-258: throw Standard_OutOfRange().
            panic!("Standard_OutOfRange");
        }
        loc_before_param(
            &self.my_points[(the_i - 1) as usize],
            the_from,
            the_tol,
            the_or,
            the_ind_from,
            the_ind_to,
        )
    }

    /// OCCT LocOpe_CSIntersector::LocalizeAfter(I, FromInd, Tol, Or, IndFrom,
    /// IndTo) (cxx L269-290) — the int FromInd overload.
    pub fn localize_after_index(
        &self,
        the_i: i32,
        the_from_ind: i32,
        the_tol: f64,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L276-279: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if the_i <= 0 || the_i > self.my_nbelem {
            // OCCT cxx L280-283: throw Standard_OutOfRange().
            panic!("Standard_OutOfRange");
        }
        loc_after_ind(
            &self.my_points[(the_i - 1) as usize],
            the_from_ind,
            the_tol,
            the_or,
            the_ind_from,
            the_ind_to,
        )
    }

    /// OCCT LocOpe_CSIntersector::LocalizeBefore(I, FromInd, Tol, Or,
    /// IndFrom, IndTo) (cxx L294-315) — the int FromInd overload; OCCT
    /// resolves the static call through the implicit int->double conversion
    /// to LocBefore(const double, ...) (architecture difference #3).
    pub fn localize_before_index(
        &self,
        the_i: i32,
        the_from_ind: i32,
        the_tol: f64,
        the_or: &mut Orientation,
        the_ind_from: &mut i32,
        the_ind_to: &mut i32,
    ) -> bool {
        if !self.my_done {
            // OCCT cxx L301-304: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        if the_i <= 0 || the_i > self.my_nbelem {
            // OCCT cxx L305-308: throw Standard_OutOfRange().
            panic!("Standard_OutOfRange");
        }
        loc_before_param(
            &self.my_points[(the_i - 1) as usize],
            the_from_ind as f64,
            the_tol,
            the_or,
            the_ind_from,
            the_ind_to,
        )
    }
}

/// OCCT static LocAfter(const NCollection_Sequence<LocOpe_PntFace>&,
/// const double From, const double Tol, ...) (cxx L319-377).
fn loc_after_param(
    the_spt: &[LocOpePntFace],
    the_from: f64,
    the_tol: f64,
    the_or: &mut Orientation,
    the_ind_from: &mut i32,
    the_ind_to: &mut i32,
) -> bool {
    // OCCT cxx L327: FMEPS = From - Tol; nbpoints = Spt.Length().
    let mut param;
    let fmeps = the_from - the_tol;
    let nbpoints = the_spt.len() as i32;
    // OCCT cxx L329-335.
    let mut ifirst = 1;
    while ifirst <= nbpoints {
        if the_spt[(ifirst - 1) as usize].parameter() >= fmeps {
            break;
        }
        ifirst += 1;
    }
    let mut ret_val = false;
    if ifirst <= nbpoints {
        // OCCT cxx L337-341.
        let mut i = ifirst;
        *the_ind_from = ifirst;
        let mut found = false;
        while !found {
            // OCCT cxx L344-346.
            *the_or = the_spt[(i - 1) as usize].orientation();
            param = the_spt[(i - 1) as usize].parameter();
            i += 1;
            // OCCT cxx L347-361.
            while i <= nbpoints {
                if the_spt[(i - 1) as usize].parameter() - param <= the_tol {
                    if *the_or != Orientation::External
                        && *the_or != the_spt[(i - 1) as usize].orientation()
                    {
                        *the_or = Orientation::External;
                    }
                    i += 1;
                } else {
                    break;
                }
            }
            // OCCT cxx L362-373.
            if *the_or == Orientation::External {
                found = i > nbpoints;
                *the_ind_from = i;
            } else {
                // on a une intersection franche
                *the_ind_to = i - 1;
                found = true;
                ret_val = true;
            }
        }
    }

    ret_val
}

/// OCCT static LocBefore(const NCollection_Sequence<LocOpe_PntFace>&,
/// const double From, const double Tol, ...) (cxx L381-438).
fn loc_before_param(
    the_spt: &[LocOpePntFace],
    the_from: f64,
    the_tol: f64,
    the_or: &mut Orientation,
    the_ind_from: &mut i32,
    the_ind_to: &mut i32,
) -> bool {
    // OCCT cxx L388: FPEPS = From + Tol.
    let mut param;
    let fpeps = the_from + the_tol;
    let nbpoints = the_spt.len() as i32;
    // OCCT cxx L390-396.
    let mut ifirst = nbpoints;
    while ifirst >= 1 {
        if the_spt[(ifirst - 1) as usize].parameter() <= fpeps {
            break;
        }
        ifirst -= 1;
    }
    let mut ret_val = false;
    if ifirst >= 1 {
        // OCCT cxx L398-402.
        let mut i = ifirst;
        *the_ind_to = ifirst;
        let mut found = false;
        while !found {
            // OCCT cxx L404-407.
            *the_or = the_spt[(i - 1) as usize].orientation();
            param = the_spt[(i - 1) as usize].parameter();
            i -= 1;
            // OCCT cxx L408-422.
            while i >= 1 {
                if param - the_spt[(i - 1) as usize].parameter() <= the_tol {
                    if *the_or != Orientation::External
                        && *the_or != the_spt[(i - 1) as usize].orientation()
                    {
                        *the_or = Orientation::External;
                    }
                    i -= 1;
                } else {
                    break;
                }
            }
            // OCCT cxx L423-434.
            if *the_or == Orientation::External {
                found = i < 1;
                *the_ind_to = i;
            } else {
                // on a une intersection franche
                *the_ind_from = i + 1;
                found = true;
                ret_val = true;
            }
        }
    }

    ret_val
}

/// OCCT static LocAfter(const NCollection_Sequence<LocOpe_PntFace>&,
/// const int FromInd, const double Tol, ...) (cxx L442-513).
fn loc_after_ind(
    the_spt: &[LocOpePntFace],
    the_from_ind: i32,
    the_tol: f64,
    the_or: &mut Orientation,
    the_ind_from: &mut i32,
    the_ind_to: &mut i32,
) -> bool {
    let nbpoints = the_spt.len() as i32;
    // OCCT cxx L449-453.
    if the_from_ind >= nbpoints {
        return false;
    }

    let mut param;
    // OCCT cxx L457-471: FMEPS is read only on the FromInd >= 1 path.
    let ifirst;
    if the_from_ind >= 1 {
        let fmeps = the_spt[(the_from_ind - 1) as usize].parameter() - the_tol;
        let mut it = the_from_ind + 1;
        loop {
            if it > nbpoints {
                break;
            }
            if the_spt[(it - 1) as usize].parameter() >= fmeps {
                break;
            }
            it += 1;
        }
        ifirst = it;
    } else {
        ifirst = 1;
    }

    let mut ret_val = false;
    if ifirst <= nbpoints {
        // OCCT cxx L474-478.
        let mut i = ifirst;
        *the_ind_from = ifirst;
        let mut found = false;
        while !found {
            // OCCT cxx L481-483.
            *the_or = the_spt[(i - 1) as usize].orientation();
            param = the_spt[(i - 1) as usize].parameter();
            i += 1;
            // OCCT cxx L484-498.
            while i <= nbpoints {
                if the_spt[(i - 1) as usize].parameter() - param <= the_tol {
                    if *the_or != Orientation::External
                        && *the_or != the_spt[(i - 1) as usize].orientation()
                    {
                        *the_or = Orientation::External;
                    }
                    i += 1;
                } else {
                    break;
                }
            }
            // OCCT cxx L499-510.
            if *the_or == Orientation::External {
                found = i > nbpoints;
                *the_ind_from = i;
            } else {
                // on a une intersection franche
                *the_ind_to = i - 1;
                found = true;
                ret_val = true;
            }
        }
    }
    ret_val
}

/// OCCT static AddPoints(IntCurvesFace_Intersector&, sequence, face)
/// (cxx L517-588).
fn add_points(
    the_int: &IntCurvesFaceIntersector,
    the_seq: &mut Vec<LocOpePntFace>,
    the_face: &Shape,
) {
    // OCCT cxx L521-523.
    let mut nbpoints = the_seq.len() as i32;
    let newpnt = the_int.nb_pnt();
    // OCCT cxx L524-587.
    for j in 1..=newpnt {
        let thept = the_int.pnt(j);
        let param = the_int.w_parameter(j);
        let paramu = the_int.u_parameter(j);
        let paramv = the_int.v_parameter(j);

        // OCCT cxx L531-567: theor from the transition (the JAG 13.09.96
        // commented orface branches are not translated).
        let theor = match the_int.transition(j) {
            TransitionOnCurve::In => Orientation::Forward,
            TransitionOnCurve::Out => Orientation::Reversed,
            TransitionOnCurve::Tangent => Orientation::External,
        };
        let newpt = LocOpePntFace::new_full(thept, the_face, theor, param, paramu, paramv);
        // OCCT cxx L569-577: insertion by increasing parameter.
        let mut k = 1;
        while k <= nbpoints {
            if the_seq[(k - 1) as usize].parameter() - param > 0.0 {
                break;
            }
            k += 1;
        }
        if k <= nbpoints {
            the_seq.insert((k - 1) as usize, newpt);
        } else {
            the_seq.push(newpt);
        }
        nbpoints += 1;
    }
}
