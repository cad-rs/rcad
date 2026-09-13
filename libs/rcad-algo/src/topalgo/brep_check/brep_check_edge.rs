//! OCCT BRepCheck_Edge (TKTopAlgo/BRepCheck).
//!
//! Source: `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Edge.cxx`
//! (L61-840) and `BRepCheck_Edge.hxx` (L28-72).
//!
//! `Tolerance()` (Edge.cxx L598-707) is not ported: it is not part of the
//! BRepCheck_Analyzer pass structure (remaining work).

use rcad_kernel::base::geom_proj_lib::project_on_plane::curve_on_plane;
use rcad_kernel::geom::{transform_curve, transform_surface, Curve2dEval, CurveEval, Surface3};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, ShapeType};

use crate::topalgo::brep_lib_validate_edge::{
    Adaptor3dCurveOnSurface, Geom2dAdaptorCurve, GeomAdaptorCurve, GeomAdaptorSurface,
    BRepLibValidateEdge,
};

use crate::brep_algo::tool::brep_tool_tolerance;
use super::brep_check_result::{
    brep_check_add, edge_curve_reps, explorer, location_matrix,
    EdgeCurveRep, BRepCheckResultBase, BRepCheckStatus,
};

/// OCCT Edge.cxx L61: `static const int NCONTROL = 23;` (used only by the
/// unported `Tolerance()`).
#[allow(dead_code)]
pub const NCONTROL: i32 = 23;

/// OCCT `Handle(Adaptor3d_Curve) myHCurve` — either a GeomAdaptor_Curve (the
/// 3D curve reference) or an Adaptor3d_CurveOnSurface (the pcurve reference).
#[derive(Debug, Clone)]
pub enum HCurveAdaptor {
    /// OCCT GeomAdaptor_Curve.
    Curve3d(GeomAdaptorCurve),
    /// OCCT Adaptor3d_CurveOnSurface.
    OnSurface(Adaptor3dCurveOnSurface),
}

impl HCurveAdaptor {
    /// OCCT Adaptor3d_Curve::FirstParameter.
    pub fn first_parameter(&self) -> f64 {
        match self {
            HCurveAdaptor::Curve3d(c) => c.first_parameter(),
            HCurveAdaptor::OnSurface(c) => c.first_parameter(),
        }
    }

    /// OCCT Adaptor3d_Curve::LastParameter.
    pub fn last_parameter(&self) -> f64 {
        match self {
            HCurveAdaptor::Curve3d(c) => c.last_parameter(),
            HCurveAdaptor::OnSurface(c) => c.last_parameter(),
        }
    }

    /// OCCT Adaptor3d_Curve::Value(U).
    pub fn value(&self, the_u: f64) -> glam::DVec3 {
        match self {
            HCurveAdaptor::Curve3d(c) => c.value(the_u),
            HCurveAdaptor::OnSurface(c) => c.value(the_u),
        }
    }
}

/// The three BRepLib_ValidateEdge call sites in InContext(FACE).
#[derive(Debug, Clone, Copy)]
enum RunKind {
    /// Edge.cxx L411-432 (the representation pcurve).
    First { closed: bool },
    /// Edge.cxx L438-453 (the seam pcurve2).
    ClosedSecond,
    /// Edge.cxx L504-512 (the on-the-fly plane projection).
    Projection,
}

/// OCCT BRepCheck_Edge (Edge.hxx L28-72).
#[derive(Debug)]
pub struct BRepCheckEdge {
    /// OCCT protected base (Result.hxx L82-90).
    pub base: BRepCheckResultBase,
    /// OCCT myCref — the reference curve representation.
    pub my_cref: Option<EdgeCurveRep>,
    /// OCCT myHCurve.
    pub my_h_curve: Option<HCurveAdaptor>,
    /// OCCT myGctrl.
    pub my_gctrl: bool,
    /// OCCT myIsExactMethod.
    pub my_is_exact_method: bool,
}

impl BRepCheckEdge {
    /// OCCT BRepCheck_Edge::BRepCheck_Edge(const TopoDS_Edge& E)
    /// (Edge.cxx L65-70).
    pub fn new(brep: &BRep, e: &Shape) -> Self {
        let mut r = BRepCheckEdge {
            base: BRepCheckResultBase::new(),
            my_cref: None,
            my_h_curve: None,
            my_gctrl: true,
            my_is_exact_method: false,
        };
        r.base.init(e);
        r.minimum(brep);
        // OCCT L68-69.
        r.my_gctrl = true;
        r.my_is_exact_method = false;
        r
    }

    /// OCCT BRepCheck_Edge::Minimum (Edge.cxx L74-259).
    pub fn minimum(&mut self, brep: &BRep) {
        if !self.base.my_min {
            // OCCT L78-80: bind a fresh list for myShape.
            self.base.my_map.bound(&self.base.my_shape);
            let my_shape = self.base.my_shape.clone();
            // OCCT L81: myCref.Nullify().
            self.my_cref = None;

            // OCCT L84: TE = the edge TShape.
            let ed_data = my_shape.as_edge();

            // OCCT L85-119: existence and uniqueness of a 3D representation.
            let reps = edge_curve_reps(brep, &my_shape);
            let mut exist = false;
            let mut unique = true;
            for cr in &reps {
                if cr.is_curve3d() {
                    if !exist {
                        exist = true;
                    } else {
                        unique = false;
                    }
                    let has_value = matches!(cr, EdgeCurveRep::Curve3D { .. });
                    if self.my_cref.is_none() && has_value {
                        self.my_cref = Some(cr.clone());
                    }
                }
            }

            // OCCT L91-97.
            let degenerated = ed_data.map(|e| e.degenerated).unwrap_or(false);
            let same_parameter = ed_data.map(|e| e.same_parameter).unwrap_or(false);
            let same_range = ed_data.map(|e| e.same_range).unwrap_or(false);
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            if !same_range && same_parameter {
                brep_check_add(lst, BRepCheckStatus::InvalidSameParameterFlag);
            }

            // OCCT L121-129.
            if !exist {
                brep_check_add(lst, BRepCheckStatus::No3DCurve);
                // myCref est nulle
            } else if !unique {
                brep_check_add(lst, BRepCheckStatus::Multiple3DCurve);
            }

            // OCCT L131-148: search for a 3D reference; if none exists, take
            // the first CurveOnSurf.
            if self.my_cref.is_none() && !degenerated {
                for cr in &reps {
                    if cr.is_curve_on_surface() {
                        self.my_cref = Some(cr.clone());
                        break;
                    }
                }
            } else if self.my_cref.is_some() && degenerated {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Minimum: myShape must be bound");
                brep_check_add(lst, BRepCheckStatus::InvalidDegeneratedFlag);
            }

            // OCCT L150-252: the range checks and the reference adaptor.
            if let Some(cref) = self.my_cref.clone() {
                let first;
                let last;
                match &cref {
                    EdgeCurveRep::Curve3D { .. } => {
                        let r = ed_data.map(|e| e.range).unwrap_or([0.0, 0.0]);
                        first = r[0];
                        last = r[1];
                    }
                    EdgeCurveRep::CurveOnSurface { range, .. }
                    | EdgeCurveRep::CurveOnClosedSurface { range, .. } => {
                        first = range[0];
                        last = range[1];
                    }
                    _ => {
                        // Curve3DUnresolved / Regularity: no usable range —
                        // GAP (no status set).
                        first = 0.0;
                        last = 0.0;
                    }
                }
                // OCCT L153: constexpr double eps = Precision::PConfusion().
                let eps = rcad_kernel::precision::PCONFUSION;
                if last <= first {
                    // OCCT L156-160.
                    self.my_cref = None;
                    let lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Minimum: myShape must be bound");
                    brep_check_add(lst, BRepCheckStatus::InvalidRange);
                } else {
                    match &cref {
                        EdgeCurveRep::Curve3D { curve, .. } => {
                            // OCCT L168: L = myShape.Location() * myCref->Location().
                            let l_mat = location_matrix(brep, my_shape.location, cref_loc(&cref));
                            // OCCT L169-170: C3d = Curve3D()->Transformed(L).
                            let c3d = transform_curve(curve, &l_mat);
                            // OCCT L171-176: IsPeriodic / aPeriod.
                            let mut is_periodic = CurveEval::is_periodic(&c3d);
                            let mut a_period = f64::MAX;
                            if is_periodic {
                                a_period = period_of_curve3d(&c3d);
                            }
                            // OCCT L177: f = C3d->FirstParameter(), l = C3d->LastParameter().
                            let dom = CurveEval::default_domain(&c3d);
                            let mut f = dom[0];
                            let mut l = dom[1];
                            // OCCT L178-189: the Geom_TrimmedCurve basis unwrap.
                            if let rcad_kernel::geom::Curve3::Trimmed(tc) = &c3d {
                                let a_c: &rcad_kernel::geom::Curve3 = &tc.curve;
                                let bdom = CurveEval::default_domain(a_c);
                                f = bdom[0];
                                l = bdom[1];
                                is_periodic = CurveEval::is_periodic(a_c);
                                if is_periodic {
                                    a_period = period_of_curve3d(a_c);
                                }
                            }
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("Minimum: myShape must be bound");
                            if is_periodic && (last - first > a_period + eps) {
                                // OCCT L190-194.
                                self.my_cref = None;
                                brep_check_add(lst, BRepCheckStatus::InvalidRange);
                            } else if !is_periodic && (first < f - eps || last > l + eps) {
                                // OCCT L195-199.
                                self.my_cref = None;
                                brep_check_add(lst, BRepCheckStatus::InvalidRange);
                            } else {
                                // OCCT L200-206: the GeomAdaptor_Curve over the
                                // transformed parameters (identity for the rigid
                                // locations rcad materialises — the kernel
                                // `transformed_parameter` default).
                                let gac = GeomAdaptorCurve::new(c3d, first, last);
                                self.my_h_curve = Some(HCurveAdaptor::Curve3d(gac));
                            }
                        }
                        EdgeCurveRep::CurveOnSurface { pcurve, .. }
                        | EdgeCurveRep::CurveOnClosedSurface { pcurve1: pcurve, .. } => {
                            // OCCT L209-249: the curve-on-surface branch.
                            // OCCT L210-212: Sref = myCref->Surface() transformed
                            // by myCref->Location(). rcad resolves the surface
                            // through the owning face.
                            let sref = self.resolve_cref_surface(brep, &cref);
                            let pc = pcurve.clone();
                            // OCCT L214-219: IsPeriodic / aPeriod.
                            let mut is_periodic = Curve2dEval::is_periodic(&pc);
                            let mut a_period = f64::MAX;
                            if is_periodic {
                                a_period = period_of_curve2d(&pc);
                            }
                            // OCCT L220.
                            let dom = Curve2dEval::default_domain(&pc);
                            let mut f = dom[0];
                            let mut l = dom[1];
                            // OCCT L221-232: the Geom2d_TrimmedCurve unwrap.
                            if let rcad_kernel::geom::Curve2d::Trimmed(tc) = &pc {
                                let a_c: &rcad_kernel::geom::Curve2d = &tc.curve;
                                let bdom = Curve2dEval::default_domain(a_c);
                                f = bdom[0];
                                l = bdom[1];
                                is_periodic = Curve2dEval::is_periodic(a_c);
                                if is_periodic {
                                    a_period = period_of_curve2d(a_c);
                                }
                            }
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("Minimum: myShape must be bound");
                            if is_periodic && (last - first > a_period + eps) {
                                // OCCT L233-237.
                                self.my_cref = None;
                                brep_check_add(lst, BRepCheckStatus::InvalidRange);
                            } else if !is_periodic && (first < f - eps || last > l + eps) {
                                // OCCT L238-242.
                                self.my_cref = None;
                                brep_check_add(lst, BRepCheckStatus::InvalidRange);
                            } else if let Some(sref) = sref {
                                // OCCT L243-249: the Adaptor3d_CurveOnSurface.
                                let gahs = GeomAdaptorSurface::new(sref);
                                let ghpc = Geom2dAdaptorCurve::new(pc, first, last);
                                let acs = Adaptor3dCurveOnSurface::new(ghpc, gahs);
                                self.my_h_curve = Some(HCurveAdaptor::OnSurface(acs));
                            } else {
                                // The owning face of the pcurve is not in this
                                // pool — GAP (no reference adaptor built).
                            }
                        }
                        _ => {
                            // Curve3DUnresolved / Regularity — GAP (no status).
                        }
                    }
                }
            }
            // OCCT L253-257.
            let lst = self
                .base
                .my_map
                .find_mut(&self.base.my_shape)
                .expect("Minimum: myShape must be bound");
            if lst.is_empty() {
                lst.push(BRepCheckStatus::NoError);
            }
            // OCCT L257.
            self.base.my_min = true;
        }
    }

    /// OCCT BRepCheck_Edge::Tolerance() (Edge.cxx L598-707).
    ///
    /// Collects every representation of the edge as a 3D-samplable adaptor
    /// (the 3D curve, the pcurves, the seam pcurves), samples NCONTROL
    /// parameters along [First, Last] and returns the maximal representation
    /// gap with the OCCT 5% margin. An edge carrying at most one
    /// representation yields `Precision::Confusion()`; an infinite sampled
    /// coordinate yields `Precision::Infinite()`.
    ///
    /// The OCCT location composition is preserved: the 3D-curve and the first
    /// pcurve use `myShape.Location() * cr->Location()`, while the seam
    /// pcurve uses `cr->Location()` alone (Edge.cxx L669).
    pub fn tolerance(&self, brep: &BRep) -> f64 {
        use rcad_kernel::core::precision::{is_infinite_value, CONFUSION, INFINITE_VALUE};
        let my_shape = self.base.my_shape.clone();
        // L601-604: nbRep = TE->Curves().Extent(); <= 1 -> Confusion().
        let reps = edge_curve_reps(brep, &my_shape);
        let mut nb_rep = reps.len() as i32;
        if nb_rep <= 1 {
            return CONFUSION;
        }
        // L606-614: First / Last (myHCurve when loaded, else the edge range).
        let (first, last) = match &self.my_h_curve {
            Some(h) => (h.first_parameter(), h.last_parameter()),
            None => match my_shape.as_edge() {
                Some(ed) => (ed.range[0], ed.range[1]),
                None => (0.0, 1.0),
            },
        };
        // L616-618: NCollection_Array1(1, nbRep*2) — 1-based slots.
        let mut the_rep: Vec<Option<HCurveAdaptor>> = vec![None; (nb_rep * 2) as usize + 1];
        let degenerated = my_shape.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
        let mut i_rep: i32 = 1;
        // L620-683: the representation walk.
        for cr in &reps {
            match cr {
                EdgeCurveRep::Curve3D { curve, location } if !degenerated => {
                    // L631-642: Loc = myShape.Location() * cr->Location().
                    let loc = location_matrix(brep, my_shape.location, *location);
                    let c3d = transform_curve(curve, &loc);
                    let gac = GeomAdaptorCurve::new(c3d, first, last);
                    // The OCCT slot dance (L636-641): a later 3D curve moves
                    // the previous slot-1 adaptor into its own slot and takes
                    // slot 1 (the reference the sampling loop reads).
                    let mut it = i_rep;
                    if i_rep > 1 {
                        the_rep[i_rep as usize] = the_rep[1].clone();
                        it = 1;
                    }
                    the_rep[it as usize] = Some(HCurveAdaptor::Curve3d(gac));
                    i_rep += 1;
                }
                EdgeCurveRep::CurveOnSurface { pcurve, location, .. }
                | EdgeCurveRep::CurveOnClosedSurface { pcurve1: pcurve, location, .. } => {
                    let Some(sref) = self.resolve_cref_surface(brep, cr) else {
                        // The owning face is not in this pool — the adaptor
                        // cannot be built (the OCCT else accounting).
                        nb_rep -= 1;
                        continue;
                    };
                    // L643-661: Sref transformed by myShape.Location() * cr->Location().
                    let loc = location_matrix(brep, my_shape.location, *location);
                    let sref = transform_surface(&sref, &loc);
                    let ghpc = Geom2dAdaptorCurve::new(pcurve.clone(), first, last);
                    let gahs = GeomAdaptorSurface::new(sref);
                    the_rep[i_rep as usize] =
                        Some(HCurveAdaptor::OnSurface(Adaptor3dCurveOnSurface::new(ghpc, gahs)));
                    i_rep += 1;
                    // L662-672: the seam pcurve2 — surface transformed by
                    // cr->Location() ALONE.
                    if let EdgeCurveRep::CurveOnClosedSurface { pcurve2, .. } = cr {
                        let Some(sref2) = self.resolve_cref_surface(brep, cr) else {
                            nb_rep -= 1;
                            continue;
                        };
                        let loc2 = location_matrix(brep, 0, *location);
                        let sref2 = transform_surface(&sref2, &loc2);
                        let ghpc2 = Geom2dAdaptorCurve::new(pcurve2.clone(), first, last);
                        let gahs2 = GeomAdaptorSurface::new(sref2);
                        the_rep[i_rep as usize] =
                            Some(HCurveAdaptor::OnSurface(Adaptor3dCurveOnSurface::new(ghpc2, gahs2)));
                        i_rep += 1;
                        nb_rep += 1;
                    }
                }
                _ => {
                    // L680-682: neither a usable 3D curve nor a curve on
                    // surface (degenerated 3D curve, regularity, unresolved).
                    nb_rep -= 1;
                }
            }
        }
        // L674-702: the NCONTROL samples; the reference is slot 1.
        let mut tol_cal = 0.0f64;
        for i in 0..NCONTROL {
            let prm = ((NCONTROL - 1 - i) as f64 * first + i as f64 * last)
                / (NCONTROL - 1) as f64;
            let Some(reference) = the_rep[1].as_ref() else {
                return INFINITE_VALUE;
            };
            let center = reference.value(prm);
            if is_infinite_value(center.x)
                || is_infinite_value(center.y)
                || is_infinite_value(center.z)
            {
                return INFINITE_VALUE;
            }
            for i_rep2 in 2..=nb_rep {
                let Some(other) = the_rep[i_rep2 as usize].as_ref() else {
                    continue;
                };
                let other_p = other.value(prm);
                if is_infinite_value(other_p.x)
                    || is_infinite_value(other_p.y)
                    || is_infinite_value(other_p.z)
                {
                    return INFINITE_VALUE;
                }
                let dist2 = center.distance_squared(other_p);
                if dist2 > tol_cal {
                    tol_cal = dist2;
                }
            }
        }
        // L705-706: the 5% margin.
        tol_cal.sqrt() * 1.05
    }

    /// The surface value of the reference pcurve representation (OCCT
    /// `myCref->Surface()`): resolved through the owning face TShape.
    fn resolve_cref_surface(&self, brep: &BRep, cref: &EdgeCurveRep) -> Option<Surface3> {
        let face_key = match cref {
            EdgeCurveRep::CurveOnSurface { face_key, .. }
            | EdgeCurveRep::CurveOnClosedSurface { face_key, .. } => *face_key,
            _ => return None,
        };
        let ts = brep
            .tshapes
            .iter()
            .find(|ts| std::sync::Arc::as_ptr(ts) as u64 == face_key.0)?;
        match ts.as_ref() {
            rcad_kernel::topods::TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        }
    }

    /// OCCT BRepCheck_Edge::InContext (Edge.cxx L263-555).
    pub fn in_context(&mut self, brep: &BRep, s: &Shape) {
        // OCCT L265-280: bound check under the (parallel) lock.
        if self.base.my_map.is_bound(s) {
            return;
        }
        self.base.my_map.bind(s.clone(), Vec::new());

        let my_shape = self.base.my_shape.clone();
        // OCCT L284-285: TE + Tol.
        let ed_data = my_shape.as_edge();
        let tol = brep_tool_tolerance(&my_shape);

        let styp = s.shape_type();
        // OCCT L289-301: look for myShape among the sub-edges of S.
        {
            let mut found = false;
            let mut more = false;
            for cur in explorer(brep, s, ShapeType::Edge) {
                more = true;
                if cur.is_same(&my_shape) {
                    found = true;
                    break;
                }
            }
            if !more || !found {
                let lst = self
                    .base
                    .my_map
                    .find_mut(s)
                    .expect("InContext: the context must be bound");
                brep_check_add(lst, BRepCheckStatus::SubshapeNotInShape);
                return;
            }
        }

        match styp {
            // OCCT L305-307: WIRE — nothing.
            ShapeType::Wire => {}
            ShapeType::Face => {
                // OCCT L309-517.
                if self.my_cref.is_some() {
                    // OCCT L313-314.
                    let same_parameter = ed_data.map(|e| e.same_parameter).unwrap_or(false);
                    let same_range = ed_data.map(|e| e.same_range).unwrap_or(false);
                    // OCCT L316-328.
                    if !same_parameter || !same_range {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(s)
                            .expect("InContext: the context must be bound");
                        if !same_parameter {
                            brep_check_add(lst, BRepCheckStatus::InvalidSameParameterFlag);
                        }
                        if !same_range {
                            brep_check_add(lst, BRepCheckStatus::InvalidSameRangeFlag);
                        }
                        return;
                    }
                    // OCCT L330-331.
                    let hcurve = self.my_h_curve.clone();
                    let first = hcurve.as_ref().map(|h| h.first_parameter()).unwrap_or(0.0);
                    let last = hcurve.as_ref().map(|h| h.last_parameter()).unwrap_or(0.0);

                    // OCCT L333-340.
                    let fd = s.as_face();
                    let floc = s.location;
                    let tfloc = fd.map(|f| f.surface_location).unwrap_or(0);
                    let su = fd.and_then(|f| f.surface.clone());
                    // OCCT L337: L = (Floc * TFloc).Predivided(myShape.Location()).
                    #[allow(unused_assignments)]
                    let l_mat = location_matrix(brep, floc, tfloc)
                        * brep.get_location(my_shape.location).inverse();
                    let _ = l_mat;
                    // OCCT L338-339: LE = myShape.Location() * myCref->Location();
                    // Etrsf = LE.Transformation().
                    let etrsf = location_matrix(
                        brep,
                        my_shape.location,
                        self.my_cref.as_ref().map(cref_loc).unwrap_or(0),
                    );
                    let mut pcurvefound = false;

                    // OCCT L342-343.
                    let eps = rcad_kernel::precision::PCONFUSION;
                    let to_run_parallel = false; // myIsParallel (single-threaded analyzer)
                    let reps = edge_curve_reps(brep, &my_shape);
                    for cr in &reps {
                        // OCCT L348: cr != myCref && cr->IsCurveOnSurface(Su, L).
                        let same_as_cref = match (&self.my_cref, cr) {
                            (Some(c), r) => cref_same(c, r),
                            _ => false,
                        };
                        if !same_as_cref
                            && cr.is_curve_on_surface()
                            && self.rep_is_on_face(brep, cr, s)
                        {
                            pcurvefound = true;
                            // OCCT L351-353: GC->Range(f, l).
                            let (f, l) = match cr {
                                EdgeCurveRep::CurveOnSurface { range, .. }
                                | EdgeCurveRep::CurveOnClosedSurface { range, .. } => {
                                    (range[0], range[1])
                                }
                                _ => (0.0, 0.0),
                            };
                            let mut ff = f;
                            let mut ll = l;
                            // OCCT L355-359.
                            let cref_is_3d = self
                                .my_cref
                                .as_ref()
                                .map(|c| c.is_curve3d())
                                .unwrap_or(false);
                            if cref_is_3d {
                                // TransformedParameter(f/l, Etrsf) — identity
                                // for the rigid locations rcad materialises
                                // (the kernel transformed_parameter default).
                                ff = f;
                                ll = l;
                                let _ = etrsf;
                            }
                            {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(s)
                                    .expect("InContext: the context must be bound");
                                // OCCT L362-367.
                                if (ff - first).abs() > eps || (ll - last).abs() > eps {
                                    brep_check_add(lst, BRepCheckStatus::InvalidSameRangeFlag);
                                    brep_check_add(
                                        lst,
                                        BRepCheckStatus::InvalidSameParameterFlag,
                                    );
                                }
                            }
                            // OCCT L369-388: pc = cr->PCurve(); periodic/range.
                            let pc: Option<rcad_kernel::geom::Curve2d> = match cr {
                                EdgeCurveRep::CurveOnSurface { pcurve, .. } => Some(pcurve.clone()),
                                EdgeCurveRep::CurveOnClosedSurface { pcurve1, .. } => {
                                    Some(pcurve1.clone())
                                }
                                _ => None,
                            };
                            if let Some(pc) = &pc {
                                let mut is_periodic = Curve2dEval::is_periodic(pc);
                                let mut a_period = f64::MAX;
                                if is_periodic {
                                    a_period = period_of_curve2d(pc);
                                }
                                let dom = Curve2dEval::default_domain(pc);
                                let mut fp = dom[0];
                                let mut lp = dom[1];
                                if let rcad_kernel::geom::Curve2d::Trimmed(tc) = pc {
                                    let a_c: &rcad_kernel::geom::Curve2d = &tc.curve;
                                    let bdom = Curve2dEval::default_domain(a_c);
                                    fp = bdom[0];
                                    lp = bdom[1];
                                    is_periodic = Curve2dEval::is_periodic(a_c);
                                    if is_periodic {
                                        a_period = period_of_curve2d(a_c);
                                    }
                                }
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(s)
                                    .expect("InContext: the context must be bound");
                                if is_periodic && (l - f > a_period + eps) {
                                    // OCCT L389-393.
                                    brep_check_add(lst, BRepCheckStatus::InvalidRange);
                                    return;
                                } else if !is_periodic && (f < fp - eps || l > lp + eps) {
                                    // OCCT L394-398.
                                    brep_check_add(lst, BRepCheckStatus::InvalidRange);
                                    return;
                                }
                            }

                            // OCCT L400-455: the geometric control.
                            if self.my_gctrl {
                                // OCCT L402-405: Sb = Su->Transformed((Floc*TFloc)).
                                let sb = su.as_ref().map(|su| {
                                    transform_surface(su, &location_matrix(brep, floc, tfloc))
                                });
                                let pc: Option<rcad_kernel::geom::Curve2d> = match cr {
                                    EdgeCurveRep::CurveOnSurface { pcurve, .. } => {
                                        Some(pcurve.clone())
                                    }
                                    EdgeCurveRep::CurveOnClosedSurface { pcurve1, .. } => {
                                        Some(pcurve1.clone())
                                    }
                                    _ => None,
                                };
                                if let (Some(sb), Some(pc)) = (sb, pc) {
                                    let gahs = GeomAdaptorSurface::new(sb);
                                    // OCCT L406-432: the first run over PCurve1.
                                    let ghpc = Geom2dAdaptorCurve::new(pc, f, l);
                                    let acs = Adaptor3dCurveOnSurface::new(ghpc, gahs.clone());
                                    self.validate_edge_run(
                                        s,
                                        hcurve.as_ref(),
                                        &acs,
                                        same_parameter,
                                        tol,
                                        to_run_parallel,
                                        RunKind::First {
                                            closed: cref_is_closed(cr),
                                        },
                                    );
                                    // OCCT L433-454: the second run over PCurve2
                                    // (same bounds) when the representation is a
                                    // closed-surface one.
                                    if cref_is_closed(cr) {
                                        if let EdgeCurveRep::CurveOnClosedSurface {
                                            pcurve2, ..
                                        } = cr
                                        {
                                            let ghpc2 =
                                                Geom2dAdaptorCurve::new(pcurve2.clone(), f, l);
                                            // OCCT L436: ACS->Load(GHPC, GAHS) —
                                            // sans doute inutile.
                                            let acs2 =
                                                Adaptor3dCurveOnSurface::new(ghpc2, gahs.clone());
                                            self.validate_edge_run(
                                                s,
                                                hcurve.as_ref(),
                                                &acs2,
                                                same_parameter,
                                                tol,
                                                to_run_parallel,
                                                RunKind::ClosedSecond,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                        // OCCT L457: itcr.Next().
                    }

                    // OCCT L460-515: no pcurve found — the on-the-fly plane
                    // projection.
                    if !pcurvefound {
                        let su_basis = match &su {
                            Some(Surface3::Trimmed(ts)) => Some((*ts.basis).clone()),
                            other => other.clone(),
                        };
                        let plane = match su_basis {
                            Some(rcad_kernel::geom::Surface3::Plane(p)) => Some(p),
                            _ => None,
                        };
                        match plane {
                            None => {
                                // OCCT L473-476: not a plane.
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(s)
                                    .expect("InContext: the context must be bound");
                                brep_check_add(lst, BRepCheckStatus::NoCurveOnSurface);
                            }
                            Some(_p) => {
                                // OCCT L480-514: on fait la projection a la volee.
                                if self.my_gctrl {
                                    // OCCT L482-483: P transformed by (Floc*TFloc).
                                    let p_t = match &su {
                                        Some(rcad_kernel::geom::Surface3::Plane(p)) => {
                                            match transform_surface(
                                                &rcad_kernel::geom::Surface3::Plane(p.clone()),
                                                &location_matrix(brep, floc, tfloc),
                                            ) {
                                                rcad_kernel::geom::Surface3::Plane(pt) => pt,
                                                _ => _p.clone(),
                                            }
                                        }
                                        _ => _p.clone(),
                                    };
                                    // OCCT L487-500: the 3D reference projected on
                                    // the plane (rcad: the BRep_Tool::CurveOnPlane
                                    // translation returns the projected pcurve
                                    // directly).
                                    let proj_pc: Option<rcad_kernel::geom::Curve2d> = match &hcurve
                                    {
                                        Some(HCurveAdaptor::Curve3d(gac)) => curve_on_plane(
                                            gac.curve(),
                                            [first, last],
                                            &p_t,
                                        ),
                                        _ => None,
                                    };
                                    match proj_pc {
                                        None => {
                                            // GAP: GeomProjLib::ProjectOnPlane on
                                            // the reference curve (Edge.cxx
                                            // L490-498) — the projection failed;
                                            // neutral (no status set).
                                        }
                                        Some(pc) => {
                                            let gahs =
                                                GeomAdaptorSurface::new(Surface3::Plane(p_t));
                                            let ghpc =
                                                Geom2dAdaptorCurve::new(pc, first, last);
                                            let acs =
                                                Adaptor3dCurveOnSurface::new(ghpc, gahs);
                                            self.validate_edge_run(
                                                s,
                                                hcurve.as_ref(),
                                                &acs,
                                                same_parameter,
                                                tol,
                                                to_run_parallel,
                                                RunKind::Projection,
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ShapeType::Solid => {
                // OCCT L518-546: the edge must be connected exactly twice.
                let mut nbconnection: i32 = 0;
                for fac in explorer(brep, s, ShapeType::Face) {
                    for e2 in explorer(brep, &fac, ShapeType::Edge) {
                        if e2.is_same(&my_shape) {
                            nbconnection += 1;
                        }
                    }
                }
                let degenerated = ed_data.map(|e| e.degenerated).unwrap_or(false);
                let lst = self
                    .base
                    .my_map
                    .find_mut(s)
                    .expect("InContext: the context must be bound");
                if nbconnection < 2 && !degenerated {
                    brep_check_add(lst, BRepCheckStatus::FreeEdge);
                } else if nbconnection > 2 {
                    brep_check_add(lst, BRepCheckStatus::InvalidMultiConnexity);
                } else {
                    brep_check_add(lst, BRepCheckStatus::NoError);
                }
            }
            _ => {
                // OCCT L548-549: default — nothing.
            }
        }
        // OCCT L551-554.
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");
        if lst.is_empty() {
            lst.push(BRepCheckStatus::NoError);
        }
    }

    /// The three OCCT call sites of BRepLib_ValidateEdge in InContext(FACE),
    /// distinguished by their failure branches:
    /// - `First` (Edge.cxx L416-432): InvalidCurveOnClosedSurface when the
    ///   representation is a seam one, else InvalidCurveOnSurface; then
    ///   InvalidSameParameterFlag unconditionally.
    /// - `ClosedSecond` (Edge.cxx L438-453): InvalidCurveOnClosedSurface;
    ///   InvalidSameParameterFlag only when SameParameter.
    /// - `Projection` (Edge.cxx L504-512): InvalidCurveOnSurface only.
    #[allow(clippy::too_many_arguments)]
    fn validate_edge_run(
        &mut self,
        s: &Shape,
        hcurve: Option<&HCurveAdaptor>,
        acs: &Adaptor3dCurveOnSurface,
        same_parameter: bool,
        tol: f64,
        to_run_parallel: bool,
        run: RunKind,
    ) {
        let Some(HCurveAdaptor::Curve3d(gac)) = hcurve else {
            // GAP: the reference curve is a pcurve adaptor (Adaptor3d_
            // CurveOnSurface) — the rcad BRepLibValidateEdge translation
            // accepts only GeomAdaptorCurve references — neutral (no status).
            return;
        };
        let mut a_validate_edge = BRepLibValidateEdge::new(gac.clone(), acs.clone(), same_parameter);
        a_validate_edge.set_exit_if_tolerance_exceeded(tol);
        a_validate_edge.set_exact_method(self.my_is_exact_method);
        a_validate_edge.set_parallel(to_run_parallel);
        a_validate_edge.process();
        let failed = !a_validate_edge.is_done() || !a_validate_edge.check_tolerance(tol);
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");
        if failed {
            match run {
                RunKind::First { closed } => {
                    // OCCT L418-425.
                    if closed {
                        brep_check_add(lst, BRepCheckStatus::InvalidCurveOnClosedSurface);
                    } else {
                        brep_check_add(lst, BRepCheckStatus::InvalidCurveOnSurface);
                    }
                    // OCCT L427.
                    brep_check_add(lst, BRepCheckStatus::InvalidSameParameterFlag);
                }
                RunKind::ClosedSecond => {
                    // OCCT L446.
                    brep_check_add(lst, BRepCheckStatus::InvalidCurveOnClosedSurface);
                    // OCCT L448-451.
                    if same_parameter {
                        brep_check_add(lst, BRepCheckStatus::InvalidSameParameterFlag);
                    }
                }
                RunKind::Projection => {
                    // OCCT L509-512.
                    brep_check_add(lst, BRepCheckStatus::InvalidCurveOnSurface);
                }
            }
        }
    }

    /// OCCT `cr->IsCurveOnSurface(Su, L)` — the representation matches the
    /// face and the composed location. rcad matches the representation's
    /// owning-face key against the OCCT pcurve key
    /// `(face TShape, L.Predivided(E.Location()))`.
    fn rep_is_on_face(&self, brep: &BRep, cr: &EdgeCurveRep, face: &Shape) -> bool {
        let face_key = match cr {
            EdgeCurveRep::CurveOnSurface { face_key, .. }
            | EdgeCurveRep::CurveOnClosedSurface { face_key, .. } => *face_key,
            _ => return false,
        };
        let expected_loc =
            brep.compose_pcurve_location(face.location, self.base.my_shape.location);
        face_key.0 == face.ptr_id() && face_key.1 == expected_loc
    }

    /// OCCT BRepCheck_Edge::Blind (Edge.cxx L559-568) — the body was removed
    /// upstream because of its uselessness.
    pub fn blind(&mut self) {
        if !self.base.my_blind {
            self.base.my_blind = true;
        }
    }

    /// OCCT BRepCheck_Edge::GeometricControls(B) (Edge.cxx L572-575).
    pub fn set_geometric_controls(&mut self, b: bool) {
        self.my_gctrl = b;
    }

    /// OCCT BRepCheck_Edge::GeometricControls() (Edge.cxx L579-582).
    pub fn geometric_controls(&self) -> bool {
        self.my_gctrl
    }

    /// OCCT BRepCheck_Edge::SetExactMethod (Edge.hxx L54).
    pub fn set_exact_method(&mut self, the_is_exact: bool) {
        self.my_is_exact_method = the_is_exact;
    }

    /// OCCT BRepCheck_Edge::SetStatus (Edge.cxx L586-594).
    pub fn set_status(&mut self, the_status: BRepCheckStatus) {
        let lst = self
            .base
            .my_map
            .find_mut(&self.base.my_shape)
            .expect("SetStatus: myShape must be bound");
        brep_check_add(lst, the_status);
    }

    /// OCCT BRepCheck_Edge::CheckPolygonOnTriangulation (Edge.cxx L711-840).
    /// rcad BRep edges carry no Poly_PolygonOnTriangulation representation,
    /// so the representation walk finds none and the function returns
    /// BRepCheck_NoError — exactly the OCCT !aHasPolygonOnTriangulation path
    /// (Edge.cxx L737-740).
    pub fn check_polygon_on_triangulation(&self, _brep: &BRep, the_edge: &Shape) -> BRepCheckStatus {
        // OCCT L713-735: the representation walk over theEdge.TShape().
        let reps = edge_curve_reps(_brep, the_edge);
        let a_has_polygon_on_triangulation = false;
        let mut a_has_curve3d = false;
        for a_cr in &reps {
            // OCCT L722: aCR->IsPolygonOnTriangulation() — no rcad
            // representation kind carries a polygon on triangulation.
            let _ = a_cr;
            // OCCT L726: aCR->IsCurve3D() && aCR->Curve3D() != NULL.
            if a_cr.is_curve3d()
                && matches!(a_cr, EdgeCurveRep::Curve3D { .. })
            {
                a_has_curve3d = true;
            }
            if a_has_polygon_on_triangulation && a_has_curve3d {
                break;
            }
        }
        // OCCT L737-740.
        if !a_has_polygon_on_triangulation || !a_has_curve3d {
            return BRepCheckStatus::NoError;
        }
        // OCCT L742-837: unreachable for rcad BRep data (no polygon on
        // triangulation representations exist).
        BRepCheckStatus::NoError
    }
}

/// The location index carried by a representation.
fn cref_loc(cr: &EdgeCurveRep) -> u32 {
    match cr {
        EdgeCurveRep::Curve3D { location, .. } => *location,
        EdgeCurveRep::CurveOnSurface { location, .. } => *location,
        EdgeCurveRep::CurveOnClosedSurface { location, .. } => *location,
        _ => 0,
    }
}

/// Representation identity (OCCT `cr != myCref` pointer comparison).
fn cref_same(a: &EdgeCurveRep, b: &EdgeCurveRep) -> bool {
    match (a, b) {
        (EdgeCurveRep::Curve3D { .. }, EdgeCurveRep::Curve3D { .. }) => true,
        (
            EdgeCurveRep::CurveOnSurface { face_key: k1, pcurve: p1, .. },
            EdgeCurveRep::CurveOnSurface { face_key: k2, pcurve: p2, .. },
        ) => k1 == k2 && pcurve_same(p1, p2),
        (
            EdgeCurveRep::CurveOnClosedSurface { face_key: k1, pcurve1: p1, .. },
            EdgeCurveRep::CurveOnClosedSurface { face_key: k2, pcurve1: p2, .. },
        ) => k1 == k2 && pcurve_same(p1, p2),
        _ => false,
    }
}

fn pcurve_same(a: &rcad_kernel::geom::Curve2d, b: &rcad_kernel::geom::Curve2d) -> bool {
    // Value identity for the pcurve (the rcad Curve2d carries no handle).
    format!("{a:?}") == format!("{b:?}")
}

/// Whether the representation is a closed-surface (seam) one.
fn cref_is_closed(cr: &EdgeCurveRep) -> bool {
    cr.is_curve_on_closed_surface()
}

/// The period of a periodic 3D curve (OCCT `C3d->Period()`): the natural
/// domain length (circle/ellipse: 2*PI).
fn period_of_curve3d(c: &rcad_kernel::geom::Curve3) -> f64 {
    let dom = CurveEval::default_domain(c);
    dom[1] - dom[0]
}

/// The period of a periodic 2D curve (OCCT `PC->Period()`).
fn period_of_curve2d(c: &rcad_kernel::geom::Curve2d) -> f64 {
    let dom = Curve2dEval::default_domain(c);
    dom[1] - dom[0]
}
