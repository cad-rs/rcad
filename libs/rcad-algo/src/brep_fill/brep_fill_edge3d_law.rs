//! OCCT BRepFill_Edge3DLaw / BRepFill_ACRLaw / BRepFill_EdgeOnSurfLaw /
//! BRepFill_DraftLaw (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_Edge3DLaw.hxx+cxx (whole files), BRepFill_ACRLaw.hxx+cxx (whole
//! files), BRepFill_EdgeOnSurfLaw.hxx+cxx (whole files),
//! BRepFill_DraftLaw.hxx+cxx (whole files).
//!
//! Architecture differences:
//! - The BRepFill_LocationLaw inheritance maps to base-struct embedding: the
//!   derived classes hold `pub base: BRepFillLocationLaw` (or the Edge3DLaw
//!   level for DraftLaw) and expose it through the
//!   [`BRepFillLocationLawOps`] trait (Rust has no inheritance; none of the
//!   four OCCT classes overrides a base method).
//! - `BRepTools_WireExplorer` maps to the wire-ordered `TWireData::edges`;
//!   `BRep_Tool::Curve / Degenerated` map to the stored `TEdgeData` fields.
//! - `handle(GeomFill_LocationLaw) Law->Copy()` maps to
//!   `Rc::new(RefCell::new(law.borrow().copy_law()))`.
//! - `GeomAdaptor_Curve(C, First, Last)` maps to a trimmed view of the rcad
//!   Curve3 (the adaptor restricts the evaluation range).
//! - `GeomFill_LocationGuide` is the batch-2 rcad type; `SetOrigine` /
//!   `Copy` map to `set_origine` / `copy_law`.
//! - The EdgeOnSurfLaw per-edge
//!   `myLaws->ChangeValue(ipath)->SetCurve(AC)` consumes an
//!   Adaptor3d_CurveOnSurface handle — the batch-1 LocationLaw trait carries
//!   the adaptor as Curve3, so the CurveOnSurface-backed set keeps the OCCT
//!   failure path at that call site (GAP, see the body).

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::topo::topods::{BRep, Orientation, Shape};

use crate::brep_fill::brep_fill_location_law::{BRepFillLocationLaw, BRepFillLocationLawOps};
use crate::brep_fill::brep_fill_section_law::{brep_tool_curve, reversed_curve};
use crate::brep_fill::compatible_wires::{compute_acr, wire_edges};
use crate::geomalgo::geomfill::curve_and_trihedron::CurveAndTrihedron;
use crate::geomalgo::geomfill::darboux::Darboux;
use crate::geomalgo::geomfill::gp_mat::GpMat;
use crate::geomalgo::geomfill::location_guide::LocationGuide;
use crate::geomalgo::geomfill::location_law::LocationLaw;
use crate::geomalgo::geomfill::trihedron_law::{curve_first_parameter, curve_last_parameter};

/// OCCT handle(GeomFill_LocationLaw) shared-cell wrap of a `Law->Copy()`:
/// the batch-1 `copy_law()` returns `Box<dyn LocationLaw>`; this forwarding
/// impl makes the boxed copy a `LocationLaw` itself, so the copy re-shares
/// as `Rc<RefCell<dyn LocationLaw>>` (the OCCT Copy-into-handle form).
impl LocationLaw for Box<dyn LocationLaw> {
    fn set_curve(&mut self, c: rcad_kernel::geom::Curve3) -> bool {
        self.as_mut().set_curve(c)
    }
    fn get_curve(&self) -> Option<rcad_kernel::geom::Curve3> {
        self.as_ref().get_curve()
    }
    fn set_trsf(&mut self, transfo: GpMat) {
        self.as_mut().set_trsf(transfo)
    }
    fn copy_law(&self) -> Box<dyn LocationLaw> {
        self.as_ref().copy_law()
    }
    fn d0(&self, param: f64, m: &mut GpMat, v: &mut glam::DVec3) -> bool {
        self.as_ref().d0(param, m, v)
    }
    fn d0_2d(&self, param: f64, m: &mut GpMat, v: &mut glam::DVec3, pnts2d: &mut [glam::DVec2]) -> bool {
        self.as_ref().d0_2d(param, m, v, pnts2d)
    }
    fn d1(
        &self,
        param: f64,
        m: &mut GpMat,
        v: &mut glam::DVec3,
        dm: &mut GpMat,
        dv: &mut glam::DVec3,
        pnts2d: &mut [glam::DVec2],
        vecs2d: &mut [glam::DVec2],
    ) -> bool {
        self.as_ref().d1(param, m, v, dm, dv, pnts2d, vecs2d)
    }
    fn nb_intervals(&self, s: rcad_kernel::math::GeomAbsShape) -> usize {
        self.as_ref().nb_intervals(s)
    }
    fn intervals(&self, t: &mut Vec<f64>, s: rcad_kernel::math::GeomAbsShape) {
        self.as_ref().intervals(t, s)
    }
    fn set_interval(&mut self, first: f64, last: f64) {
        self.as_mut().set_interval(first, last)
    }
    fn get_interval(&self, first: &mut f64, last: &mut f64) {
        self.as_ref().get_interval(first, last)
    }
    fn get_domain(&self, first: &mut f64, last: &mut f64) {
        self.as_ref().get_domain(first, last)
    }
    fn get_maximal_norm(&self) -> f64 {
        self.as_ref().get_maximal_norm()
    }
    fn get_average_law(&self, am: &mut GpMat, av: &mut glam::DVec3) {
        self.as_ref().get_average_law(am, av)
    }
}

/// OCCT GeomAdaptor_Curve(C, First, Last) — the rcad adaptor restricts the
/// evaluation range through a trimmed view.
fn restrict_to_range(
    c: &rcad_kernel::geom::Curve3,
    first: f64,
    last: f64,
) -> rcad_kernel::geom::Curve3 {
    rcad_kernel::geom::Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
        curve: Box::new(c.clone()),
        first,
        last,
    })
}

/// OCCT gp_Dir::IsParallel(Other, AngularTolerance).
fn gp_dir_is_parallel(a: glam::DVec3, b: glam::DVec3, ang_tol: f64) -> bool {
    let la = a.length();
    let lb = b.length();
    if la < 1e-300 || lb < 1e-300 {
        return false;
    }
    let cross = a.cross(b).length() / (la * lb);
    cross <= ang_tol.abs().sin() + 1e-18
}

/// The empty base state before Init (the default-constructed OCCT members).
fn empty_base() -> BRepFillLocationLaw {
    BRepFillLocationLaw {
        my_path: Shape::null(),
        my_tol: 0.0,
        my_laws: Vec::new(),
        my_length: Vec::new(),
        my_edges: Vec::new(),
        my_disc: None,
        my_type: 0,
    }
}

// ---------------------------------------------------------------------------
// BRepFill_Edge3DLaw
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Edge3DLaw (hxx L19-30): Build Location Law, with a Wire.
pub struct BRepFillEdge3DLaw {
    /// OCCT BRepFill_LocationLaw base sub-object.
    pub base: BRepFillLocationLaw,
}

impl BRepFillEdge3DLaw {
    /// OCCT BRepFill_Edge3DLaw(Path, Law) (cxx L43-79).
    pub fn new(brep: &BRep, path: &Shape, law: Rc<RefCell<dyn LocationLaw>>) -> Self {
        let mut this = BRepFillEdge3DLaw { base: empty_base() };
        this.base.init(brep, path);

        let mut ipath = 0usize; // OCCT ipath starts at 0 and is pre-incremented

        for e in wire_edges(brep, path) {
            // E = wexp.Current(); if (!BRep_Tool::Degenerated(E))
            if !brep.edge(e.clone()).degenerated {
                ipath += 1;
                this.base.my_edges[ipath - 1] = e.clone();
                let (mut c, mut first, mut last) = match brep_tool_curve(brep, &e) {
                    Some(v) => v,
                    None => continue,
                };
                let or = e.orientation;
                if or == Orientation::Reversed {
                    // Geom_TrimmedCurve CBis(C, First, Last); CBis->Reverse();
                    // ("Pour eviter de deteriorer la topologie")
                    c = reversed_curve(&c);
                    // First = C->FirstParameter(); Last = C->LastParameter();
                    first = curve_first_parameter(&c);
                    last = curve_last_parameter(&c);
                }

                // AC = new GeomAdaptor_Curve(C, First, Last);
                // myLaws->SetValue(ipath, Law->Copy());
                // myLaws->ChangeValue(ipath)->SetCurve(AC);
                let copy: Rc<RefCell<dyn LocationLaw>> = Rc::new(RefCell::new(law.borrow().copy_law()));
                copy.borrow_mut().set_curve(restrict_to_range(&c, first, last));
                this.base.my_laws[ipath - 1] = copy;
            }
        }
        this
    }
}

/// OCCT BRepFill_LocationLaw inheritance (no override).
impl BRepFillLocationLawOps for BRepFillEdge3DLaw {
    fn base(&self) -> &BRepFillLocationLaw {
        &self.base
    }
    fn base_mut(&mut self) -> &mut BRepFillLocationLaw {
        &mut self.base
    }
}

// ---------------------------------------------------------------------------
// BRepFill_ACRLaw
// ---------------------------------------------------------------------------

/// OCCT BRepFill_ACRLaw (hxx L29-40): Build Location Law, with a Wire. In
/// the case of guided contour and trihedron by reduced curvilinear abscissa.
pub struct BRepFillACRLaw {
    /// OCCT BRepFill_LocationLaw base sub-object.
    pub base: BRepFillLocationLaw,
    /// OCCT handle(NCollection_HArray1<double>) OrigParam (private, the
    /// 0-based array of length NbEdge + 1).
    pub orig_param: Vec<f64>,
}

impl BRepFillACRLaw {
    /// OCCT BRepFill_ACRLaw(Path, theLaw) (cxx L41-103).
    pub fn new(brep: &BRep, path: &Shape, the_law: &Rc<RefCell<LocationGuide>>) -> Self {
        let mut this = BRepFillACRLaw {
            base: empty_base(),
            orig_param: Vec::new(),
        };
        this.base.init(brep, path);

        // calculate the nb of edge of the path (L52-58)
        let nb_edge = wire_edges(brep, path).len();

        // tab to memorize ACR for each edge
        // NCollection_Array1<double> Orig(0, NbEdge);
        // BRepFill::ComputeACR(Path, Orig);
        let orig = compute_acr(brep, path);

        // OrigParam->SetValue(0, 0);
        // for (ipath = 1; ipath <= NbEdge; ipath++)
        //   OrigParam->SetValue(ipath, Orig(ipath));
        this.orig_param = vec![0.0; nb_edge + 1];
        for ipath in 1..=nb_edge.min(orig.len().saturating_sub(1)) {
            this.orig_param[ipath] = orig[ipath];
        }

        // process each edge of the trajectory (L82-101)
        let mut ipath = 0usize;
        for e in wire_edges(brep, path) {
            if !brep.edge(e.clone()).degenerated {
                ipath += 1;
                this.base.my_edges[ipath - 1] = e.clone();
                let (mut c, mut first, mut last) = match brep_tool_curve(brep, &e) {
                    Some(v) => v,
                    None => continue,
                };
                let or = e.orientation;
                if or == Orientation::Reversed {
                    // Geom_TrimmedCurve CBis(C, First, Last); CBis->Reverse();
                    // ("To avoid damaging the topology")
                    c = reversed_curve(&c);
                    first = curve_first_parameter(&c);
                    last = curve_last_parameter(&c);
                }

                // Set the parameters for the case multi-edges
                // double t1 = OrigParam->Value(ipath - 1);
                // double t2 = OrigParam->Value(ipath);
                // Loc->SetOrigine(t1, t2);
                let t1 = this.orig_param[ipath - 1];
                let t2 = this.orig_param[ipath];
                the_law.borrow_mut().set_origine(t1, t2);

                // myLaws->SetValue(ipath, Loc->Copy());
                // myLaws->ChangeValue(ipath)->SetCurve(AC);
                let copy: Rc<RefCell<dyn LocationLaw>> = Rc::new(RefCell::new(the_law.borrow().copy_law()));
                copy.borrow_mut().set_curve(restrict_to_range(&c, first, last));
                this.base.my_laws[ipath - 1] = copy;
            }
        }
        this
    }
}

/// OCCT BRepFill_LocationLaw inheritance (no override).
impl BRepFillLocationLawOps for BRepFillACRLaw {
    fn base(&self) -> &BRepFillLocationLaw {
        &self.base
    }
    fn base_mut(&mut self) -> &mut BRepFillLocationLaw {
        &mut self.base
    }
}

// ---------------------------------------------------------------------------
// BRepFill_EdgeOnSurfLaw
// ---------------------------------------------------------------------------

/// OCCT BRepFill_EdgeOnSurfLaw (hxx L23-35): Build Location Law, with a
/// Wire and a Surface.
pub struct BRepFillEdgeOnSurfLaw {
    /// OCCT BRepFill_LocationLaw base sub-object.
    pub base: BRepFillLocationLaw,
    /// OCCT bool hasresult (private).
    pub hasresult: bool,
}

impl BRepFillEdgeOnSurfLaw {
    /// OCCT BRepFill_EdgeOnSurfLaw(Path, Surf) (cxx L36-95).
    pub fn new(brep: &BRep, path: &Shape, surf: &Shape) -> Self {
        let mut this = BRepFillEdgeOnSurfLaw {
            base: empty_base(),
            hasresult: true,
        };
        this.base.init(brep, path);

        // GeomFill_Darboux TLaw = new GeomFill_Darboux();
        // GeomFill_CurveAndTrihedron Law = new GeomFill_CurveAndTrihedron(TLaw);
        let t_law = Darboux::new();
        let law = Rc::new(RefCell::new(CurveAndTrihedron::new(Box::new(t_law))));

        let mut ipath = 0usize;
        for e in wire_edges(brep, path) {
            if !brep.edge(e.clone()).degenerated {
                ipath += 1;
                this.base.my_edges[ipath - 1] = e.clone();

                // for (Trouve = false, exp.Init(Surf, TopAbs_FACE);
                //      exp.More() && !Trouve; exp.Next())
                //   C = BRep_Tool::CurveOnSurface(E, F, First, Last);
                let mut trouve = false;
                for f in explorer_faces(surf) {
                    if brep_tool_curve_on_surface(brep, &e, &f).is_some() {
                        trouve = true;
                        break;
                    }
                }
                if !trouve {
                    // Impossible to construct the law.
                    this.hasresult = false;
                    this.base.my_laws.clear();
                    return this;
                }

                // (the REVERSED Geom2d_TrimmedCurve reversal of cxx L78-84 is
                // carried by the CurveOnSurface GAP below)

                // AC2d = new Geom2dAdaptor_Curve(C, First, Last);
                // AC = new Adaptor3d_CurveOnSurface(AC2d, AS);
                // myLaws->SetValue(ipath, Law->Copy());
                // myLaws->ChangeValue(ipath)->SetCurve(AC);
                let copy: Rc<RefCell<dyn LocationLaw>> = Rc::new(RefCell::new(law.borrow().copy_law()));
                gap_curve_on_surface_set(&copy);
                this.base.my_laws[ipath - 1] = copy;
            }
        }
        this
    }

    /// OCCT HasResult (cxx L98-101).
    pub fn has_result(&self) -> bool {
        self.hasresult
    }
}

/// OCCT TopExp_Explorer(Surf, TopAbs_FACE) — the faces of the support.
fn explorer_faces(surf: &Shape) -> Vec<Shape> {
    use crate::feat::brep_feat_builder::explorer;
    use rcad_kernel::topo::topods::ShapeType;
    explorer(surf, ShapeType::Face, ShapeType::Shape)
}

/// OCCT BRep_Tool::CurveOnSurface(E, F, First, Last) — the stored pcurve.
fn brep_tool_curve_on_surface(
    brep: &BRep,
    e: &Shape,
    f: &Shape,
) -> Option<(rcad_kernel::geom::Curve2d, f64, f64)> {
    let ed = brep.edge(e.clone());
    let key = brep.pcurve_key(f.index);
    ed.pcurves
        .get(&key)
        .map(|(c, f0, l0)| (c.clone(), *f0, *l0))
}

/// GAP: GeomFill_CurveAndTrihedron::SetCurve(handle(Adaptor3d_Curve)) with a
/// CurveOnSurface — the batch-1 LocationLaw trait carries the adaptor as
/// Curve3 and the CurveOnSurface-backed set (BRepFill_EdgeOnSurfLaw.cxx
/// L92) keeps the OCCT failure path (the darboux.rs set_curve_on_surface
/// route is the pending re-home).
fn gap_curve_on_surface_set(_law: &Rc<RefCell<dyn LocationLaw>>) {
    panic!(
        "GAP: Adaptor3d_CurveOnSurface SetCurve (TKG3d) is not translated — \
         BRepFill_EdgeOnSurfLaw ctor (BRepFill_EdgeOnSurfLaw.cxx L92)"
    )
}

/// OCCT BRepFill_LocationLaw inheritance (no override).
impl BRepFillLocationLawOps for BRepFillEdgeOnSurfLaw {
    fn base(&self) -> &BRepFillLocationLaw {
        &self.base
    }
    fn base_mut(&mut self) -> &mut BRepFillLocationLaw {
        &mut self.base
    }
}

// ---------------------------------------------------------------------------
// BRepFill_DraftLaw
// ---------------------------------------------------------------------------

/// OCCT BRepFill_DraftLaw (hxx L19-29): Build Location Law, with a Wire
/// (the draft specialization of BRepFill_Edge3DLaw).
pub struct BRepFillDraftLaw {
    /// OCCT BRepFill_Edge3DLaw sub-object (the inheritance chain
    /// BRepFill_DraftLaw -> BRepFill_Edge3DLaw -> BRepFill_LocationLaw).
    pub edge3d: BRepFillEdge3DLaw,
}

/// OCCT static ToG0 (BRepFill_DraftLaw.cxx L29-33) — Calculate a
/// transformation T tq T.M2 = M1 (the file-local duplicate of the
/// BRepFill_LocationLaw static).
fn to_g0(m1: &GpMat, m2: &GpMat) -> GpMat {
    // T = M2.Inverted(); T *= M1;
    let mut res = GpMat::identity();
    crate::brep_fill::brep_fill_location_law::to_g0(m1, m2, &mut res);
    res
}

impl BRepFillDraftLaw {
    /// OCCT BRepFill_DraftLaw(Path, Law) (cxx L41-44) — forwards to the
    /// BRepFill_Edge3DLaw constructor.
    pub fn new(brep: &BRep, path: &Shape, law: Rc<RefCell<dyn LocationLaw>>) -> Self {
        BRepFillDraftLaw {
            edge3d: BRepFillEdge3DLaw::new(brep, path, law),
        }
    }

    /// OCCT CleanLaw (cxx L46-80) — To clean the little discontinuities.
    pub fn clean_law(&mut self, tol_angular: f64) {
        let mut first = 0.0;
        let mut last = 0.0;
        let mut trsf = GpMat::identity();
        let mut m1 = GpMat::identity();
        let mut m2 = GpMat::identity();
        let mut v = glam::DVec3::ZERO;
        let mut t1;
        let mut t2;
        let mut n1;
        let mut n2;

        let laws_len = self.edge3d.base.my_laws.len();
        self.edge3d.base.my_laws[0]
            .borrow()
            .get_domain(&mut first, &mut last);

        let mut ipath = 2usize;
        while ipath <= laws_len {
            self.edge3d.base.my_laws[ipath - 2].borrow().d0(last, &mut m1, &mut v);
            self.edge3d.base.my_laws[ipath - 1]
                .borrow()
                .get_domain(&mut first, &mut last);
            self.edge3d.base.my_laws[ipath - 1].borrow().d0(first, &mut m2, &mut v);
            t1 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m1, 3);
            t2 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m2, 3);
            n1 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m1, 1);
            n2 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m2, 1);
            let _ = (t1, t2);
            if gp_dir_is_parallel(n1, n2, tol_angular) {
                // Correction G0 des normales...
                trsf = to_g0(&m1, &m2);
                self.edge3d.base.my_laws[ipath - 1].borrow_mut().set_trsf(trsf);
            }
            ipath += 1;
        }
    }
}

/// OCCT BRepFill_LocationLaw inheritance (no override).
impl BRepFillLocationLawOps for BRepFillDraftLaw {
    fn base(&self) -> &BRepFillLocationLaw {
        &self.edge3d.base
    }
    fn base_mut(&mut self) -> &mut BRepFillLocationLaw {
        &mut self.edge3d.base
    }
}
