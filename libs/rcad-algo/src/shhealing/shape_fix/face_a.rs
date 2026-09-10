//! 1:1 translation of OCCT `ShapeFix_Face` — part A
//! (`TKShHealing/ShapeFix/ShapeFix_Face.hxx` L17-301 +
//! `ShapeFix_Face.lxx` L18-119 + `ShapeFix_Face.cxx` L94-770).
//!
//! The impl submodules follow the OCCT continued-file convention:
//! [`super::face_b`] carries `FixAddNaturalBound` / `FixOrientation` (cxx
//! L769-1646) and [`super::face_c`] the rest of the fixing methods (cxx
//! L1650-3259).
//!
//! Architecture bridges:
//! 1. `ShapeFix_Root` inheritance — the rcad composition field `base` (the
//!    ShapeFix_Wire precedent).
//! 2. `BRep` pool argument — the rcad equivalents take `brep: &mut BRep`.
//! 3. `occ::handle<ShapeAnalysis_Surface> mySurf` -> `Option<...>` (the null
//!    surface of `Init(face)` maps to `None`); the OCCT null-handle pass to
//!    `ShapeFix_Wire::SetFace(face, null)` maps to the plain `SetFace(face)`
//!    overload (annotated at the site).
//! 4. `NCollection_DataMap<TopoDS_Shape, ...>` -> the maps keyed by the
//!    (TShape pointer, location) IsSame identity.
//! 5. `Bnd_Box2d` -> `rcad_kernel::math::bnd::BndBox2d`;
//!    `BndLib_Add2dCurve::Add` -> the [`super::intersection_tool::
//!    bnd_lib_add2d_curve`] sampler; `BRepTopAdaptor_FClass2d` -> the
//!    `topalgo::brep_top_adaptor::FClass2dTopol` translation.
//! 6. `BRepBuilderAPI_MakeFace(surf, tol)` (the natural-bound face with
//!    wires) — GAP: the kernel `add_tface` creates the face with the natural
//!    restriction flag and no wires; the natural-bound wire construction
//!    (BRepLib_MakeFace) is kernel scope, not translated (module doc of
//!    `face_c.rs`).

use std::collections::HashMap;

use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType};

use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes};
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::root::ShapeFixRoot;
use crate::shhealing::shape_fix::shape_fix::{least_edge_size, MessageProgressRange};
use crate::shhealing::shape_fix::wire::ShapeFixWire;

/// OCCT `NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
/// TopTools_ShapeMapHasher>` — the rcad wire-to-inner-wires map (the
/// IsSame identity key).
pub(crate) type WireListMap = HashMap<(u64, u32), Vec<Shape>>;

/// OCCT `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape,
/// TopTools_ShapeMapHasher>` — the rcad shape-to-shape map.
pub(crate) type ShapeShapeMap = HashMap<(u64, u32), Shape>;

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L94-101 — IsSurfaceUVInfinite (static).
// ---------------------------------------------------------------------------

/// OCCT static IsSurfaceUVInfinite (cxx L94-101).
pub(crate) fn is_surface_uv_infinite(the_surf: &rcad_kernel::geom::Surface3) -> bool {
    // OCCT L96-97: theSurf->Bounds(UMin, UMax, VMin, VMax).
    let dom = rcad_kernel::geom::SurfaceEval::default_domain(the_surf);
    let (u_min, u_max, v_min, v_max) = (dom[0], dom[1], dom[2], dom[3]);

    // OCCT L99-100: Precision::IsInfinite.
    let is_inf = |v: f64| v.abs() >= rcad_kernel::precision::INFINITE_VALUE;
    is_inf(u_min) || is_inf(u_max) || is_inf(v_min) || is_inf(v_max)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L103-107 — IsSurfaceUVPeriodic (static).
// ---------------------------------------------------------------------------

/// OCCT static IsSurfaceUVPeriodic (cxx L103-107): the GeomAdaptor_Surface
/// periodicity check re-hosted over the kernel surface evaluation.
pub(crate) fn is_surface_uv_periodic(the_surf: &rcad_kernel::geom::Surface3) -> bool {
    (rcad_kernel::geom::SurfaceEval::is_u_periodic(the_surf)
        && rcad_kernel::geom::SurfaceEval::is_v_periodic(the_surf))
        || matches!(the_surf, rcad_kernel::geom::Surface3::Sphere(_))
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L111-117 / hxx L256-295 — the constructor and the
// member fields.
// ---------------------------------------------------------------------------

/// OCCT ShapeFix_Face (hxx L52-296): performs various fixes on a face and
/// its wires — the ShapeFix_Wire fixes, wire orientation, natural bounds,
/// missing seam edge, and null-area wire detection.
pub struct ShapeFixFace {
    /// The OCCT ShapeFix_Root base subobject (bridge #1).
    pub base: ShapeFixRoot,
    /// OCCT mySurf (hxx L256).
    pub(crate) my_surf: Option<ShapeAnalysisSurface>,
    /// OCCT myFace (hxx L257).
    pub(crate) my_face: Shape,
    /// OCCT myResult (hxx L258).
    pub(crate) my_result: Shape,
    /// OCCT myFixWire (hxx L259).
    pub(crate) my_fix_wire: ShapeFixWire,
    /// OCCT myFwd (hxx L260).
    pub(crate) my_fwd: bool,
    /// OCCT myStatus (hxx L261).
    pub(crate) my_status: i32,
    /// OCCT myFixWireMode (hxx L285).
    pub(crate) my_fix_wire_mode: i32,
    /// OCCT myFixOrientationMode (hxx L286).
    pub(crate) my_fix_orientation_mode: i32,
    /// OCCT myFixAddNaturalBoundMode (hxx L287).
    pub(crate) my_fix_add_natural_bound_mode: i32,
    /// OCCT myFixMissingSeamMode (hxx L288).
    pub(crate) my_fix_missing_seam_mode: i32,
    /// OCCT myFixSmallAreaWireMode (hxx L289).
    pub(crate) my_fix_small_area_wire_mode: i32,
    /// OCCT myRemoveSmallAreaFaceMode (hxx L290).
    pub(crate) my_remove_small_area_face_mode: i32,
    /// OCCT myFixLoopWiresMode (hxx L291).
    pub(crate) my_fix_loop_wires_mode: i32,
    /// OCCT myFixIntersectingWiresMode (hxx L292).
    pub(crate) my_fix_intersecting_wires_mode: i32,
    /// OCCT myFixSplitFaceMode (hxx L293).
    pub(crate) my_fix_split_face_mode: i32,
    /// OCCT myAutoCorrectPrecisionMode (hxx L294).
    pub(crate) my_auto_correct_precision_mode: i32,
    /// OCCT myFixPeriodicDegenerated (hxx L295).
    pub(crate) my_fix_periodic_degenerated: i32,
}

impl Default for ShapeFixFace {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixFace {
    // OCCT ShapeFix_Face.cxx L111-117 — ShapeFix_Face().
    /// OCCT ShapeFix_Face::ShapeFix_Face() (cxx L111-117): creates an empty
    /// tool.
    pub fn new() -> Self {
        let mut this = ShapeFixFace {
            base: ShapeFixRoot::new(),
            my_surf: None,
            my_face: Shape::null(),
            my_result: Shape::null(),
            my_fix_wire: ShapeFixWire::new(),
            my_fwd: true,
            my_status: 0,
            my_fix_wire_mode: -1,
            my_fix_orientation_mode: -1,
            my_fix_add_natural_bound_mode: -1,
            my_fix_missing_seam_mode: -1,
            my_fix_small_area_wire_mode: -1,
            my_remove_small_area_face_mode: -1,
            my_fix_loop_wires_mode: -1,
            my_fix_intersecting_wires_mode: -1,
            my_fix_split_face_mode: -1,
            my_auto_correct_precision_mode: 1,
            my_fix_periodic_degenerated: -1,
        };
        this.clear_modes();
        this
    }

    // OCCT ShapeFix_Face.cxx L121-128 — ShapeFix_Face(face).
    /// OCCT ShapeFix_Face::ShapeFix_Face(face) (cxx L121-128): creates a
    /// tool and loads a face.
    pub fn with_face(brep: &mut BRep, face: &Shape) -> Self {
        let mut this = ShapeFixFace {
            base: ShapeFixRoot::new(),
            my_surf: None,
            my_face: Shape::null(),
            my_result: Shape::null(),
            my_fix_wire: ShapeFixWire::new(),
            my_fwd: true,
            my_status: 0,
            my_fix_wire_mode: -1,
            my_fix_orientation_mode: -1,
            my_fix_add_natural_bound_mode: -1,
            my_fix_missing_seam_mode: -1,
            my_fix_small_area_wire_mode: -1,
            my_remove_small_area_face_mode: -1,
            my_fix_loop_wires_mode: -1,
            my_fix_intersecting_wires_mode: -1,
            my_fix_split_face_mode: -1,
            my_auto_correct_precision_mode: 1,
            my_fix_periodic_degenerated: -1,
        };
        this.clear_modes();
        this.init_face(brep, face);
        this
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.lxx L18-116 — the inline accessors.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixWireMode (lxx L18-21).
    pub fn fix_wire_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_wire_mode
    }

    /// OCCT ShapeFix_Face::FixOrientationMode (lxx L25-28).
    pub fn fix_orientation_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_orientation_mode
    }

    /// OCCT ShapeFix_Face::FixAddNaturalBoundMode (lxx L32-35).
    pub fn fix_add_natural_bound_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_add_natural_bound_mode
    }

    /// OCCT ShapeFix_Face::FixMissingSeamMode (lxx L39-42).
    pub fn fix_missing_seam_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_missing_seam_mode
    }

    /// OCCT ShapeFix_Face::FixSmallAreaWireMode (lxx L46-49).
    pub fn fix_small_area_wire_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_small_area_wire_mode
    }

    /// OCCT ShapeFix_Face::RemoveSmallAreaFaceMode (lxx L53-56).
    pub fn remove_small_area_face_mode(&mut self) -> &mut i32 {
        &mut self.my_remove_small_area_face_mode
    }

    /// OCCT ShapeFix_Face::FixIntersectingWiresMode (lxx L60-63).
    pub fn fix_intersecting_wires_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_intersecting_wires_mode
    }

    /// OCCT ShapeFix_Face::FixLoopWiresMode (lxx L67-70).
    pub fn fix_loop_wires_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_loop_wires_mode
    }

    /// OCCT ShapeFix_Face::FixSplitFaceMode (lxx L74-77).
    pub fn fix_split_face_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_split_face_mode
    }

    /// OCCT ShapeFix_Face::AutoCorrectPrecisionMode (lxx L81-84).
    pub fn auto_correct_precision_mode(&mut self) -> &mut i32 {
        &mut self.my_auto_correct_precision_mode
    }

    /// OCCT ShapeFix_Face::FixPeriodicDegeneratedMode (lxx L88-91).
    pub fn fix_periodic_degenerated_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_periodic_degenerated
    }

    /// OCCT ShapeFix_Face::Face (lxx L95-98): returns a face which
    /// corresponds to the current state.
    pub fn face(&self) -> Shape {
        self.my_face.clone()
    }

    /// OCCT ShapeFix_Face::Result (lxx L102-105): returns the resulting
    /// shape (Face or Shell if split).
    pub fn result(&self) -> Shape {
        self.my_result.clone()
    }

    /// OCCT ShapeFix_Face::Status (lxx L109-112): the status of the last
    /// call to Perform.
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    /// OCCT ShapeFix_Face::FixWireTool (lxx L116-119): returns the tool for
    /// fixing wires.
    pub fn fix_wire_tool(&mut self) -> &mut ShapeFixWire {
        &mut self.my_fix_wire
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L132-145 — ClearModes.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::ClearModes (cxx L132-145): sets all modes to
    /// default.
    pub fn clear_modes(&mut self) {
        self.my_fix_wire_mode = -1;
        self.my_fix_orientation_mode = -1;
        self.my_fix_add_natural_bound_mode = -1;
        self.my_fix_missing_seam_mode = -1;
        self.my_fix_small_area_wire_mode = -1;
        self.my_remove_small_area_face_mode = -1;
        self.my_fix_intersecting_wires_mode = -1;
        self.my_fix_loop_wires_mode = -1;
        self.my_fix_split_face_mode = -1;
        self.my_auto_correct_precision_mode = 1;
        self.my_fix_periodic_degenerated = -1;
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L149-153 — SetMsgRegistrator.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::SetMsgRegistrator (cxx L149-153).
    pub fn set_msg_registrator(&mut self, msgreg: crate::shhealing::shape_fix::root::MsgRegistratorHandle) {
        self.base.set_msg_registrator(msgreg.clone());
        self.my_fix_wire.base.set_msg_registrator(msgreg);
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L157-161 — SetPrecision.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::SetPrecision (cxx L157-161): sets the basic
    /// precision value (also to FixWireTool).
    pub fn set_precision(&mut self, preci: f64) {
        self.base.set_precision(preci);
        self.my_fix_wire.set_precision(preci);
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L165-169 — SetMinTolerance.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::SetMinTolerance (cxx L165-169).
    pub fn set_min_tolerance(&mut self, mintol: f64) {
        self.base.set_min_tolerance(mintol);
        self.my_fix_wire.base.set_min_tolerance(mintol);
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L173-177 — SetMaxTolerance.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::SetMaxTolerance (cxx L173-177).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.base.set_max_tolerance(maxtol);
        self.my_fix_wire.base.set_max_tolerance(maxtol);
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L181-186 — Init(surf, preci, fwd).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::Init(surf, preci, fwd) (cxx L181-186): starts the
    /// creation of the face; by default it will be FORWARD, or REVERSED if
    /// `fwd` is False.
    pub fn init_surface(
        &mut self,
        brep: &mut BRep,
        surf: &rcad_kernel::geom::Surface3,
        preci: f64,
        fwd: bool,
    ) {
        self.my_status = 0;
        let sas = ShapeAnalysisSurface::new(surf.clone());
        self.init_analyzer_surface(brep, &sas, preci, fwd);
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L190-205 — Init(sas, preci, fwd).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::Init(surf(ShapeAnalysis_Surface), preci, fwd)
    /// (cxx L190-205).
    pub fn init_analyzer_surface(
        &mut self,
        brep: &mut BRep,
        surf: &ShapeAnalysisSurface,
        preci: f64,
        fwd: bool,
    ) {
        self.my_status = 0;
        self.my_surf = Some(sas_dup(surf));
        self.set_precision(preci);
        // OCCT L197-198: BRep_Builder B; B.MakeFace(myFace, mySurf->Surface(),
        // Precision::Confusion()) — the natural-restriction face with no
        // wires.
        let face = brep.add_tface(
            Some(surf.surface().clone()),
            Shape::null(),
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );
        self.my_face = face;
        let my_face = self.my_face.clone();
        self.base.my_shape = my_face;
        self.my_fwd = fwd;
        if !fwd {
            self.my_face.orientation = Orientation::Reversed;
        }
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L209-221 — Init(face).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::Init(face) (cxx L209-221): loads a whole face
    /// already created, with its wires, sense and location.
    pub fn init_face(&mut self, brep: &mut BRep, the_face: &Shape) {
        self.my_status = 0;
        // OCCT L212-217: check if surface is null — no
        // ShapeAnalysis_Surface in that case.
        let a_surface = super::split_tool::brep_tool_surface(brep, the_face);
        if a_surface.is_some() {
            self.my_surf = Some(ShapeAnalysisSurface::new(a_surface.unwrap()));
        } else {
            self.my_surf = None;
        }
        self.my_fwd = the_face.orientation != Orientation::Reversed;
        self.my_face = the_face.clone();
        let my_face = self.my_face.clone();
        self.base.my_shape = my_face;
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L225-235 — Add.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::Add (cxx L225-235): adds a wire to the current
    /// face using BRep_Builder, without taking into account the orientation
    /// of the face (as if the face were FORWARD).
    pub fn add(&mut self, brep: &mut BRep, wire: &Shape) {
        if wire.is_null() {
            return;
        }
        // OCCT L233: fc = myFace.Oriented(TopAbs_FORWARD) — :l2 abv 10 Jan 99.
        let fc = shape_oriented_fwd(&self.my_face);
        builder_add(brep, &fc, wire);
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L345-770 — Perform.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::Perform (cxx L345-770): performs all the fixes,
    /// depending on the modes.
    pub fn perform(&mut self, brep: &mut BRep, the_progress: MessageProgressRange) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        // OCCT L348: myFixWire->SetContext(Context()).
        self.my_fix_wire.base.my_context = self.base.my_context.clone();
        // OCCT L349-353: theAdvFixWire (the rcad value is never null — the
        // handle-null branch is dead).
        let a_init_face = self.my_face.clone();
        // perform first part of fixes on wires
        let mut isfix_reorder = false;
        let mut is_replaced = false;

        // OCCT L362: gka fix in order to avoid lost messages (OCC21771).
        let mut a_map_reordered_wires: ShapeShapeMap = HashMap::new();

        let a_sav_preci = self.base.my_precision;
        if ShapeFixRoot::need_fix(self.my_fix_wire_mode, true) {
            // OCCT L367: SetFace(myFace, mySurf).
            match self.my_surf.as_ref() {
                Some(sas) => {
                    let sas = sas_dup(sas);
                    self.my_fix_wire.set_face_with_surface(brep, &self.my_face, sas)
                }
                // The OCCT null-handle form of SetFace(face, null) — the
                // analyzer surface is cleared (bridge #3).
                None => self.my_fix_wire.set_face(brep, &self.my_face),
            }

            // OCCT L369-374.
            let us_fix_lacking_mode = *self.my_fix_wire.fix_lacking_mode();
            let us_fix_self_intersection_mode = *self.my_fix_wire.fix_self_intersection_mode();
            *self.my_fix_wire.fix_lacking_mode() = 0;
            *self.my_fix_wire.fix_self_intersection_mode() = 0;

            let mut fixed = false;
            let mut s = self.my_face.clone();
            if self.base.my_context.is_some() {
                let f = self.my_face.clone();
                s = self.context_apply(brep, &f);
            }
            // OCCT L382-384.
            let empty_copied = brep.empty_copied(&s);
            let mut tmp_face = empty_copied;
            tmp_face.orientation = Orientation::Forward;

            // OCCT L386-397: skl 29.03.2010 (OCC21623).
            if self.my_auto_correct_precision_mode != 0 {
                let size = least_edge_size(brep, &s);
                let mut newpreci = a_sav_preci.min(size / 2.0);
                newpreci *= 1.00001f64;
                if a_sav_preci > newpreci && newpreci > CONFUSION {
                    self.set_precision(newpreci);
                    self.my_fix_wire.set_precision(newpreci);
                }
            }

            isfix_reorder = false;

            // OCCT L401-455.
            for value in iter_subshapes(brep, &s, false, true) {
                if value.shape_type() != ShapeType::Wire {
                    builder_add(brep, &tmp_face, &value);
                    continue;
                }
                let mut wire = value;
                self.my_fix_wire.load(brep, &wire);
                if self.my_fix_wire.nb_edges() == 0 {
                    if self.my_fix_wire.wire_data().unwrap().nb_nonmanifold_edges() != 0 {
                        builder_add(brep, &tmp_face, &wire);
                    } else {
                        fixed = true;
                        self.my_status |= encode_status(ShapeExtendStatus::Done5);
                    }
                    continue;
                }
                if self.my_fix_wire.perform(brep, the_progress) {
                    // fixed = true;
                    isfix_reorder =
                        self.my_fix_wire.status_reorder(ShapeExtendStatus::Done) || isfix_reorder;
                    fixed = self.my_fix_wire.status_small(ShapeExtendStatus::Done)
                        || self.my_fix_wire.status_connected(ShapeExtendStatus::Done)
                        || self.my_fix_wire.status_edge_curves(ShapeExtendStatus::Done)
                        || self.my_fix_wire.status_notches(ShapeExtendStatus::Done)
                        || self.my_fix_wire.status_fix_tails(ShapeExtendStatus::Done)
                        || self.my_fix_wire.status_degenerated(ShapeExtendStatus::Done)
                        || self.my_fix_wire.status_closed(ShapeExtendStatus::Done);
                    let w = self.my_fix_wire.wire(brep);
                    if fixed {
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &wire, &w);
                        }
                        if self.my_fix_wire.nb_edges() == 0 {
                            self.my_status |= encode_status(ShapeExtendStatus::Done5);
                            continue;
                        }
                    } else if !wire.is_same(&w) {
                        a_map_reordered_wires.insert((wire.ptr_id(), wire.location), w.clone());
                    }

                    wire = w;
                }
                builder_add(brep, &tmp_face, &wire);
            }

            // OCCT L457-462.
            *self.my_fix_wire.fix_lacking_mode() = us_fix_lacking_mode;
            *self.my_fix_wire.fix_self_intersection_mode() = us_fix_self_intersection_mode;
            if !self.my_fwd {
                tmp_face.orientation = Orientation::Reversed;
            }

            // OCCT L464-471.
            if fixed {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.replace(brep, &s, &tmp_face);
                }
                is_replaced = true;
            }
            // OCCT L472-479.
            if fixed || isfix_reorder {
                self.my_face = tmp_face;
                if !self.my_fix_wire.status_reorder(ShapeExtendStatus::Done5) {
                    self.my_status |= encode_status(ShapeExtendStatus::Done1);
                }
            }
        }

        // OCCT L482-483.
        self.my_result = self.my_face.clone();
        let sav_shape = self.my_face.clone(); // gka BUG 6555

        // OCCT L485-489: specific case for conic surfaces.
        if ShapeFixRoot::need_fix(self.my_fix_periodic_degenerated, true) {
            self.fix_periodic_degenerated(brep);
        }

        // OCCT L491-498: fix missing seam.
        if ShapeFixRoot::need_fix(self.my_fix_missing_seam_mode, true) {
            if self.fix_missing_seam(brep) {
                self.my_status |= encode_status(ShapeExtendStatus::Done3);
            }
        }

        // OCCT L500-502: cycle by all possible faces coming from
        // FixMissingSeam; each face is processed as if it was single.
        for exp in crate::shhealing::shape_build::brep_tool::topexp_explorer(
            brep,
            &self.my_result,
            ShapeType::Face,
        ) {
            self.my_face = exp;
            let mut need_check_split_wire = false;

            // OCCT L508: perform second part of fixes on wires.
            if ShapeFixRoot::need_fix(self.my_fix_wire_mode, true) {
                // OCCT L511.
                match self.my_surf.as_ref() {
                    Some(sas) => {
                        let sas = sas_dup(sas);
                        self.my_fix_wire.set_face_with_surface(brep, &self.my_face, sas)
                    }
                    None => self.my_fix_wire.set_face(brep, &self.my_face),
                }

                // OCCT L513-520.
                let us_fix_small_mode = *self.my_fix_wire.fix_small_mode();
                let us_fix_connected_mode = *self.my_fix_wire.fix_connected_mode();
                let us_fix_edge_curves_mode = *self.my_fix_wire.fix_edge_curves_mode();
                let us_fix_degenerated_mode = *self.my_fix_wire.fix_degenerated_mode();
                *self.my_fix_wire.fix_small_mode() = 0;
                *self.my_fix_wire.fix_connected_mode() = 0;
                *self.my_fix_wire.fix_edge_curves_mode() = 0;
                *self.my_fix_wire.fix_degenerated_mode() = 0;

                let mut fixed = false;
                let mut s = self.my_face.clone();
                if self.base.my_context.is_some() {
                    let f = self.my_face.clone();
                    s = self.context_apply(brep, &f);
                }
                // OCCT L528-530.
                let empty_copied = brep.empty_copied(&s);
                let mut tmp_face = empty_copied;
                tmp_face.orientation = Orientation::Forward;
                for value in iter_subshapes(brep, &s, false, true) {
                    if value.shape_type() != ShapeType::Wire {
                        builder_add(brep, &tmp_face, &value);
                        continue;
                    }

                    let mut wire = value;
                    self.my_fix_wire.load(brep, &wire);
                    if self.my_fix_wire.nb_edges() == 0 {
                        if self.my_fix_wire.wire_data().unwrap().nb_nonmanifold_edges() != 0 {
                            builder_add(brep, &tmp_face, &wire);
                        } else {
                            fixed = true;
                            self.my_status |= encode_status(ShapeExtendStatus::Done5);
                        }
                        continue;
                    }
                    if self.my_fix_wire.perform(brep, the_progress) {
                        isfix_reorder = self.my_fix_wire.status_reorder(ShapeExtendStatus::Done);
                        fixed = self.my_fix_wire.status_lacking(ShapeExtendStatus::Done)
                            || self.my_fix_wire.status_self_intersection(ShapeExtendStatus::Done)
                            || self.my_fix_wire.status_notches(ShapeExtendStatus::Done)
                            || self.my_fix_wire.status_fix_tails(ShapeExtendStatus::Done);
                        let w = self.my_fix_wire.wire(brep);
                        if fixed {
                            if let Some(ctx) = self.base.my_context.as_mut() {
                                ctx.replace(brep, &wire, &w);
                            }
                        } else if !wire.is_same(&w) {
                            a_map_reordered_wires
                                .insert((wire.ptr_id(), wire.location), w.clone());
                        }

                        wire = w;
                    }
                    // OCCT L576-579.
                    if self.my_fix_wire.status_removed_segment() {
                        need_check_split_wire = true;
                    }

                    // OCCT L581: fix for loop of wire.
                    let mut a_loop_wires: Vec<Shape> = Vec::new();
                    if ShapeFixRoot::need_fix(self.my_fix_loop_wires_mode, true)
                        && self.fix_loop_wire(brep, &mut a_loop_wires)
                    {
                        if a_loop_wires.len() > 1 {
                            // Wire was split on several wires
                            let wire_for_msg = wire.clone();
                            self.base.send_warning(
                                &wire_for_msg,
                                &MessageMsg::from_key("FixAdvFace.FixLoopWire.MSG0"),
                            );
                        }
                        self.my_status |= encode_status(ShapeExtendStatus::Done7);
                        fixed = true;
                        for k in 1..=(a_loop_wires.len() as i32) {
                            builder_add(brep, &tmp_face, &a_loop_wires[(k - 1) as usize]);
                        }
                    } else {
                        builder_add(brep, &tmp_face, &wire);
                    }
                }

                // OCCT L605-608.
                *self.my_fix_wire.fix_small_mode() = us_fix_small_mode;
                *self.my_fix_wire.fix_connected_mode() = us_fix_connected_mode;
                *self.my_fix_wire.fix_edge_curves_mode() = us_fix_edge_curves_mode;
                *self.my_fix_wire.fix_degenerated_mode() = us_fix_degenerated_mode;

                // OCCT L610-626.
                if fixed {
                    if !self.my_fwd {
                        tmp_face.orientation = Orientation::Reversed;
                    }
                    if !is_replaced
                        && !occt_is_same(&a_init_face, &self.my_result)
                        && self.base.my_context.is_some()
                    {
                        // gka 06.09.04 BUG 6555
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &a_init_face, &sav_shape);
                        }
                    }
                    if self.base.my_context.is_some() {
                        let s2 = self.context_apply(brep, &s);
                        let _ = s2;
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &s, &tmp_face);
                        }
                    }
                    self.my_face = tmp_face;
                    self.my_status |= encode_status(ShapeExtendStatus::Done1);
                }
            }

            // OCCT L629-673.
            if need_check_split_wire {
                // try to split wire - it is needed if some segments were
                // removed in ShapeFix_Wire::FixSelfIntersection()
                let mut s = self.my_face.clone();
                if self.base.my_context.is_some() {
                    let f = self.my_face.clone();
                    s = self.context_apply(brep, &f);
                }
                // OCCT L638-640.
                let empty_copied = brep.empty_copied(&s);
                let mut tmp_face = empty_copied;
                tmp_face.orientation = Orientation::Forward;
                let mut a_wires: Vec<Shape> = Vec::new();
                let mut nbw = 0i32;
                for value in iter_subshapes(brep, &s, false, true) {
                    if value.shape_type() != ShapeType::Wire {
                        builder_add(brep, &tmp_face, &value);
                        continue;
                    }
                    if value.orientation != Orientation::Forward
                        && value.orientation != Orientation::Reversed
                    {
                        builder_add(brep, &tmp_face, &value);
                        continue;
                    }
                    nbw += 1;
                    let wire = value;
                    super::face_b::split_wire(brep, &tmp_face, &wire, &mut a_wires);
                }
                // OCCT L660-672.
                if nbw < a_wires.len() as i32 {
                    for iw in 1..=(a_wires.len() as i32) {
                        builder_add(brep, &tmp_face, &a_wires[(iw - 1) as usize]);
                    }
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &s, &tmp_face);
                    }
                    self.my_status |= encode_status(ShapeExtendStatus::Done8);
                    self.my_face = tmp_face;
                }
            }

            // OCCT L675-686.
            if self.fix_wires_two_coinc_edges(brep) {
                self.my_status |= encode_status(ShapeExtendStatus::Done7);
            }
            if ShapeFixRoot::need_fix(self.my_fix_intersecting_wires_mode, true) {
                if self.fix_intersecting_wires(brep) {
                    self.my_status |= encode_status(ShapeExtendStatus::Done6);
                }
            }

            // OCCT L688-698: fix orientation.
            let mut map_wires: WireListMap = WireListMap::new();
            map_wires.clear();
            if ShapeFixRoot::need_fix(self.my_fix_orientation_mode, true) {
                if self.fix_orientation_map(brep, &mut map_wires) {
                    self.my_status |= encode_status(ShapeExtendStatus::Done2);
                }
            }

            // OCCT L700.
            let my_face = self.my_face.clone();
            super::split_tool::brep_tools_update(brep, &my_face);

            // OCCT L702-708: fix natural bounds.
            let mut need_split = true;
            if self.fix_add_natural_bound(brep) {
                need_split = false;
                self.my_status |= encode_status(ShapeExtendStatus::Done5);
            }

            // OCCT L710-717: split face.
            if ShapeFixRoot::need_fix(self.my_fix_split_face_mode, true)
                && need_split
                && map_wires.len() > 1
            {
                if self.fix_split_face(brep, &map_wires) {
                    self.my_status |= encode_status(ShapeExtendStatus::Done8);
                }
            }
        }

        // OCCT L720-722: return the original preci.
        self.set_precision(a_sav_preci);
        self.my_fix_wire.set_precision(a_sav_preci);

        // OCCT L724-739: cycle by all possible faces coming from
        // FixAddNaturalBound; each face is processed as if it was single.
        for exp in crate::shhealing::shape_build::brep_tool::topexp_explorer(
            brep,
            &self.my_result,
            ShapeType::Face,
        ) {
            self.my_face = exp;

            // fix small-area wires
            if ShapeFixRoot::need_fix(self.my_fix_small_area_wire_mode, false) {
                let is_remove_face = ShapeFixRoot::need_fix(self.my_remove_small_area_face_mode, false);
                if self.fix_small_area_wire(brep, is_remove_face) {
                    self.my_status |= encode_status(ShapeExtendStatus::Done4);
                }
            }
        }

        // OCCT L741-767.
        if self.base.my_context.is_some() {
            if self.status(ShapeExtendStatus::Done)
                && !is_replaced
                && !shapes_is_same(&a_init_face, &sav_shape)
            {
                // gka fix in order to avoid lost messages (OCC21771)
                if !a_map_reordered_wires.is_empty() {
                    for a_cur_w0 in iter_subshapes(brep, &a_init_face, false, true) {
                        let mut a_cur_w = a_cur_w0;
                        while let Some(a_fix_w) =
                            a_map_reordered_wires.get(&(a_cur_w.ptr_id(), a_cur_w.location))
                        {
                            let a_fix_w = a_fix_w.clone();
                            if let Some(ctx) = self.base.my_context.as_mut() {
                                ctx.replace(brep, &a_cur_w, &a_fix_w);
                            }
                            a_cur_w = a_fix_w;
                        }
                    }
                }
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.replace(brep, &a_init_face, &sav_shape);
                }
            }
            let ctx = self.base.my_context.as_mut().unwrap();
            self.my_result = ctx.apply(brep, &a_init_face, ShapeType::Shape); // gka 06.09.04
        } else if !self.status(ShapeExtendStatus::Done) {
            self.my_result = a_init_face;
        }

        // OCCT L769.
        self.status(ShapeExtendStatus::Done)
    }

    // -----------------------------------------------------------------------
    // Internal helpers shared by the impl submodules.
    // -----------------------------------------------------------------------

    /// OCCT `Context()->Apply(shape)` — the root-context application; a None
    /// context returns the shape unchanged (the callers guard with
    /// `Context().IsNull()` exactly as OCCT does).
    pub(crate) fn context_apply(&mut self, brep: &mut BRep, shape: &Shape) -> Shape {
        match self.base.my_context.as_mut() {
            Some(ctx) => ctx.apply(brep, shape, ShapeType::Shape),
            None => shape.clone(),
        }
    }
}

/// OCCT `TopoDS_Shape::Oriented(o)` — a copy with the given orientation.
pub(crate) fn shape_oriented_fwd(s: &Shape) -> Shape {
    let mut r = s.clone();
    r.orientation = Orientation::Forward;
    r
}

/// OCCT `TopoDS_Shape::operator==` — IsEqual (TShape + location +
/// orientation); the brep-free form for the same-pool shapes.
pub(crate) fn occt_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location && a.orientation == b.orientation
}

/// OCCT `TopoDS_Shape::IsSame` — the TShape + location identity.
pub(crate) fn shapes_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT `handle<ShapeAnalysis_Surface>` — the shared-handle copy re-hosted
/// through the landed `init_from_other` (the rcad value model has no shared
/// handles; the analyzer state is duplicated deterministically).
pub(crate) fn sas_dup(other: &ShapeAnalysisSurface) -> ShapeAnalysisSurface {
    let mut s = ShapeAnalysisSurface::new(other.surface().clone());
    s.init_from_other(other);
    s
}
