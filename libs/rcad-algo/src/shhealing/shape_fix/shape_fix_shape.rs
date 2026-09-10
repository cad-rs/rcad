//! 1:1 translation of OCCT `ShapeFix_Shape` — the top-level fixing driver
//! (`TKShHealing/ShapeFix/ShapeFix_Shape.hxx` L17-155 + `ShapeFix_Shape.lxx`
//! L17-104 + `ShapeFix_Shape.cxx` L1-358, docket row `ShapeFix_Shape` —
//! W3 tranche 4).
//!
//! Function-count equation (OCCT ShapeFix_Shape.cxx members + lxx inlines =
//! rcad shape_fix_shape.rs): `ShapeFix_Shape()` (cxx L39-50) /
//! `ShapeFix_Shape(shape)` (cxx L54-66) / `Init` (cxx L70-79) / `Perform`
//! (cxx L83-303) / `SameParameter` (cxx L307-312, protected) / `Shape`
//! (cxx L316-319) / `SetMsgRegistrator` (cxx L323-327) / `SetPrecision`
//! (cxx L331-335) / `SetMinTolerance` (cxx L339-343) / `SetMaxTolerance`
//! (cxx L347-351) / `Status` (cxx L355-358) — 11 members = 11 rcad methods;
//! the lxx inlines `FixSolidTool` (L24-27) / `FixShellTool` (L31-34) /
//! `FixFaceTool` (L38-41) / `FixWireTool` (L45-48) / `FixEdgeTool` (L52-55) /
//! `FixSolidMode` (L59-62) / `FixFreeShellMode` (L66-69) / `FixFreeFaceMode`
//! (L73-76) / `FixFreeWireMode` (L80-83) / `FixSameParameterMode` (L87-90) /
//! `FixVertexPositionMode` (L94-97) / `FixVertexTolMode` (L101-104) —
//! 12 = 12.
//!
//! Architecture bridges:
//! 1. `ShapeFix_Root` inheritance — the rcad composition field `base` (the
//!    ShapeFix_Shell precedent).
//! 2. `occ::handle<T>` tools — rcad value fields; the tool accessors return
//!    `&mut` and are never null, so the OCCT handle-null branches are dead.
//! 3. `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` — a
//!    `HashSet<(u64, u32)>` keyed by (TShape pointer, location index) (the
//!    shape_fix.rs `shape_key` bridge; the null location is index 0).
//! 4. `Message_ProgressScope` — the no-abort bridge (shape_fix.rs): `More()`
//!    is always true, so OCCT's non-aborted control flow is preserved.
//! 5. `myFixSolid` — OCCT `ShapeFix_Solid` is a later W3 docket row; the
//!    [`super::shape_fix_gap_deps::ShapeFixSolidGap`] GAP carrier hosts the
//!    consumed surface (constructor / Init / inline FixShellTool / the Set*
//!    forwarders to the shell tool) with the real `ShapeFixShell` embedded,
//!    so the FixShellTool/FixFaceTool/FixWireTool/FixEdgeTool chain of
//!    ShapeFix_Shape stays real; the carrier `Perform` keeps OCCT's "nothing
//!    fixed" path.

use std::collections::HashSet;

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, ShapeType, tshape_flags};

use crate::shhealing::shape_build::brep_tool::{iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::edge::ShapeFixEdge;
use crate::shhealing::shape_fix::face_a::ShapeFixFace;
use crate::shhealing::shape_fix::root::{MsgRegistratorHandle, ShapeFixRoot};
use crate::shhealing::shape_fix::shell::ShapeFixShell;
use crate::shhealing::shape_fix::wire::ShapeFixWire;

use super::shape_fix::{MessageProgressRange, MessageProgressScope, shape_key};
use super::shape_fix_gap_deps::ShapeFixSolidGap;

/// OCCT ShapeFix_Shape (hxx L44-151): fixing shape in general.
pub struct ShapeFixShape {
    /// The OCCT ShapeFix_Root base subobject (bridge 1; carries the
    /// inherited myShape / myContext / myMsgReg / myPrecision / myMinTol /
    /// myMaxTol).
    pub base: ShapeFixRoot,
    /// OCCT myResult (hxx L140).
    pub(crate) my_result: Shape,
    /// OCCT myFixSolid (hxx L141) — the ShapeFix_Solid GAP carrier (bridge
    /// 5).
    pub(crate) my_fix_solid: ShapeFixSolidGap,
    /// OCCT myMapFixingShape (hxx L142) — the assembly-sharing guard
    /// (bridge 3).
    pub(crate) my_map_fixing_shape: HashSet<(u64, u32)>,
    /// OCCT myFixSolidMode (hxx L143).
    pub(crate) my_fix_solid_mode: i32,
    /// OCCT myFixShellMode (hxx L144).
    pub(crate) my_fix_shell_mode: i32,
    /// OCCT myFixFaceMode (hxx L145).
    pub(crate) my_fix_face_mode: i32,
    /// OCCT myFixWireMode (hxx L146).
    pub(crate) my_fix_wire_mode: i32,
    /// OCCT myFixSameParameterMode (hxx L147).
    pub(crate) my_fix_same_parameter_mode: i32,
    /// OCCT myFixVertexPositionMode (hxx L148).
    pub(crate) my_fix_vertex_position_mode: i32,
    /// OCCT myFixVertexTolMode (hxx L149).
    pub(crate) my_fix_vertex_tol_mode: i32,
    /// OCCT myStatus (hxx L150).
    pub(crate) my_status: i32,
}

impl Default for ShapeFixShape {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixShape {
    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.cxx L39-50 — ShapeFix_Shape().
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::ShapeFix_Shape() (cxx L39-50): empty constructor.
    pub fn new() -> Self {
        ShapeFixShape {
            base: ShapeFixRoot::new(),
            // OCCT L140: myResult — the default-constructed null handle.
            my_result: Shape::null(),
            // OCCT L49: myFixSolid = new ShapeFix_Solid.
            my_fix_solid: ShapeFixSolidGap::new(),
            my_map_fixing_shape: HashSet::new(),
            // OCCT L42-48.
            my_fix_solid_mode: -1,
            my_fix_shell_mode: -1,
            my_fix_face_mode: -1,
            my_fix_wire_mode: -1,
            my_fix_same_parameter_mode: -1,
            my_fix_vertex_position_mode: 0,
            my_fix_vertex_tol_mode: -1,
            // OCCT L41: myStatus = EncodeStatus(ShapeExtend_OK).
            my_status: encode_status(ShapeExtendStatus::Ok),
        }
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.cxx L54-66 — ShapeFix_Shape(shape).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::ShapeFix_Shape(shape) (cxx L54-66): initializes
    /// by a shape.
    pub fn with_shape(shape: &Shape) -> Self {
        let mut this = ShapeFixShape {
            base: ShapeFixRoot::new(),
            // OCCT L140: myResult — the default-constructed null handle.
            my_result: Shape::null(),
            // OCCT L62: myFixSolid = new ShapeFix_Solid.
            my_fix_solid: ShapeFixSolidGap::new(),
            my_map_fixing_shape: HashSet::new(),
            // OCCT L57-61 and L63-64.
            my_fix_solid_mode: -1,
            my_fix_shell_mode: -1,
            my_fix_face_mode: -1,
            my_fix_wire_mode: -1,
            my_fix_same_parameter_mode: -1,
            my_fix_vertex_position_mode: 0,
            my_fix_vertex_tol_mode: -1,
            // OCCT L56: myStatus = EncodeStatus(ShapeExtend_OK).
            my_status: encode_status(ShapeExtendStatus::Ok),
        };
        // OCCT L65: Init(shape).
        this.init(shape);
        this
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.lxx L24-55 — the tool accessors.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::FixSolidTool (lxx L24-27): returns the tool for
    /// fixing solids (the ShapeFix_Solid GAP carrier).
    pub fn fix_solid_tool(&mut self) -> &mut ShapeFixSolidGap {
        &mut self.my_fix_solid
    }

    /// OCCT ShapeFix_Shape::FixShellTool (lxx L31-34): returns
    /// myFixSolid->FixShellTool().
    pub fn fix_shell_tool(&mut self) -> &mut ShapeFixShell {
        self.my_fix_solid.fix_shell_tool()
    }

    /// OCCT ShapeFix_Shape::FixFaceTool (lxx L38-41): returns myFixSolid->
    /// FixShellTool()->FixFaceTool().
    pub fn fix_face_tool(&mut self) -> &mut ShapeFixFace {
        self.fix_shell_tool().fix_face_tool()
    }

    /// OCCT ShapeFix_Shape::FixWireTool (lxx L45-48): returns myFixSolid->
    /// FixShellTool()->FixFaceTool()->FixWireTool().
    pub fn fix_wire_tool(&mut self) -> &mut ShapeFixWire {
        self.fix_face_tool().fix_wire_tool()
    }

    /// OCCT ShapeFix_Shape::FixEdgeTool (lxx L52-55): returns myFixSolid->
    /// FixShellTool()->FixFaceTool()->FixWireTool()->FixEdgeTool().
    pub fn fix_edge_tool(&mut self) -> &mut ShapeFixEdge {
        self.fix_wire_tool().fix_edge_tool()
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.lxx L59-104 — the mode accessors.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::FixSolidMode (lxx L59-62).
    pub fn fix_solid_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_solid_mode
    }

    /// OCCT ShapeFix_Shape::FixFreeShellMode (lxx L66-69).
    pub fn fix_free_shell_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_shell_mode
    }

    /// OCCT ShapeFix_Shape::FixFreeFaceMode (lxx L73-76).
    pub fn fix_free_face_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_face_mode
    }

    /// OCCT ShapeFix_Shape::FixFreeWireMode (lxx L80-83).
    pub fn fix_free_wire_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_wire_mode
    }

    /// OCCT ShapeFix_Shape::FixSameParameterMode (lxx L87-90).
    pub fn fix_same_parameter_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_same_parameter_mode
    }

    /// OCCT ShapeFix_Shape::FixVertexPositionMode (lxx L94-97).
    pub fn fix_vertex_position_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_vertex_position_mode
    }

    /// OCCT ShapeFix_Shape::FixVertexTolMode (lxx L101-104).
    pub fn fix_vertex_tol_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_vertex_tol_mode
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.cxx L70-79 — Init.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::Init (cxx L70-79): initializes by a shape.
    pub fn init(&mut self, shape: &Shape) {
        // OCCT L72: myShape = shape.
        self.base.my_shape = shape.clone();
        // OCCT L73-77.
        if self.base.my_context.is_none() {
            // OCCT L75: SetContext(new ShapeBuild_ReShape).
            self.base.set_context(ShapeBuildReShape::new());
            // OCCT L76: Context()->ModeConsiderLocation() = true.
            let ctx = self.base.my_context.as_mut().unwrap();
            *ctx.mode_consider_location() = true;
        }
        // OCCT L78: myResult = myShape.
        self.my_result = self.base.my_shape.clone();
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.cxx L83-303 — Perform.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::Perform (cxx L83-303): iterates on sub-shapes
    /// and performs fixes.  `the_progress` is the rcad no-abort bridge (the
    /// OCCT default `Message_ProgressRange` argument).
    pub fn perform(&mut self, brep: &mut BRep, the_progress: MessageProgressRange) -> bool {
        // OCCT L85-86: int savFixSmallAreaWireMode = 0 (the initial value is
        // overwritten by the L90 assignment — the rcad form merges the
        // declaration into it); savFixVertexTolMode = myFixVertexTolMode.
        let sav_fix_vertex_tol_mode = self.my_fix_vertex_tol_mode;
        // OCCT L87-95: fft = FixFaceTool(); the rcad value tool is never
        // null (the handle-null branch is dead);
        // savFixSmallAreaWireMode = fft->FixSmallAreaWireMode().
        let mut sav_fix_small_area_wire_mode = *self.fix_face_tool().fix_small_area_wire_mode();
        if sav_fix_small_area_wire_mode == -1
            && self.base.my_shape.shape_type() == ShapeType::Face
        {
            // OCCT L93: fft->FixSmallAreaWireMode() = true (bool assigned to
            // the int& mode).
            *self.fix_face_tool().fix_small_area_wire_mode() = 1;
        }

        // OCCT L97-99: myStatus = EncodeStatus(ShapeExtend_OK); status;
        // TopAbs_ShapeEnum st.
        self.my_status = encode_status(ShapeExtendStatus::Ok);
        let mut status = false;

        // OCCT L101: gka fix for sharing assembly.
        // OCCT L102-104: TopLoc_Location nullLoc, L; L = myShape.Location()
        // — the location is re-assigned to its own value below (no-op in the
        // rcad value model).
        let mut a_shape_null_loc = self.base.my_shape.clone();
        // OCCT L105: aIsRecorded = Context()->IsNewShape(myShape).
        let a_is_recorded = self
            .base
            .my_context
            .as_ref()
            .unwrap()
            .is_new_shape(&self.base.my_shape);
        // OCCT L106: aShapeNullLoc.Location(nullLoc) — the identity location
        // (the rcad identity-location index is 0).
        a_shape_null_loc.location = 0;
        // OCCT L107-113.
        if a_is_recorded || self.my_map_fixing_shape.contains(&shape_key(&a_shape_null_loc)) {
            // OCCT L109: myShape.Location(L, false) — re-assigns its own
            // location (no-op).
            let my_shape = self.base.my_shape.clone();
            // OCCT L110: myResult = Context()->Apply(myShape).
            self.my_result = self.context_apply(brep, &my_shape);
            // OCCT L111-112: status = true; return status.
            status = true;
            return status;
        }
        // OCCT L114: myMapFixingShape.Add(aShapeNullLoc).
        self.my_map_fixing_shape.insert(shape_key(&a_shape_null_loc));
        // OCCT L116: myShape.Location(L, false) — no-op.
        // OCCT L117: S = Context()->Apply(myShape).
        let my_shape = self.base.my_shape.clone();
        let mut s = self.context_apply(brep, &my_shape);
        // OCCT L118-121.
        if ShapeFixRoot::need_fix(self.my_fix_vertex_position_mode, true) {
            // OCCT L120: ShapeFix::FixVertexPosition(S, Precision(),
            // Context()).
            let precision = self.base.precision();
            if self.base.my_context.is_some() {
                let ctx = self.base.my_context.as_mut().unwrap();
                crate::shhealing::shape_fix::shape_fix::fix_vertex_position(
                    brep, &mut s, precision, ctx,
                );
            }
        }

        // OCCT L123: st = S.ShapeType().
        let st = s.shape_type();

        // OCCT L125-128: the progress scope for the fix stages ("Fixing
        // stage", 2 steps) — the no-abort bridge.
        let a_ps = MessageProgressScope;

        // OCCT L130-251: the switch on st.
        match st {
            // OCCT L132-133: case TopAbs_COMPOUND / TopAbs_COMPSOLID.
            ShapeType::Compound | ShapeType::CompSolid => {
                // OCCT L134-137.
                let shape = self.base.my_shape.clone();
                let sav_fix_same_parameter_mode = self.my_fix_same_parameter_mode;
                self.my_fix_same_parameter_mode = 0; // OCCT: = false
                self.my_fix_vertex_tol_mode = 0; // OCCT: = false
                // OCCT L138: aShapesNb = S.NbChildren() — the sub-shape
                // progress-scope step count (the no-abort bridge carries no
                // step counter).
                let a_shapes_nb = brep.nb_children(s.clone());
                let _ = a_shapes_nb;

                // OCCT L140-141: the progress scope for sub-shape fixing
                // ("Fixing sub-shape", aShapesNb steps) — the no-abort
                // bridge.
                let a_ps_sub_shape = MessageProgressScope;
                // OCCT L142-149: for (TopoDS_Iterator anIter(S); anIter.More()
                // && aPSSubShape.More(); anIter.Next()).
                for an_iter_value in iter_subshapes(brep, &s, true, true) {
                    // OCCT L144: myShape = anIter.Value().
                    self.base.my_shape = an_iter_value;
                    // OCCT L145: Perform(aPSSubShape.Next()).
                    if self.perform(brep, the_progress) {
                        status = true;
                    }
                }
                // OCCT L150-153: if (!aPSSubShape.More()) return false; //
                // aborted execution — the no-abort bridge never halts.
                if !a_ps_sub_shape.more() {
                    return false; // aborted execution
                }

                // OCCT L155-157.
                self.my_fix_same_parameter_mode = sav_fix_same_parameter_mode;
                self.my_fix_vertex_tol_mode = sav_fix_vertex_tol_mode;
                self.base.my_shape = shape;
            }
            // OCCT L160-175: case TopAbs_SOLID.
            ShapeType::Solid => {
                // OCCT L161-164: if (!NeedFix(myFixSolidMode)) break.
                if ShapeFixRoot::need_fix(self.my_fix_solid_mode, true) {
                    // OCCT L165: myFixSolid->Init(TopoDS::Solid(S)).
                    self.my_fix_solid.init(brep, &s);
                    // OCCT L166: myFixSolid->SetContext(Context()).
                    self.my_fix_solid.base.my_context = self.base.my_context.clone();

                    // OCCT L168-171: myFixSolid->Perform(aPS.Next()).
                    if self.my_fix_solid.perform(brep, the_progress) {
                        status = true;
                    }

                    // OCCT L173.
                    self.my_status |= encode_status(ShapeExtendStatus::Done4);
                }
            }
            // OCCT L176-192: case TopAbs_SHELL.
            ShapeType::Shell => {
                // OCCT L177-180: if (!NeedFix(myFixShellMode)) break.
                if ShapeFixRoot::need_fix(self.my_fix_shell_mode, true) {
                    // OCCT L181: sfsh = FixShellTool() — the rcad value tool
                    // (never null).
                    let context = self.base.my_context.clone();
                    let sfsh = self.fix_shell_tool();
                    // OCCT L182: sfsh->Init(TopoDS::Shell(S)).
                    sfsh.init(brep, &s);
                    // OCCT L183: sfsh->SetContext(Context()).
                    sfsh.base.my_context = context;

                    // OCCT L185-188: sfsh->Perform(aPS.Next()).
                    if sfsh.perform(brep, the_progress) {
                        status = true;
                    }

                    // OCCT L190.
                    self.my_status |= encode_status(ShapeExtendStatus::Done4);
                }
            }
            // OCCT L193-211: case TopAbs_FACE.
            ShapeType::Face => {
                // OCCT L194-197: if (!NeedFix(myFixFaceMode)) break.
                if ShapeFixRoot::need_fix(self.my_fix_face_mode, true) {
                    // OCCT L198-200.
                    let sav_topo_mode = *self
                        .fix_face_tool()
                        .fix_wire_tool()
                        .modify_topology_mode();
                    *self
                        .fix_face_tool()
                        .fix_wire_tool()
                        .modify_topology_mode() = true;
                    // OCCT L201-202.
                    let context = self.base.my_context.clone();
                    let sff = self.fix_face_tool();
                    sff.init_face(brep, &s);
                    sff.base.my_context = context;

                    // OCCT L204-207: sff->Perform() — the OCCT no-arg
                    // default-progress form.
                    if sff.perform(brep, the_progress) {
                        status = true;
                    }
                    // OCCT L208.
                    *self
                        .fix_face_tool()
                        .fix_wire_tool()
                        .modify_topology_mode() = sav_topo_mode;
                    // OCCT L209.
                    self.my_status |= encode_status(ShapeExtendStatus::Done3);
                }
            }
            // OCCT L212-237: case TopAbs_WIRE.
            ShapeType::Wire => {
                // OCCT L213-216: if (!NeedFix(myFixWireMode)) break.
                if ShapeFixRoot::need_fix(self.my_fix_wire_mode, true) {
                    // OCCT L217-220.
                    let sav_topo_mode;
                    let sav_closed_mode;
                    {
                        let sfw_modes = self.fix_wire_tool();
                        sav_topo_mode = *sfw_modes.modify_topology_mode();
                        sav_closed_mode = *sfw_modes.closed_wire_mode();
                        *sfw_modes.modify_topology_mode() = true;
                    }
                    // OCCT L221-224.
                    {
                        let sfw_modes = self.fix_wire_tool();
                        if !brep.has_flag(s.clone(), tshape_flags::CLOSED) {
                            *sfw_modes.closed_wire_mode() = false;
                        }
                    }
                    // OCCT L225-227.
                    let context = self.base.my_context.clone();
                    let mut wire_replaced: Option<Shape> = None;
                    {
                        let sfw = self.fix_wire_tool();
                        // OCCT L225: sfw->SetFace(TopoDS_Face()) — the null
                        // face.
                        sfw.set_face(brep, &Shape::null());
                        // OCCT L226: sfw->Load(TopoDS::Wire(S)).
                        sfw.load(brep, &s);
                        // OCCT L227: sfw->SetContext(Context()).
                        sfw.base.my_context = context;
                        // OCCT L228-232: sfw->Perform() — the OCCT no-arg
                        // default-progress form.
                        if sfw.perform(brep, the_progress) {
                            status = true;
                            wire_replaced = Some(sfw.wire(brep));
                        }
                    }
                    // OCCT L231: Context()->Replace(S, sfw->Wire()) —
                    // replace for wire only.
                    if let Some(sfww) = wire_replaced {
                        if self.base.my_context.is_some() {
                            let ctx = self.base.my_context.as_mut().unwrap();
                            ctx.replace(brep, &s, &sfww);
                        }
                    }
                    // OCCT L233-234.
                    {
                        let sfw_modes = self.fix_wire_tool();
                        *sfw_modes.modify_topology_mode() = sav_topo_mode;
                        *sfw_modes.closed_wire_mode() = sav_closed_mode;
                    }
                    // OCCT L235.
                    self.my_status |= encode_status(ShapeExtendStatus::Done2);
                }
            }
            // OCCT L238-246: case TopAbs_EDGE.
            ShapeType::Edge => {
                // OCCT L239-240: sfe = FixEdgeTool(); sfe->SetContext(
                // Context()).
                let context = self.base.my_context.clone();
                let fixed = {
                    let sfe = self.fix_edge_tool();
                    sfe.set_context(context.unwrap());
                    // OCCT L241-242: sfe->FixVertexTolerance(TopoDS::Edge(S)).
                    sfe.fix_vertex_tolerance(brep, &s)
                };
                // OCCT L243-244.
                if fixed {
                    self.my_status |= encode_status(ShapeExtendStatus::Done1);
                }
            }
            // OCCT L247-250: case TopAbs_VERTEX / TopAbs_SHAPE / default.
            ShapeType::Vertex | ShapeType::Shape => {}
        }
        // OCCT L252-255: if (!aPS.More()) return false; // aborted execution
        // — the no-abort bridge never halts.
        if !a_ps.more() {
            return false; // aborted execution
        }

        // OCCT L257: myResult = Context()->Apply(S).
        self.my_result = self.context_apply(brep, &s);

        // OCCT L259-266.
        if ShapeFixRoot::need_fix(self.my_fix_same_parameter_mode, true) {
            // OCCT L261: SameParameter(myResult, false, aPS.Next()).
            let my_result = self.my_result.clone();
            self.same_parameter(brep, &my_result, false, the_progress);
            // OCCT L262-265: if (!aPS.More()) return false; // aborted
            // execution — the no-abort bridge never halts.
            if !a_ps.more() {
                return false; // aborted execution
            }
        }
        // OCCT L267-293.
        if ShapeFixRoot::need_fix(self.my_fix_vertex_tol_mode, true) {
            // OCCT L269-274.
            let mut nb_f = 0;
            let an_exp_f = topexp_explorer(brep, &self.my_result, ShapeType::Face);
            let mut an_exp_f_pos = 0usize;
            while an_exp_f_pos < an_exp_f.len() && nb_f <= 1 {
                nb_f += 1;
                an_exp_f_pos += 1;
            }
            if nb_f > 1 {
                // OCCT L277-281: fix for bug 0025455 — for the case when a
                // vertex belongs to different faces it is necessary to check
                // the vertex tolerances after all fixes (e.g. for the case
                // when a cutting edge was performed).
                // OCCT L282: sfe = FixEdgeTool().
                let sfe = self.fix_edge_tool();
                // OCCT L283: anExpF.ReInit() — the second pass over the
                // faces.
                for an_exp_f_cur in &an_exp_f {
                    // OCCT L285: aF = TopoDS::Face(anExpF.Current()).
                    let a_f = an_exp_f_cur.clone();
                    // OCCT L286-288: TopExp_Explorer anExpE(aF,
                    // TopAbs_EDGE).
                    for an_exp_e_cur in topexp_explorer(brep, &a_f, ShapeType::Edge) {
                        // OCCT L289: sfe->FixVertexTolerance(
                        // TopoDS::Edge(anExpE.Current()), aF).
                        sfe.fix_vertex_tolerance_face(brep, &an_exp_e_cur, &a_f);
                    }
                }
            }
        }

        // OCCT L295: myResult = Context()->Apply(myResult).
        let my_result = self.my_result.clone();
        self.my_result = self.context_apply(brep, &my_result);

        // OCCT L297-300: fft->FixSmallAreaWireMode() =
        // savFixSmallAreaWireMode.
        *self.fix_face_tool().fix_small_area_wire_mode() = sav_fix_small_area_wire_mode;

        // OCCT L302: return status.
        status
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.cxx L307-312 — SameParameter (protected).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::SameParameter (cxx L307-312): fixes the same
    /// parameterization problem on the passed shape by updating the
    /// tolerances of the corresponding topological entities.
    pub(crate) fn same_parameter(
        &mut self,
        brep: &mut BRep,
        sh: &Shape,
        enforce: bool,
        the_progress: MessageProgressRange,
    ) {
        // OCCT L311: ShapeFix::SameParameter(sh, enforce, 0.0, theProgress)
        // — the static's theMsgReg argument keeps its null-handle default.
        crate::shhealing::shape_fix::shape_fix::same_parameter(
            brep,
            sh,
            enforce,
            0.0,
            the_progress,
            None,
        );
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Shape.cxx L316-358 — the accessors and the setters.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Shape::Shape (cxx L316-319): returns the resulting
    /// shape.
    pub fn shape_result(&self) -> Shape {
        self.my_result.clone()
    }

    /// OCCT ShapeFix_Shape::SetMsgRegistrator (cxx L323-327): sets the
    /// message registrator.
    pub fn set_msg_registrator(&mut self, msgreg: MsgRegistratorHandle) {
        // OCCT L325: ShapeFix_Root::SetMsgRegistrator(msgreg).
        self.base.set_msg_registrator(msgreg.clone());
        // OCCT L326: myFixSolid->SetMsgRegistrator(msgreg).
        self.my_fix_solid.set_msg_registrator(msgreg);
    }

    /// OCCT ShapeFix_Shape::SetPrecision (cxx L331-335): sets the basic
    /// precision value (also to FixSolidTool).
    pub fn set_precision(&mut self, preci: f64) {
        // OCCT L333.
        self.base.set_precision(preci);
        // OCCT L334.
        self.my_fix_solid.set_precision(preci);
    }

    /// OCCT ShapeFix_Shape::SetMinTolerance (cxx L339-343): sets the minimal
    /// allowed tolerance (also to FixSolidTool).
    pub fn set_min_tolerance(&mut self, mintol: f64) {
        // OCCT L341.
        self.base.set_min_tolerance(mintol);
        // OCCT L342.
        self.my_fix_solid.set_min_tolerance(mintol);
    }

    /// OCCT ShapeFix_Shape::SetMaxTolerance (cxx L347-351): sets the maximal
    /// allowed tolerance (also to FixSolidTool).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        // OCCT L349.
        self.base.set_max_tolerance(maxtol);
        // OCCT L350.
        self.my_fix_solid.set_max_tolerance(maxtol);
    }

    /// OCCT ShapeFix_Shape::Status (cxx L355-358): the status of the last
    /// Fix (a combination of the ShapeExtend_DONE1..DONE6 flags: DONE1 —
    /// some free edges were fixed, DONE2 — some free wires, DONE3 — some
    /// free faces, DONE4 — some free shells/solids).
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    // -----------------------------------------------------------------------
    // Internal helpers shared with the module.
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
