//! 1:1 translation of OCCT `ShapeFix_Root`
//! (`TKShHealing/ShapeFix/ShapeFix_Root.hxx` L17-126 + `ShapeFix_Root.cxx`
//! L1-96 + `ShapeFix_Root.lxx` L1-105, docket row `ShapeFix_Root` — the base
//! class of every ShapeFix fixing operation).
//!
//! Function-count equation (OCCT ShapeFix_Root cxx + lxx = rcad root.rs):
//! `ShapeFix_Root()` / `Set` / `SetContext` / `SetMsgRegistrator` /
//! `SetPrecision` / `SetMinTolerance` / `SetMaxTolerance` / `SendMsg` (cxx
//! L24-95) + `Context` / `MsgRegistrator` / `Precision` / `MinTolerance` /
//! `MaxTolerance` / `LimitTolerance` / `SendMsg(myShape)` / `SendWarning` x2 /
//! `SendFail` x2 / `NeedFix` (lxx L18-104) — 18 OCCT functions = 18 rcad
//! functions.
//!
//! Architecture bridges:
//! - `occ::handle<ShapeExtend_BasicMsgRegistrator>` — a shared
//!   `Rc<RefCell<dyn BasicMsgRegistrator>>` (the handle copy in `Set` and
//!   `SetMsgRegistrator` keeps OCCT's shared-handle semantics; the null
//!   handle is `None`).
//! - `occ::handle<ShapeBuild_ReShape>` — `Option<ShapeBuildReShape>`.
//! - OCCT overloads `SendMsg(shape, msg, grav)` / `SendMsg(msg, grav)` (and
//!   the Warning/Fail pairs) — Rust has no overloading; the myShape forms
//!   carry the `_own` suffix.

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::core::message::MessageGravity;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;

use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::basic_msg_registrator::BasicMsgRegistrator;
use crate::shhealing::shape_extend::msg::MessageMsg;

/// OCCT `handle(ShapeExtend_BasicMsgRegistrator)` — the shared-handle
/// mapping (module doc, bridge 1).
pub type MsgRegistratorHandle = Rc<RefCell<dyn BasicMsgRegistrator>>;

/// OCCT ShapeFix_Root (ShapeFix_Root.hxx L39-126): root class for fixing
/// operations; provides the reshape context, the basic precision value and
/// the limit (minimal and maximal) values for tolerances, and the message
/// registrator.
pub struct ShapeFixRoot {
    /// OCCT myShape (hxx L114, protected).
    pub(crate) my_shape: Shape,
    /// OCCT myContext (hxx L117).
    pub(crate) my_context: Option<ShapeBuildReShape>,
    /// OCCT myMsgReg (hxx L118).
    pub(crate) my_msg_reg: Option<MsgRegistratorHandle>,
    /// OCCT myPrecision (hxx L119).
    pub(crate) my_precision: f64,
    /// OCCT myMinTol (hxx L120).
    pub(crate) my_min_tol: f64,
    /// OCCT myMaxTol (hxx L121).
    pub(crate) my_max_tol: f64,
}

impl Default for ShapeFixRoot {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixRoot {
    /// OCCT ShapeFix_Root::ShapeFix_Root() (cxx L24-28): empty constructor
    /// (no context is created); the tolerances start at Precision::Confusion
    /// and a fresh basic message registrator is allocated.
    pub fn new() -> Self {
        ShapeFixRoot {
            my_shape: Shape::null(),
            my_context: None,
            my_msg_reg: Some(Rc::new(RefCell::new(
                crate::shhealing::shape_extend::basic_msg_registrator::BasicRegistrator::new(),
            ))),
            my_precision: CONFUSION,
            my_min_tol: CONFUSION,
            my_max_tol: CONFUSION,
        }
    }

    /// OCCT ShapeFix_Root::Set (cxx L32-40): copy all fields from another
    /// Root object.
    pub fn set(&mut self, root: &ShapeFixRoot) {
        self.my_context = match &root.my_context {
            Some(c) => Some(c.clone()),
            None => None,
        };
        self.my_msg_reg = root.my_msg_reg.clone();
        self.my_precision = root.my_precision;
        self.my_min_tol = root.my_min_tol;
        self.my_max_tol = root.my_max_tol;
        self.my_shape = root.my_shape.clone();
    }

    /// OCCT ShapeFix_Root::SetContext (cxx L44-47): sets the context.
    pub fn set_context(&mut self, context: ShapeBuildReShape) {
        self.my_context = Some(context);
    }

    /// OCCT ShapeFix_Root::Context (lxx L18-21): returns the context.
    pub fn context(&self) -> Option<&ShapeBuildReShape> {
        self.my_context.as_ref()
    }

    /// OCCT ShapeFix_Root::Context — the mutable form (the OCCT handle is
    /// shared and mutable through it; the rcad value model needs the
    /// explicit accessor).
    pub fn context_mut(&mut self) -> Option<&mut ShapeBuildReShape> {
        self.my_context.as_mut()
    }

    /// OCCT ShapeFix_Root::SetMsgRegistrator (cxx L51-54): sets the message
    /// registrator.
    pub fn set_msg_registrator(&mut self, msgreg: MsgRegistratorHandle) {
        self.my_msg_reg = Some(msgreg);
    }

    /// OCCT ShapeFix_Root::MsgRegistrator (lxx L25-28): returns the message
    /// registrator.
    pub fn msg_registrator(&self) -> &Option<MsgRegistratorHandle> {
        &self.my_msg_reg
    }

    /// OCCT ShapeFix_Root::SetPrecision (cxx L58-69): sets the basic
    /// precision value; the limit tolerances are widened to cover it.
    pub fn set_precision(&mut self, preci: f64) {
        self.my_precision = preci;
        if self.my_max_tol < self.my_precision {
            self.my_max_tol = self.my_precision;
        }
        if self.my_min_tol > self.my_precision {
            self.my_min_tol = self.my_precision;
        }
    }

    /// OCCT ShapeFix_Root::Precision (lxx L32-35): returns the basic
    /// precision value.
    pub fn precision(&self) -> f64 {
        self.my_precision
    }

    /// OCCT ShapeFix_Root::SetMinTolerance (cxx L73-76): sets the minimal
    /// allowed tolerance.
    pub fn set_min_tolerance(&mut self, mintol: f64) {
        self.my_min_tol = mintol;
    }

    /// OCCT ShapeFix_Root::MinTolerance (lxx L39-42): returns the minimal
    /// allowed tolerance.
    pub fn min_tolerance(&self) -> f64 {
        self.my_min_tol
    }

    /// OCCT ShapeFix_Root::SetMaxTolerance (cxx L80-83): sets the maximal
    /// allowed tolerance.
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        self.my_max_tol = maxtol;
    }

    /// OCCT ShapeFix_Root::MaxTolerance (lxx L46-49): returns the maximal
    /// allowed tolerance.
    pub fn max_tolerance(&self) -> f64 {
        self.my_max_tol
    }

    /// OCCT ShapeFix_Root::LimitTolerance (lxx L53-57): the tolerance
    /// limited by [myMinTol, myMaxTol] — only the maximal restriction is
    /// implemented.
    pub fn limit_tolerance(&self, toler: f64) -> f64 {
        if self.my_max_tol < toler {
            self.my_max_tol
        } else {
            toler
        }
    }

    /// OCCT ShapeFix_Root::SendMsg (cxx L87-95): sends a message to be
    /// attached to the shape (calls the corresponding message of the message
    /// registrator).
    pub fn send_msg(&self, shape: &Shape, message: &MessageMsg, gravity: MessageGravity) {
        if let Some(reg) = &self.my_msg_reg {
            reg.borrow_mut().send_shape(shape, message, gravity);
        }
    }

    /// OCCT ShapeFix_Root::SendMsg (lxx L61-64): the myShape form — calls
    /// the previous method (Rust overload suffix, module doc).
    pub fn send_msg_own(&self, message: &MessageMsg, gravity: MessageGravity) {
        self.send_msg(&self.my_shape, message, gravity);
    }

    /// OCCT ShapeFix_Root::SendWarning (lxx L68-71): calls SendMsg with
    /// gravity set to Message_Warning.
    pub fn send_warning(&self, shape: &Shape, message: &MessageMsg) {
        self.send_msg(shape, message, MessageGravity::Alert);
    }

    /// OCCT ShapeFix_Root::SendWarning (lxx L75-78): the myShape form.
    pub fn send_warning_own(&self, message: &MessageMsg) {
        self.send_warning(&self.my_shape, message);
    }

    /// OCCT ShapeFix_Root::SendFail (lxx L82-85): calls SendMsg with gravity
    /// set to Message_Fail.
    pub fn send_fail(&self, shape: &Shape, message: &MessageMsg) {
        self.send_msg(shape, message, MessageGravity::Fail);
    }

    /// OCCT ShapeFix_Root::SendFail (lxx L89-92): the myShape form.
    pub fn send_fail_own(&self, message: &MessageMsg) {
        self.send_fail(&self.my_shape, message);
    }

    /// OCCT ShapeFix_Root::NeedFix (lxx L101-104, static): defines if the
    /// fixing method needs to be called according to its specific flag (when
    /// set) or to the additional criteria (when the flag is default).
    pub fn need_fix(flag: i32, need: bool) -> bool {
        if flag < 0 {
            need
        } else {
            flag > 0
        }
    }
}
