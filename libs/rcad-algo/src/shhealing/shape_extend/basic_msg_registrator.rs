//! OCCT ShapeExtend_BasicMsgRegistrator (TKShHealing): `.hxx` L27-58 and
//! `.cxx` L26-51 — abstract class used for attaching messages to objects
//! (e.g. shapes) during Shape Healing.
//!
//! Architecture mapping: the OCCT class hierarchy
//! (BasicMsgRegistrator <- MsgRegistrator, virtual `Send` dispatch) maps to
//! a Rust trait with default (empty) method bodies; consumers hold
//! `&mut dyn BasicMsgRegistrator` where OCCT holds
//! `handle(ShapeExtend_BasicMsgRegistrator)`.
//!
//! GAP carrier: `Message_Msg` (TKMessage) lives in [`super::msg`];
//! `handle(Standard_Transient)` maps to `Option<u64>` (the transient
//! pointer identity; `None` is the null handle).

use super::msg::MessageMsg;
use rcad_kernel::core::message::MessageGravity;

/// OCCT ShapeExtend_BasicMsgRegistrator (ShapeExtend_BasicMsgRegistrator.hxx
/// L36-58).
pub trait BasicMsgRegistrator {
    /// OCCT ShapeExtend_BasicMsgRegistrator() — empty constructor
    /// (ShapeExtend_BasicMsgRegistrator.cxx L26).
    fn new() -> Self
    where
        Self: Sized;

    /// OCCT Send(handle(Standard_Transient), Message_Msg, Message_Gravity)
    /// (ShapeExtend_BasicMsgRegistrator.cxx L30-34): sends a message to be
    /// attached to the object; the base implementation is empty.
    fn send_transient(
        &mut self,
        _object: Option<u64>,
        _message: &MessageMsg,
        _gravity: MessageGravity,
    ) {
    }

    /// OCCT Send(TopoDS_Shape, Message_Msg, Message_Gravity)
    /// (ShapeExtend_BasicMsgRegistrator.cxx L38-42): sends a message to be
    /// attached to the shape; the base implementation is empty.
    fn send_shape(&mut self, _shape: &rcad_kernel::topo::topods::Shape, _message: &MessageMsg, _gravity: MessageGravity) {
    }

    /// OCCT Send(Message_Msg, Message_Gravity)
    /// (ShapeExtend_BasicMsgRegistrator.cxx L46-51): calls Send with a Null
    /// Transient.
    fn send(&mut self, message: &MessageMsg, gravity: MessageGravity) {
        let dummy: Option<u64> = None;
        self.send_transient(dummy, message, gravity);
    }
}
