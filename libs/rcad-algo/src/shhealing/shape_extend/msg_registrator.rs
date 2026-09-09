//! OCCT ShapeExtend_MsgRegistrator (TKShHealing): `.hxx` L42-78, `.cxx`
//! L26-80 and `.lxx` L20-34 — attaches messages to the objects (generic
//! Transient or shape).
//!
//! Architecture mapping:
//! - `NCollection_DataMap(handle(Standard_Transient),
//!   NCollection_List<Message_Msg>)` -> `HashMap<u64, Vec<MessageMsg>>`
//!   (the transient pointer identity is the key; the list maps to a Vec).
//! - `NCollection_DataMap<TopoDS_Shape, ..., TopTools_ShapeMapHasher>`
//!   (IsSame identity: TShape + Location) ->
//!   `HashMap<(u64, u32), Vec<MessageMsg>>` (TShape pointer + location).
//!
//! GAP carrier: `Message_Msg` (TKMessage) lives in [`super::msg`].

use super::basic_msg_registrator::BasicMsgRegistrator;
use super::msg::MessageMsg;
use rcad_kernel::core::message::MessageGravity;
use rcad_kernel::topo::topods::Shape;
use std::collections::HashMap;

/// OCCT ShapeExtend_MsgRegistrator (ShapeExtend_MsgRegistrator.hxx L42-78).
#[derive(Debug, Clone, Default)]
pub struct MsgRegistrator {
    /// OCCT myMapTransient (hxx L74-75).
    my_map_transient: HashMap<u64, Vec<MessageMsg>>,
    /// OCCT myMapShape (hxx L76-77).
    my_map_shape: HashMap<(u64, u32), Vec<MessageMsg>>,
}

impl MsgRegistrator {
    /// OCCT ShapeExtend_MsgRegistrator() (ShapeExtend_MsgRegistrator.cxx
    /// L26-28): creates an object.
    pub fn new() -> Self {
        MsgRegistrator {
            my_map_transient: HashMap::new(),
            my_map_shape: HashMap::new(),
        }
    }

    /// OCCT MapTransient() (ShapeExtend_MsgRegistrator.lxx L20-24): returns
    /// a Map of objects and message list.
    pub fn map_transient(&self) -> &HashMap<u64, Vec<MessageMsg>> {
        &self.my_map_transient
    }

    /// OCCT MapShape() (ShapeExtend_MsgRegistrator.lxx L28-34): returns a
    /// Map of shapes and message list.
    pub fn map_shape(&self) -> &HashMap<(u64, u32), Vec<MessageMsg>> {
        &self.my_map_shape
    }
}

impl BasicMsgRegistrator for MsgRegistrator {
    fn new() -> Self {
        MsgRegistrator::new()
    }

    /// OCCT Send(handle(Standard_Transient), Message_Msg, Message_Gravity)
    /// (ShapeExtend_MsgRegistrator.cxx L32-54): if the object is in the map
    /// then the message is added to the list, otherwise the object is
    /// firstly added to the map.  A null object returns early (the
    /// OCCT_DEBUG warning is compiled out).
    fn send_transient(&mut self, object: Option<u64>, message: &MessageMsg, _gravity: MessageGravity) {
        let Some(object) = object else {
            return;
        };
        match self.my_map_transient.get_mut(&object) {
            Some(list) => list.push(message.clone()),
            None => {
                self.my_map_transient.insert(object, vec![message.clone()]);
            }
        }
    }

    /// OCCT Send(TopoDS_Shape, Message_Msg, Message_Gravity)
    /// (ShapeExtend_MsgRegistrator.cxx L58-80): if the shape is in the map
    /// then the message is added to the list, otherwise the shape is firstly
    /// added to the map.  A null shape returns early (the OCCT_DEBUG warning
    /// is compiled out).
    fn send_shape(&mut self, shape: &Shape, message: &MessageMsg, _gravity: MessageGravity) {
        if shape.is_null() {
            return;
        }
        // TopTools_ShapeMapHasher: IsSame identity (TShape + Location).
        let key = (shape.ptr_id(), shape.location);
        match self.my_map_shape.get_mut(&key) {
            Some(list) => list.push(message.clone()),
            None => {
                self.my_map_shape.insert(key, vec![message.clone()]);
            }
        }
    }
}
