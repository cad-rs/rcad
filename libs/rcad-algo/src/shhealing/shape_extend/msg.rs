//! GAP carrier for OCCT `Message_Msg` (TKMessage, `Message_Msg.hxx`).
//!
//! The ShapeExtend registrators store and pass messages by value
//! (`NCollection_List<Message_Msg>`); the message machinery itself
//! (parameter substitution, `Message_MsgFile` lookup by key) belongs to the
//! untranslated TKMessage package.  The carrier keeps the two identifying
//! fields the healing stack relies on (the original message key and the
//! description string) and is clonable like the OCCT handle-semantics copy.

/// OCCT Message_Msg (Message_Msg.hxx L40-70) — GAP carrier: a message
/// identified by its original key, with its description string.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MessageMsg {
    /// OCCT `myOriginalMessageName` — the message key (e.g.
    /// "ShapeFix.FixSmallSolid.MSG0").
    pub original_message_name: String,
    /// OCCT `myDescriptionString` — the message text.
    pub description_string: String,
}

impl MessageMsg {
    /// OCCT Message_Msg() — the empty (NoMsg) message.
    pub fn new() -> Self {
        MessageMsg {
            original_message_name: String::new(),
            description_string: String::new(),
        }
    }

    /// OCCT Message_Msg(theMsgDescr) — a message from its descriptor key.
    pub fn from_key(the_msg_descr: &str) -> Self {
        MessageMsg {
            original_message_name: the_msg_descr.to_string(),
            description_string: String::new(),
        }
    }

    /// OCCT Message_Msg::OriginalMessageName() — the message key.
    pub fn original_message_name(&self) -> &str {
        &self.original_message_name
    }

    /// OCCT Message_Msg::GetDescriptionString() — the message text.
    pub fn description_string(&self) -> &str {
        &self.description_string
    }
}

impl Default for MessageMsg {
    fn default() -> Self {
        Self::new()
    }
}
