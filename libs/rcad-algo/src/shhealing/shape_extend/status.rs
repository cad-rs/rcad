//! OCCT ShapeExtend_Status + ShapeExtend::EncodeStatus/DecodeStatus
//! (TKShHealing ShapeExtend package).
//!
//! 1:1 translation of `ShapeExtend_Status.hxx` (L37-64) and the two static
//! methods in `ShapeExtend.cxx`. The status is a bit field shared by the whole
//! ShapeFix/ShapeAnalysis/ShapeBuild healing stack.

/// OCCT ShapeExtend_Status (ShapeExtend_Status.hxx L37-64).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShapeExtendStatus {
    /// Nothing done, everything OK.
    Ok,
    /// Something was done, case 1.
    Done1,
    /// Something was done, case 2.
    Done2,
    /// Something was done, case 3.
    Done3,
    /// Something was done, case 4.
    Done4,
    /// Something was done, case 5.
    Done5,
    /// Something was done, case 6.
    Done6,
    /// Something was done, case 7.
    Done7,
    /// Something was done, case 8.
    Done8,
    /// Something was done (any of DONE#).
    Done,
    /// The method failed, case 1.
    Fail1,
    /// The method failed, case 2.
    Fail2,
    /// The method failed, case 3.
    Fail3,
    /// The method failed, case 4.
    Fail4,
    /// The method failed, case 5.
    Fail5,
    /// The method failed, case 6.
    Fail6,
    /// The method failed, case 7.
    Fail7,
    /// The method failed, case 8.
    Fail8,
    /// The method failed (any of FAIL# occurred).
    Fail,
}

/// OCCT ShapeExtend::EncodeStatus (ShapeExtend.cxx L38-80): enumeration to a
/// bit flag. OK encodes to 0x0000; DONE# to bit #; DONE to the whole low byte;
/// FAIL# to bit # raised by one byte; FAIL to the whole high byte.
pub fn encode_status(status: ShapeExtendStatus) -> i32 {
    match status {
        ShapeExtendStatus::Ok => 0x0000,
        ShapeExtendStatus::Done1 => 0x0001,
        ShapeExtendStatus::Done2 => 0x0002,
        ShapeExtendStatus::Done3 => 0x0004,
        ShapeExtendStatus::Done4 => 0x0008,
        ShapeExtendStatus::Done5 => 0x0010,
        ShapeExtendStatus::Done6 => 0x0020,
        ShapeExtendStatus::Done7 => 0x0040,
        ShapeExtendStatus::Done8 => 0x0080,
        ShapeExtendStatus::Done => 0x00ff,
        ShapeExtendStatus::Fail1 => 0x0100,
        ShapeExtendStatus::Fail2 => 0x0200,
        ShapeExtendStatus::Fail3 => 0x0400,
        ShapeExtendStatus::Fail4 => 0x0800,
        ShapeExtendStatus::Fail5 => 0x1000,
        ShapeExtendStatus::Fail6 => 0x2000,
        ShapeExtendStatus::Fail7 => 0x4000,
        ShapeExtendStatus::Fail8 => 0x8000,
        ShapeExtendStatus::Fail => 0xff00,
    }
}

/// OCCT ShapeExtend::DecodeStatus (ShapeExtend.cxx L82-90): OK is true only
/// when the flag is completely empty; any other status is a bit test.
pub fn decode_status(flag: i32, status: ShapeExtendStatus) -> bool {
    if status == ShapeExtendStatus::Ok {
        return flag == 0;
    }
    (flag & encode_status(status)) != 0
}
