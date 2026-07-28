//! TKernel core types: precision constants, unit conversion, color management, message/progress.
//!
//! | Module       | OCCT         | Description                       |
//! |--------------|--------------|-----------------------------------|
//! | `precision`  | `Precision`  | Tolerance constants and helpers   |
//! | `units`      | `Units/API`  | Unit conversion and systems       |
//! | `color`      | `Quantity_*` | Color spaces and named colors     |
//! | `message`    | `Message`    | Progress indication and messaging |

pub mod color;
pub mod message;
pub mod precision;
pub mod units;
