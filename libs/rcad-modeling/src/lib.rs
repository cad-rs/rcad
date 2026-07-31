pub mod prim;
pub mod sewing;

pub use prim::primapi::MakeBox;
pub use prim::primapi::MakeCylinder;
pub use prim::primapi::MakeCone;
pub use prim::primapi::MakeSphere;
pub use prim::primapi::MakeTorus;
pub use prim::primapi::{
    box_brep, cone_brep, cylinder_brep, make_box_brep, make_cone_brep, make_cylinder_brep,
    make_sphere_brep, make_torus_brep, sphere_brep, torus_brep,
};
pub use sewing::{SewingResult, sew_shells};

/// Error type returned by the primitive constructors.
/// (Previously defined in the removed `builder` module; kept at the crate root.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BuildError {
    NonFiniteValue(&'static str),
    NonPositiveValue(&'static str),
    ZeroVector(&'static str),
    ParallelVectors(&'static str, &'static str),
    DegenerateGeometry(&'static str),
    InvalidIndex(usize),
}

impl std::fmt::Display for BuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteValue(name) => write!(f, "{name} must be finite"),
            Self::NonPositiveValue(name) => write!(f, "{name} must be > 0"),
            Self::ZeroVector(name) => write!(f, "{name} must be non-zero"),
            Self::ParallelVectors(a, b) => write!(f, "{a} must not be parallel to {b}"),
            Self::DegenerateGeometry(msg) => write!(f, "degenerate geometry: {msg}"),
            Self::InvalidIndex(idx) => write!(f, "invalid index: {idx}"),
        }
    }
}

impl std::error::Error for BuildError {}
