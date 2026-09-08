//! OCCT TKTopAlgo/MAT2d — the bisecting locus machinery on a set of
//! 2d lines (engine of BRepFill_OffsetWire), 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/MAT2d/
pub mod mat2d_bi_int;
pub mod mat2d_circuit;
pub mod mat2d_connexion;
pub mod mat2d_cut_curve;
pub mod mat2d_mat2d;
pub mod mat2d_mat2d_b;
pub mod mat2d_mini_path;
pub mod mat2d_tool2d;

pub use mat2d_bi_int::Mat2dBiInt;
pub use mat2d_connexion::{HandleMat2dConnexion, Mat2dConnexion};
pub use mat2d_circuit::Mat2dCircuit;
pub use mat2d_cut_curve::{CurAndInf2d, Mat2dCutCurve};
pub use mat2d_mat2d::Mat2dMat2d;
pub use mat2d_mini_path::{ExtCC2d, Mat2dMiniPath};
pub use mat2d_tool2d::Mat2dTool2d;

// Architecture note (OCCT MAT_Side.hxx, TKTopAlgo/MAT): the MAT package
// translation owns the plain MAT_Side enum; re-exported here because
// MAT2d_Tool2d::Sense / MAT2d_CutCurve consume it.
pub use crate::topalgo::mat::MatSide;

use rcad_kernel::geom::Curve2d;

/// Architecture note (Rust enum vs C++ inheritance):
/// OCCT MAT2d passes the elements of the figure as
/// `occ::handle<Geom2d_Geometry>` whose dynamic type is either a
/// Geom2d_CartesianPoint (a Geom2d_Point) or a Geom2d_Curve (in MAT2d use:
/// a Geom2d_TrimmedCurve).  The `STANDARD_TYPE(Geom2d_CartesianPoint)`
/// dynamic-type tests scattered through MAT2d_Circuit / MAT2d_Tool2d /
/// MAT2d_MiniPath map to `matches!(g, Geom2dGeometry::Point(_))`.
#[derive(Debug, Clone)]
pub enum Geom2dGeometry {
    /// OCCT Geom2d_CartesianPoint — its Pnt2d().
    Point(glam::DVec2),
    /// OCCT Geom2d_Curve (Geom2d_TrimmedCurve in MAT2d use).
    Curve(Curve2d),
}

impl Geom2dGeometry {
    /// OCCT `Type == STANDARD_TYPE(Geom2d_CartesianPoint)`.
    pub fn is_point(&self) -> bool {
        matches!(self, Geom2dGeometry::Point(_))
    }

    /// OCCT `down_cast<Geom2d_Point>(...)->Pnt2d()` — panics on a curve
    /// (OCCT raises Standard_NoSuchObject via the null down-cast).
    pub fn pnt2d(&self) -> glam::DVec2 {
        match self {
            Geom2dGeometry::Point(p) => *p,
            Geom2dGeometry::Curve(_) => panic!("Geom2dGeometry::pnt2d on a curve"),
        }
    }

    /// OCCT `down_cast<Geom2d_Curve>(...)` — panics on a point.
    pub fn curve(&self) -> &Curve2d {
        match self {
            Geom2dGeometry::Curve(c) => c,
            Geom2dGeometry::Point(_) => panic!("Geom2dGeometry::curve on a point"),
        }
    }
}
