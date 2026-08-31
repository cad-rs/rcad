//! OCCT math::GaussPoints / GaussWeights (TKMath math.cxx L1958-1996).
//!
//! Returns Gauss-Legendre abscissae/weights for `Index` points; the static
//! tables only store the positive half (values repeat mirrored).

use super::gauss_tables::{POINT, WEIGHT};
use super::VecD;

/// OCCT math::GaussPointsMax.
pub fn gauss_points_max() -> usize {
    61
}

/// OCCT math::GaussPoints(Index, GPoint) — fills `GPoint` (1-based) with the
/// Index Gauss points in DECREASING order (as OCCT stores them).
pub fn gauss_points(index: usize, gpoint: &mut VecD) {
    let mut som = 0usize;
    for i in 1..index {
        som += (i + 1) >> 1;
    }
    let ind = (index + 1) >> 1;

    for i in 1..=ind {
        gpoint.set(i, POINT[som + i]);
        if i + ind <= index {
            gpoint.set(i + ind, -POINT[som + i]);
        }
    }
}

/// OCCT math::GaussWeights(Index, GWeight).
pub fn gauss_weights(index: usize, gweight: &mut VecD) {
    let mut som = 0usize;
    for i in 1..index {
        som += (i + 1) >> 1;
    }
    let ind = (index + 1) >> 1;

    for i in 1..=ind {
        gweight.set(i, WEIGHT[som + i]);
        if i + ind <= index {
            gweight.set(i + ind, WEIGHT[som + i]);
        }
    }
}
