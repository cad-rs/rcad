//! OCCT GeomFill_Filling (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Filling.hxx L33-53 + GeomFill_Filling.cxx L26-60, the base part
//! shared by the Stretch / Coons / Curved concrete fillings.
//!
//! Mapping: `NCollection_Array2<gp_Pnt>` -> `Vec<Vec<DVec3>>` indexed
//! `[u - 1][v - 1]` inside `1..=n` loops; the `myWeights` null handle ->
//! an empty `Vec` until a rational Init fills it.

use glam::DVec3;

/// OCCT GeomFill_Filling base part.
#[derive(Debug, Clone, Default)]
pub struct FillingBase {
    pub(crate) is_rational: bool,
    pub(crate) poles: Vec<Vec<DVec3>>,
    pub(crate) weights: Vec<Vec<f64>>,
}

impl FillingBase {
    /// OCCT GeomFill_Filling() — IsRational = false, no poles yet.
    pub fn new() -> Self {
        FillingBase {
            is_rational: false,
            poles: Vec::new(),
            weights: Vec::new(),
        }
    }

    /// OCCT NbUPoles() == myPoles->ColLength().
    pub fn nb_u_poles(&self) -> usize {
        self.poles.len()
    }

    /// OCCT NbVPoles() == myPoles->RowLength().
    pub fn nb_v_poles(&self) -> usize {
        self.poles[0].len()
    }

    /// OCCT isRational().
    pub fn is_rational(&self) -> bool {
        self.is_rational
    }

    /// OCCT Poles(Poles) — the pole grid `[u][v]`.
    pub fn poles(&self) -> &Vec<Vec<DVec3>> {
        &self.poles
    }

    /// OCCT Weights(Weights) — the weight grid `[u][v]` (empty when the
    /// filling is not rational).
    pub fn weights(&self) -> &Vec<Vec<f64>> {
        &self.weights
    }
}
