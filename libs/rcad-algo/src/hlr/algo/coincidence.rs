//! OCCT HLRAlgo_Coincidence (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_Coincidence.hxx` (L38-80, fully inline).

use rcad_kernel::topods::State;

/// OCCT HLRAlgo_Coincidence — used in an Interference to store information
/// on the "hiding" edge: 2D data (tangent parameter of the projection at
/// the intersection point) and 3D data (state before / after the face).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coincidence {
    my_fe: i32,
    my_param: f64,
    my_st_bef: State,
    my_st_aft: State,
}

impl Default for Coincidence {
    fn default() -> Self {
        Coincidence {
            my_fe: 0,
            my_param: 0.0,
            my_st_bef: State::In,
            my_st_aft: State::In,
        }
    }
}

impl Coincidence {
    /// OCCT HLRAlgo_Coincidence() — hxx L43-49.
    pub fn new() -> Self {
        Self::default()
    }

    /// OCCT Set2D(FE, Param) — hxx L51-55.
    pub fn set_2d(&mut self, fe: i32, param: f64) {
        self.my_fe = fe;
        self.my_param = param;
    }

    /// OCCT SetState3D(stbef, staft) — hxx L57-61.
    pub fn set_state_3d(&mut self, st_bef: State, st_aft: State) {
        self.my_st_bef = st_bef;
        self.my_st_aft = st_aft;
    }

    /// OCCT Value2D(FE, Param) — hxx L63-67.
    pub fn value_2d(&self) -> (i32, f64) {
        (self.my_fe, self.my_param)
    }

    /// OCCT State3D(stbef, staft) — hxx L69-73.
    pub fn state_3d(&self) -> (State, State) {
        (self.my_st_bef, self.my_st_aft)
    }
}
