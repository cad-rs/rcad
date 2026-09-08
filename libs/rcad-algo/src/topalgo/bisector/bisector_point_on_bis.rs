//! OCCT Bisector_PointOnBis — a triplet (param on C1, param on C2, param on
//! the bisector) plus distance/point, 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_PointOnBis.hxx (L27-75) / .cxx (L22-139).

use glam::DVec2;

/// OCCT Bisector_PointOnBis (Bisector_PointOnBis.hxx L27-75).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PointOnBis {
    param1: f64,
    param2: f64,
    param_bis: f64,
    distance: f64,
    infinite: bool,
    point: DVec2,
}

impl Default for PointOnBis {
    /// OCCT Bisector_PointOnBis() (L22-29).
    fn default() -> Self {
        PointOnBis {
            param1: 0.0,
            param2: 0.0,
            param_bis: 0.0,
            distance: 0.0,
            infinite: false,
            point: DVec2::ZERO,
        }
    }
}

impl PointOnBis {
    /// OCCT Bisector_PointOnBis() (L22-29).
    pub fn new() -> Self {
        PointOnBis::default()
    }

    /// OCCT Bisector_PointOnBis(Param1, Param2, ParamBis, Distance, P) (L33-45).
    pub fn new_full(param1: f64, param2: f64, param_bis: f64, distance: f64, p: DVec2) -> Self {
        PointOnBis {
            param1,
            param2,
            param_bis,
            distance,
            infinite: false,
            point: p,
        }
    }

    /// OCCT ParamOnC1(const double Param) (L49-52) — setter.
    pub fn set_param_on_c1(&mut self, param: f64) {
        self.param1 = param;
    }

    /// OCCT ParamOnC2(const double Param) (L56-59) — setter.
    pub fn set_param_on_c2(&mut self, param: f64) {
        self.param2 = param;
    }

    /// OCCT ParamOnBis(const double Param) (L63-66) — setter.
    pub fn set_param_on_bis(&mut self, param: f64) {
        self.param_bis = param;
    }

    /// OCCT Distance(const double Distance) (L70-73) — setter.
    pub fn set_distance(&mut self, distance: f64) {
        self.distance = distance;
    }

    /// OCCT Point(const gp_Pnt2d& P) (L77-80) — setter.
    pub fn set_point(&mut self, p: DVec2) {
        self.point = p;
    }

    /// OCCT IsInfinite(const bool Infinite) (L84-87) — setter.
    pub fn set_infinite(&mut self, infinite: bool) {
        self.infinite = infinite;
    }

    /// OCCT ParamOnC1() const (L91-94) — getter.
    pub fn param_on_c1(&self) -> f64 {
        self.param1
    }

    /// OCCT ParamOnC2() const (L98-101) — getter.
    pub fn param_on_c2(&self) -> f64 {
        self.param2
    }

    /// OCCT ParamOnBis() const (L105-108) — getter.
    pub fn param_on_bis(&self) -> f64 {
        self.param_bis
    }

    /// OCCT Distance() const (L112-115) — getter.
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// OCCT Point() const (L119-122) — getter.
    pub fn point(&self) -> DVec2 {
        self.point
    }

    /// OCCT IsInfinite() const (L126-129) — getter.
    pub fn is_infinite(&self) -> bool {
        self.infinite
    }

    /// OCCT Dump() const (L133-139).
    pub fn dump(&self) {
        println!("Param1    :{}", self.param1);
        println!("Param2    :{}", self.param2);
        println!("Param Bis :{}", self.param_bis);
        println!("Distance  :{}", self.distance);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // OCCT Bisector_PointOnBis — direct-API regression (not run in CI).
    #[test]
    #[ignore]
    fn point_on_bis_defaults_and_accessors() {
        let mut p = PointOnBis::new();
        assert!(!p.is_infinite());
        p.set_param_on_c1(1.0);
        p.set_param_on_bis(3.0);
        assert_eq!(p.param_on_c1(), 1.0);
        assert_eq!(p.param_on_bis(), 3.0);
    }
}
