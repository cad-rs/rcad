//! OCCT Bisector_PolyBis — polygon of PointOnBis, 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/Bisector/
//!         Bisector_PolyBis.hxx (L29-55) / .cxx (L24-131).

use super::bisector_point_on_bis::PointOnBis;

#[cfg(test)]
use glam::DVec2;

/// OCCT gp::Resolution() (gp.hxx L102).
const GP_RESOLUTION: f64 = 1e-15;

/// OCCT Bisector_PolyBis (Bisector_PolyBis.hxx L29-55).
///
/// OCCT stores `Bisector_PointOnBis thePoints[30]` and uses 1-based indices
/// (`Append` writes `thePoints[++nbPoints]`).  The Rust mirror keeps the same
/// 1-based layout with slot 0 unused, hence the +1 capacity; writing past
/// index 30 is a hard failure here (in C++ it would be out-of-bounds UB).
#[derive(Debug, Clone)]
pub struct PolyBis {
    the_points: [PointOnBis; 31],
    nb_points: i32,
}

impl Default for PolyBis {
    /// OCCT Bisector_PolyBis() (L24-27).
    fn default() -> Self {
        PolyBis {
            the_points: [PointOnBis::new(); 31],
            nb_points: 0,
        }
    }
}

impl PolyBis {
    /// OCCT Bisector_PolyBis() (L24-27).
    pub fn new() -> Self {
        PolyBis::default()
    }

    /// OCCT Append(const Bisector_PointOnBis& P) (L31-35).
    pub fn append(&mut self, p: PointOnBis) {
        self.nb_points += 1;
        assert!(
            (1..=30).contains(&self.nb_points),
            "Bisector_PolyBis::Append: thePoints[30] overflow (OCCT out-of-bounds)"
        );
        self.the_points[self.nb_points as usize] = p;
    }

    /// OCCT Length() const (L39-42).
    pub fn length(&self) -> i32 {
        self.nb_points
    }

    /// OCCT IsEmpty() const (L46-49).
    pub fn is_empty(&self) -> bool {
        self.nb_points == 0
    }

    /// OCCT Value(const int Index) const (L53-56) — 1-based Index.
    pub fn value(&self, index: i32) -> &PointOnBis {
        &self.the_points[index as usize]
    }

    /// OCCT First() const (L60-63).
    pub fn first(&self) -> &PointOnBis {
        &self.the_points[1]
    }

    /// OCCT Last() const (L67-70).
    pub fn last(&self) -> &PointOnBis {
        &self.the_points[self.nb_points as usize]
    }

    /// OCCT Interval(const double U) const (L81-119).
    pub fn interval(&self, u: f64) -> i32 {
        if self.last().param_on_bis() - u < GP_RESOLUTION {
            return self.nb_points - 1;
        }
        let d_u = (self.last().param_on_bis() - self.first().param_on_bis())
            / (self.nb_points - 1) as f64;
        if d_u <= GP_RESOLUTION {
            return 1;
        }

        // C++ `int(...)` truncates toward zero; `as i32` matches.
        let mut int_u = ((u - self.first().param_on_bis()).abs() / d_u) as i32;
        int_u += 1;

        if self.the_points[int_u as usize].param_on_bis() >= u {
            let mut i = int_u;
            while i >= 1 {
                if self.the_points[i as usize].param_on_bis() <= u {
                    int_u = i;
                    break;
                }
                i -= 1;
            }
        } else {
            let mut i = int_u;
            while i <= self.nb_points - 1 {
                if self.the_points[i as usize].param_on_bis() >= u {
                    int_u = i - 1;
                    break;
                }
                i += 1;
            }
        }
        int_u
    }

    /// OCCT Transform(const gp_Trsf2d& T) (L123-131).
    ///
    /// Architecture difference: OCCT gp_Trsf2d maps to `glam::DAffine2`
    /// (see rcad-kernel base/gc/transforms.rs).
    pub fn transform(&mut self, t: &glam::DAffine2) {
        let mut i = 1;
        while i <= self.nb_points {
            let mut p = self.the_points[i as usize].point();
            p = t.transform_point2(p);
            self.the_points[i as usize].set_point(p);
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // OCCT Bisector_PolyBis — direct-API regression.
    #[test]
    fn poly_bis_append_and_interval() {
        let mut poly = PolyBis::new();
        for i in 0..5 {
            let mut p = PointOnBis::new();
            p.set_param_on_bis(i as f64);
            p.set_point(DVec2::new(i as f64, 0.0));
            poly.append(p);
        }
        assert_eq!(poly.length(), 5);
        assert!(!poly.is_empty());
        // OCCT Interval(2.5) on params 0..4 (PolyBis.cxx L81-119): dU = 1,
        // IntU = int(2.5) + 1 = 3; thePoints[3].ParamOnBis() = 2 < 2.5 so
        // the forward walk runs and i=4 (param 3 >= 2.5) sets IntU = 3.
        assert_eq!(poly.interval(2.5), 3);
        // U beyond/at the last parameter returns nbPoints - 1 (L82-85).
        assert_eq!(poly.interval(4.0), 4);
        assert_eq!(poly.interval(5.0), 4);
    }
}
