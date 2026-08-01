// OCCT CSLib_Class2d (CSLib_Class2d.hxx / .cxx)
// Low-level 2D point-in-polygon classification.
//
// OCCT CSLib_Class2d.cxx L43-98 (init), L141-187 (SiDans), L191-230
// (SiDans_OnMode), L234-275 (internalSiDans), L279-342
// (internalSiDansOuOn).
//
// Determines whether a 2D point lies inside / outside / on the boundary of a
// closed polygon using ray casting. The polygon is internally normalized to
// [0,1] x [0,1] for numerical stability.

use glam::DVec2;
use rcad_kernel::PCONFUSION;

/// OCCT CSLib_Class2d::Result — point-in-polygon classification result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Class2dResult {
    /// Point is strictly inside the polygon.
    Inside = 1,
    /// Point is strictly outside the polygon.
    Outside = -1,
    /// Point is on boundary or classification is uncertain.
    Uncertain = 0,
}

// OCCT CSLib_Class2d.cxx L30-38: transformToNormalized
fn transform_to_normalized(the_u: f64, the_u_min: f64, the_u_range: f64) -> f64 {
    const THE_MIN_RANGE: f64 = 1e-10;
    if the_u_range > THE_MIN_RANGE {
        return (the_u - the_u_min) / the_u_range;
    }
    the_u
}

/// OCCT CSLib_Class2d — low-level 2D point-in-polygon classifier.
pub struct Class2d {
    /// X coordinates (normalized), with the closing point appended at
    /// `points_count` (myPnts2dX).
    px: Vec<f64>,
    /// Y coordinates (normalized), with the closing point appended at
    /// `points_count` (myPnts2dY).
    py: Vec<f64>,
    /// Tolerance in U direction (normalized) (myTolU).
    tol_u: f64,
    /// Tolerance in V direction (normalized) (myTolV).
    tol_v: f64,
    /// Number of polygon vertices (myPointsCount).
    points_count: usize,
    /// Original minimum U bound (myUMin).
    u_min: f64,
    /// Original minimum V bound (myVMin).
    v_min: f64,
    /// Original maximum U bound (myUMax).
    u_max: f64,
    /// Original maximum V bound (myVMax).
    v_max: f64,
}

impl Class2d {
    /// OCCT CSLib_Class2d::init (CSLib_Class2d.cxx L43-98) — construct from
    /// polygon vertices and UV bounds. The polygon is closed automatically
    /// (the first point is repeated at the end).
    pub fn new(
        the_pnts2d: &[DVec2],
        the_tol_u: f64,
        the_tol_v: f64,
        the_u_min: f64,
        the_v_min: f64,
        the_u_max: f64,
        the_v_max: f64,
    ) -> Self {
        let mut c = Class2d {
            px: Vec::new(),
            py: Vec::new(),
            tol_u: the_tol_u,
            tol_v: the_tol_v,
            points_count: 0,
            u_min: the_u_min,
            v_min: the_v_min,
            u_max: the_u_max,
            v_max: the_v_max,
        };
        // Validate input parameters.
        if the_u_max <= the_u_min || the_v_max <= the_v_min || the_pnts2d.len() < 3 {
            c.points_count = 0;
            return c;
        }
        c.points_count = the_pnts2d.len();
        c.tol_u = the_tol_u;
        c.tol_v = the_tol_v;
        // Allocate arrays with one extra element for closing the polygon.
        c.px.resize(c.points_count + 1, 0.0);
        c.py.resize(c.points_count + 1, 0.0);
        let a_du = the_u_max - the_u_min;
        let a_dv = the_v_max - the_v_min;
        // Transform points to normalized coordinates.
        for i in 0..c.points_count {
            let a_p2d = the_pnts2d[i];
            c.px[i] = transform_to_normalized(a_p2d.x, the_u_min, a_du);
            c.py[i] = transform_to_normalized(a_p2d.y, the_v_min, a_dv);
        }
        // Close the polygon by copying first point to last position.
        c.px[c.points_count] = c.px[0];
        c.py[c.points_count] = c.py[0];
        // Normalize tolerances.
        const THE_MIN_RANGE: f64 = 1e-10;
        if a_du > THE_MIN_RANGE {
            c.tol_u /= a_du;
        }
        if a_dv > THE_MIN_RANGE {
            c.tol_v /= a_dv;
        }
        c
    }

    /// OCCT CSLib_Class2d::SiDans (CSLib_Class2d.cxx L141-187) — classify a
    /// point with the construction tolerances.
    pub fn si_dans(&self, the_point: DVec2) -> Class2dResult {
        if self.points_count == 0 {
            return Class2dResult::Uncertain;
        }
        let mut a_x = the_point.x;
        let mut a_y = the_point.y;
        // Compute tolerance in original coordinate space.
        let a_tol_u = self.tol_u * (self.u_max - self.u_min);
        let a_tol_v = self.tol_v * (self.v_max - self.v_min);
        // Quick rejection test for points clearly outside the bounding box.
        if a_x < (self.u_min - a_tol_u)
            || a_x > (self.u_max + a_tol_u)
            || a_y < (self.v_min - a_tol_v)
            || a_y > (self.v_max + a_tol_v)
        {
            return Class2dResult::Outside;
        }
        // Transform to normalized coordinates.
        a_x = transform_to_normalized(a_x, self.u_min, self.u_max - self.u_min);
        a_y = transform_to_normalized(a_y, self.v_min, self.v_max - self.v_min);
        // Perform classification with ON detection.
        let a_result = self.internal_si_dans_ou_on(a_x, a_y);
        if a_result == Class2dResult::Uncertain {
            return Class2dResult::Uncertain; // ON boundary
        }
        // Check corner points with tolerance for boundary detection.
        if self.tol_u > 0.0 || self.tol_v > 0.0 {
            let is_inside = a_result == Class2dResult::Inside;
            if is_inside != self.internal_si_dans(a_x - self.tol_u, a_y - self.tol_v)
                || is_inside != self.internal_si_dans(a_x + self.tol_u, a_y - self.tol_v)
                || is_inside != self.internal_si_dans(a_x - self.tol_u, a_y + self.tol_v)
                || is_inside != self.internal_si_dans(a_x + self.tol_u, a_y + self.tol_v)
            {
                return Class2dResult::Uncertain; // Near boundary
            }
        }
        a_result
    }

    /// OCCT CSLib_Class2d::SiDans_OnMode (CSLib_Class2d.cxx L191-230) —
    /// classify a point with an explicit ON tolerance.
    pub fn si_dans_on_mode(&self, the_point: DVec2, the_tol: f64) -> Class2dResult {
        if self.points_count == 0 {
            return Class2dResult::Uncertain;
        }
        let mut a_x = the_point.x;
        let mut a_y = the_point.y;
        // Quick rejection test.
        if a_x < (self.u_min - the_tol)
            || a_x > (self.u_max + the_tol)
            || a_y < (self.v_min - the_tol)
            || a_y > (self.v_max + the_tol)
        {
            return Class2dResult::Outside;
        }
        // Transform to normalized coordinates.
        a_x = transform_to_normalized(a_x, self.u_min, self.u_max - self.u_min);
        a_y = transform_to_normalized(a_y, self.v_min, self.v_max - self.v_min);
        // Perform classification with ON detection.
        let a_result = self.internal_si_dans_ou_on(a_x, a_y);
        // Check corner points with tolerance.
        if the_tol > 0.0 {
            let is_inside = a_result == Class2dResult::Inside;
            if is_inside != self.internal_si_dans(a_x - the_tol, a_y - the_tol)
                || is_inside != self.internal_si_dans(a_x + the_tol, a_y - the_tol)
                || is_inside != self.internal_si_dans(a_x - the_tol, a_y + the_tol)
                || is_inside != self.internal_si_dans(a_x + the_tol, a_y + the_tol)
            {
                return Class2dResult::Uncertain;
            }
        }
        a_result
    }

    /// OCCT CSLib_Class2d::internalSiDans (CSLib_Class2d.cxx L234-275) —
    /// ray-casting algorithm in normalized coordinates: count edge crossings
    /// with a horizontal ray from (Px, Py) to +infinity.
    fn internal_si_dans(&self, the_px: f64, the_py: f64) -> bool {
        let mut a_nb_crossings = 0usize;
        let mut a_prev_dx = self.px[0] - the_px;
        let mut a_prev_dy = self.py[0] - the_py;
        let mut a_prev_y_is_negative = a_prev_dy < 0.0;
        for a_next_idx in 1..=self.points_count {
            let a_curr_dx = self.px[a_next_idx] - the_px;
            let a_curr_dy = self.py[a_next_idx] - the_py;
            let a_curr_y_is_negative = a_curr_dy < 0.0;
            // Check for edge crossing when Y changes sign.
            if a_curr_y_is_negative != a_prev_y_is_negative {
                if a_prev_dx > 0.0 && a_curr_dx > 0.0 {
                    // Both endpoints are to the right of the test point.
                    a_nb_crossings += 1;
                } else if a_prev_dx > 0.0 || a_curr_dx > 0.0 {
                    // Compute X intersection with horizontal line Y = 0.
                    let a_x_intersect =
                        a_prev_dx - a_prev_dy * (a_curr_dx - a_prev_dx) / (a_curr_dy - a_prev_dy);
                    if a_x_intersect > 0.0 {
                        a_nb_crossings += 1;
                    }
                }
                a_prev_y_is_negative = a_curr_y_is_negative;
            }
            a_prev_dx = a_curr_dx;
            a_prev_dy = a_curr_dy;
        }
        // Odd number of crossings means inside.
        (a_nb_crossings & 1) != 0
    }

    /// OCCT CSLib_Class2d::internalSiDansOuOn (CSLib_Class2d.cxx L279-342) —
    /// ray-casting algorithm with ON detection, in normalized coordinates.
    fn internal_si_dans_ou_on(&self, the_px: f64, the_py: f64) -> Class2dResult {
        let mut a_nb_crossings = 0usize;
        let mut a_prev_dx = self.px[0] - the_px;
        let mut a_prev_dy = self.py[0] - the_py;
        let mut a_prev_y_is_negative = a_prev_dy < 0.0;
        for a_next_idx in 1..=self.points_count {
            let a_prev_idx = a_next_idx - 1;
            let a_curr_dx = self.px[a_next_idx] - the_px;
            let a_curr_dy = self.py[a_next_idx] - the_py;
            // Check if point is very close to current vertex.
            if a_curr_dx < self.tol_u
                && a_curr_dx > -self.tol_u
                && a_curr_dy < self.tol_v
                && a_curr_dy > -self.tol_v
            {
                return Class2dResult::Uncertain; // ON boundary (at vertex)
            }
            // Check if point is ON the edge by computing Y at the test point's X.
            // Skip interpolation for nearly vertical edges to avoid division
            // instability. For vertical edges, ON detection is handled by the
            // tolerance check above.
            let a_edge_dx = self.px[a_next_idx] - self.px[a_prev_idx];
            if (self.px[a_prev_idx] - the_px) * a_curr_dx < 0.0
                && a_edge_dx.abs() > PCONFUSION
            {
                let a_interp_y = self.py[a_next_idx]
                    - (self.py[a_next_idx] - self.py[a_prev_idx]) / a_edge_dx * a_curr_dx;
                let a_delta_y = a_interp_y - the_py;
                if a_delta_y >= -self.tol_v && a_delta_y <= self.tol_v {
                    return Class2dResult::Uncertain; // ON boundary (on edge)
                }
            }
            let a_curr_y_is_negative = a_curr_dy < 0.0;
            if a_curr_y_is_negative != a_prev_y_is_negative {
                if a_prev_dx > 0.0 && a_curr_dx > 0.0 {
                    a_nb_crossings += 1;
                } else if a_prev_dx > 0.0 || a_curr_dx > 0.0 {
                    let a_x_intersect =
                        a_prev_dx - a_prev_dy * (a_curr_dx - a_prev_dx) / (a_curr_dy - a_prev_dy);
                    if a_x_intersect > 0.0 {
                        a_nb_crossings += 1;
                    }
                }
                a_prev_y_is_negative = a_curr_y_is_negative;
            }
            a_prev_dx = a_curr_dx;
            a_prev_dy = a_curr_dy;
        }
        // Odd number of crossings means inside.
        if (a_nb_crossings & 1) != 0 {
            Class2dResult::Inside
        } else {
            Class2dResult::Outside
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unit square polygon (CCW): (0,0),(1,0),(1,1),(0,1).
    fn square() -> Class2d {
        Class2d::new(
            &[
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 0.0),
                DVec2::new(1.0, 1.0),
                DVec2::new(0.0, 1.0),
            ],
            1e-7,
            1e-7,
            0.0,
            0.0,
            1.0,
            1.0,
        )
    }

    #[test]
    fn inside_point() {
        assert_eq!(square().si_dans(DVec2::new(0.5, 0.5)), Class2dResult::Inside);
    }

    #[test]
    fn outside_point() {
        assert_eq!(square().si_dans(DVec2::new(2.0, 2.0)), Class2dResult::Outside);
        // Inside the bounding box but outside the polygon (above the square).
        assert_eq!(square().si_dans(DVec2::new(0.5, 1.5)), Class2dResult::Outside);
    }

    #[test]
    fn on_edge_uncertain() {
        // Point on the bottom edge y=0 is within tolerance → Uncertain.
        assert_eq!(
            square().si_dans(DVec2::new(0.5, 1e-8)),
            Class2dResult::Uncertain
        );
    }

    #[test]
    fn vertex_uncertain() {
        assert_eq!(
            square().si_dans(DVec2::new(1e-8, 1e-8)),
            Class2dResult::Uncertain
        );
    }

    #[test]
    fn degenerate_polygon_uncertain() {
        // < 3 points or zero-size bounds → myPointsCount == 0 → Uncertain.
        let c = Class2d::new(&[DVec2::ZERO, DVec2::new(1.0, 0.0)], 1e-7, 1e-7, 0.0, 0.0, 1.0, 1.0);
        assert_eq!(c.si_dans(DVec2::new(0.5, 0.5)), Class2dResult::Uncertain);
        let c2 = Class2d::new(&square_poly(), 1e-7, 1e-7, 0.0, 0.0, 0.0, 1.0);
        assert_eq!(c2.si_dans(DVec2::new(0.5, 0.5)), Class2dResult::Uncertain);
    }

    #[test]
    fn on_mode() {
        // SiDans_OnMode with a tolerance large enough to reach the boundary.
        let c = square();
        assert_eq!(c.si_dans_on_mode(DVec2::new(0.5, 1e-6), 1e-5), Class2dResult::Uncertain);
    }

    fn square_poly() -> Vec<DVec2> {
        vec![
            DVec2::new(0.0, 0.0),
            DVec2::new(1.0, 0.0),
            DVec2::new(1.0, 1.0),
            DVec2::new(0.0, 1.0),
        ]
    }
}
