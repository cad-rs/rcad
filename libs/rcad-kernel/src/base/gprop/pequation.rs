//! OCCT GProp_PEquation (TKGeomBase GProp package).
//!
//! 1:1 translation of `GProp_PEquation.hxx` (L40-116) +
//! `GProp_PEquation.cxx` (L27-224).  `MathLin::Jacobi` maps to the ported
//! `math_Jacobi` (same cyclic-Jacobi eigen decomposition, OCCT keeps both
//! entry points); `GProp_PGProps` is the sibling `pg_props` module.

use glam::DVec3;

use crate::base::gprop::pg_props::GPropPGProps;
use crate::core::precision::{REAL_FIRST, REAL_LAST};
use crate::geom::Plane;
use crate::math::gp::Lin;
use crate::math::math_jacobi::MathJacobi;

/// OCCT enum class GProp_PEquation::Type (hxx L46-53).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquationType {
    /// Not yet computed.
    None,
    /// All points are coincident within tolerance.
    Point,
    /// Points are collinear within tolerance.
    Line,
    /// Points are coplanar within tolerance.
    Plane,
    /// Points span 3D space.
    Space,
}

/// OCCT GProp_PEquation — analyzes a collection of 3D points to decide
/// whether they are coincident, collinear, coplanar, or span 3D space.
#[derive(Debug, Clone, Copy)]
pub struct PEquation {
    my_type: EquationType,
    /// Type-specific point (mean point / box corner).
    my_g: DVec3,
    /// Plane normal / line direction / box edge 1.
    my_v1: DVec3,
    /// Box edge 2 (Space type only).
    my_v2: DVec3,
    /// Box edge 3 (Space type only).
    my_v3: DVec3,
    /// Centre of mass, always valid.
    my_barycentre: DVec3,
    /// Principal axis unit vectors.
    my_axes: [DVec3; 3],
    /// Extent along each principal axis.
    my_extents: [f64; 3],
}

impl PEquation {
    /// OCCT GProp_PEquation(thePnts, theTol) — cxx L27-177.
    pub fn new(the_pnts: &[DVec3], the_tol: f64) -> Self {
        let a_props = GPropPGProps::from_points(the_pnts);
        let my_barycentre = a_props.centre_of_mass();
        let my_g = my_barycentre;

        let an_inertia = a_props.matrix_of_inertia();

        // math_Matrix aMat(1, 3, 1, 3) filled from the inertia matrix.
        let mut a_mat = crate::math::MatD::new(3, 3);
        for i in 1..=3 {
            for j in 1..=3 {
                a_mat.set(i, j, an_inertia[i - 1][j - 1]);
            }
        }

        // const MathUtils::EigenResult anEigen = MathLin::Jacobi(aMat, true)
        // — the ported math_Jacobi performs the same cyclic Jacobi solve
        // with eigenvalue ordering.
        let an_eigen = MathJacobi::new(&a_mat);
        let mut eq = PEquation {
            my_type: EquationType::None,
            my_g,
            my_v1: DVec3::ZERO,
            my_v2: DVec3::ZERO,
            my_v3: DVec3::ZERO,
            my_barycentre,
            my_axes: [DVec3::X, DVec3::Y, DVec3::Z],
            my_extents: [0.0; 3],
        };
        if !an_eigen.is_done() {
            // OCCT L47-54: identity axes, zero extents, Type::None.
            eq.my_axes = [DVec3::X, DVec3::Y, DVec3::Z];
            eq.my_extents = [0.0; 3];
            return eq;
        }

        let a_vecs = an_eigen.vectors();
        // Eigenvectors are the COLUMNS of aVecs (cxx L57-60), 1-based.
        eq.my_axes[0] = DVec3::new(a_vecs.get(1, 1), a_vecs.get(2, 1), a_vecs.get(3, 1));
        eq.my_axes[1] = DVec3::new(a_vecs.get(1, 2), a_vecs.get(2, 2), a_vecs.get(3, 2));
        eq.my_axes[2] = DVec3::new(a_vecs.get(1, 3), a_vecs.get(2, 3), a_vecs.get(3, 3));

        let (a_xg, a_yg, a_zg) = (my_g.x, my_g.y, my_g.z);

        let mut a_max1 = REAL_FIRST;
        let mut a_min1 = REAL_LAST;
        let mut a_max2 = REAL_FIRST;
        let mut a_min2 = REAL_LAST;
        let mut a_max3 = REAL_FIRST;
        let mut a_min3 = REAL_LAST;

        for p in the_pnts {
            let a_dx = p.x - a_xg;
            let a_dy = p.y - a_yg;
            let a_dz = p.z - a_zg;

            let a_d1 = a_dx * eq.my_axes[0].x + a_dy * eq.my_axes[0].y + a_dz * eq.my_axes[0].z;
            if a_d1 > a_max1 {
                a_max1 = a_d1;
            }
            if a_d1 < a_min1 {
                a_min1 = a_d1;
            }

            let a_d2 = a_dx * eq.my_axes[1].x + a_dy * eq.my_axes[1].y + a_dz * eq.my_axes[1].z;
            if a_d2 > a_max2 {
                a_max2 = a_d2;
            }
            if a_d2 < a_min2 {
                a_min2 = a_d2;
            }

            let a_d3 = a_dx * eq.my_axes[2].x + a_dy * eq.my_axes[2].y + a_dz * eq.my_axes[2].z;
            if a_d3 > a_max3 {
                a_max3 = a_d3;
            }
            if a_d3 < a_min3 {
                a_min3 = a_d3;
            }
        }

        eq.my_extents[0] = a_max1 - a_min1;
        eq.my_extents[1] = a_max2 - a_min2;
        eq.my_extents[2] = a_max3 - a_min3;

        let mut a_dimension = 3;
        let mut a_dim_code = 0;
        if eq.my_extents[0].abs() <= the_tol {
            a_dimension -= 1;
            a_dim_code = 1;
        }
        if eq.my_extents[1].abs() <= the_tol {
            a_dimension -= 1;
            a_dim_code = 2 * (a_dim_code + 1);
        }
        if eq.my_extents[2].abs() <= the_tol {
            a_dimension -= 1;
            a_dim_code = 3 * (a_dim_code + 1);
        }

        match a_dimension {
            0 => {
                eq.my_type = EquationType::Point;
            }
            1 => {
                eq.my_type = EquationType::Line;
                if a_dim_code == 4 {
                    eq.my_v1 = eq.my_axes[2];
                } else if a_dim_code == 6 {
                    eq.my_v1 = eq.my_axes[1];
                } else {
                    eq.my_v1 = eq.my_axes[0];
                }
            }
            2 => {
                eq.my_type = EquationType::Plane;
                if a_dim_code == 1 {
                    eq.my_v1 = eq.my_axes[0];
                } else if a_dim_code == 2 {
                    eq.my_v1 = eq.my_axes[1];
                } else {
                    eq.my_v1 = eq.my_axes[2];
                }
            }
            3 => {
                eq.my_type = EquationType::Space;
                eq.my_g = eq.my_g
                    + a_min1 * eq.my_axes[0]
                    + a_min2 * eq.my_axes[1]
                    + a_min3 * eq.my_axes[2];
                eq.my_v1 = eq.my_extents[0] * eq.my_axes[0];
                eq.my_v2 = eq.my_extents[1] * eq.my_axes[1];
                eq.my_v3 = eq.my_extents[2] * eq.my_axes[2];
            }
            _ => {}
        }
        eq
    }

    /// OCCT GetType() — hxx L61.
    pub fn get_type(&self) -> EquationType {
        self.my_type
    }

    /// OCCT IsPlanar() — hxx L64.
    pub fn is_planar(&self) -> bool {
        self.my_type == EquationType::Plane
    }

    /// OCCT IsLinear() — hxx L67.
    pub fn is_linear(&self) -> bool {
        self.my_type == EquationType::Line
    }

    /// OCCT IsPoint() — hxx L70.
    pub fn is_point(&self) -> bool {
        self.my_type == EquationType::Point
    }

    /// OCCT IsSpace() — hxx L73.
    pub fn is_space(&self) -> bool {
        self.my_type == EquationType::Space
    }

    /// OCCT Plane() — cxx L181-188 (mean plane through myG with normal
    /// myV1).
    pub fn plane(&self) -> Plane {
        if !self.is_planar() {
            panic!("GProp_PEquation::Plane: result is not planar");
        }
        Plane::new(self.my_g, self.my_v1)
    }

    /// OCCT Line() — cxx L192-199.
    pub fn line(&self) -> Lin {
        if !self.is_linear() {
            panic!("GProp_PEquation::Line: result is not linear");
        }
        Lin::from_pnt_dir(self.my_g, self.my_v1)
    }

    /// OCCT Point() — cxx L203-210.
    pub fn point(&self) -> DVec3 {
        if !self.is_point() {
            panic!("GProp_PEquation::Point: result is not a point");
        }
        self.my_g
    }

    /// OCCT Box(theP, theV1, theV2, theV3) — cxx L214-224: box corner and
    /// edge vectors along the principal axes.
    pub fn box_parts(&self) -> (DVec3, DVec3, DVec3, DVec3) {
        if !self.is_space() {
            panic!("GProp_PEquation::Box: result is not a 3D space");
        }
        (self.my_g, self.my_v1, self.my_v2, self.my_v3)
    }

    /// OCCT Barycentre() — hxx L96.
    pub fn barycentre(&self) -> DVec3 {
        self.my_barycentre
    }

    /// OCCT PrincipalAxis(theIndex) — hxx L100 (1-based).
    pub fn principal_axis(&self, the_index: usize) -> DVec3 {
        self.my_axes[the_index - 1]
    }

    /// OCCT Extent(theIndex) — hxx L104 (1-based).
    pub fn extent(&self, the_index: usize) -> f64 {
        self.my_extents[the_index - 1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Collinear cloud -> Line type with the direction along X.
    #[test]
    fn collinear_points_classified_as_line() {
        let pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(2.5, 0.0, 0.0),
            DVec3::new(-1.0, 0.0, 0.0),
        ];
        let eq = PEquation::new(&pts, 1e-6);
        assert!(eq.is_linear());
        let l = eq.line();
        assert!((l.dir.x.abs() - 1.0).abs() < 1e-9);
    }

    /// Coincident cloud -> Point type; Planar cloud -> Plane.
    #[test]
    fn point_and_plane_classification() {
        let one = DVec3::new(2.0, -3.0, 5.0);
        let eq_pt = PEquation::new(&[one, one, one], 1e-9);
        assert!(eq_pt.is_point());
        assert_eq!(eq_pt.point(), one);

        let plane_pts = vec![
            DVec3::new(0.0, 0.0, 1.0),
            DVec3::new(4.0, 1.0, 1.0),
            DVec3::new(-2.0, 3.0, 1.0),
            DVec3::new(1.0, -2.0, 1.0),
            DVec3::new(2.0, 2.5, 1.0),
        ];
        let eq_pl = PEquation::new(&plane_pts, 1e-6);
        assert!(eq_pl.is_planar());
        let pl = eq_pl.plane();
        assert!((pl.normal.z.abs() - 1.0).abs() < 1e-6);

        let space_pts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(3.0, 0.0, 0.0),
            DVec3::new(0.0, 4.0, 0.0),
            DVec3::new(0.0, 0.0, 5.0),
        ];
        let eq_sp = PEquation::new(&space_pts, 1e-6);
        assert!(eq_sp.is_space());
        let (_p, v1, v2, v3) = eq_sp.box_parts();
        assert!(v1.length() > 0.0 && v2.length() > 0.0 && v3.length() > 0.0);
    }
}
