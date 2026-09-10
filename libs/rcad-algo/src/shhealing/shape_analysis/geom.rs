//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_Geom`
//! (`ShapeAnalysis_Geom.hxx` L17-60 + `ShapeAnalysis_Geom.cxx` L1-178).
//!
//! Two static utilities: `NearestPlane` (the mean plane of a point cloud via
//! GProp_PEquation) and `PositionTrsf` (the placement transform decoded from
//! a 3x4 transformation coefficient matrix).
//!
//! Architecture bridges:
//! 1. `NCollection_Array1<gp_Pnt>` -> `&[DVec3]` (the 0-based slice; the
//!    OCCT Lower()/Upper() 1-based walk maps to the 0..len index walk).
//! 2. `NCollection_HArray2<double>` -> the local [`HArray2F`] re-host with
//!    the OCCT 1-based (row, col) accessors.
//! 3. `gp_Pln` -> `rcad_kernel::geom::Plane` (`Plane::new(origin, normal)`
//!    is the documented rcad equivalent of `gp_Pln(gp_Pnt, gp_Dir)`).
//! 4. `gp_Trsf` / `gp_Ax3` / `gp_XYZ` -> `rcad_kernel::math::gp::{Trsf,
//!    Ax3}` / `DVec3`; the `gp_GTrsf` coefficient matrix -> the local
//!    [`GTrsf`] re-host (only the Value(i, j) / TranslationPart members the
//!    OCCT body consumes).

use glam::DVec3;
use rcad_kernel::base::gprop::pequation::{EquationType, PEquation};
use rcad_kernel::geom::Plane;
use rcad_kernel::math::gp::{Ax3, Trsf};

// ---------------------------------------------------------------------------
// OCCT re-hosts (architecture bridges #2 and #4).
// ---------------------------------------------------------------------------

/// OCCT NCollection_HArray2<double> — the 1-based 2D array of coefficients.
pub struct HArray2F {
    lower_row: usize,
    #[allow(dead_code)]
    upper_row: usize,
    lower_col: usize,
    upper_col: usize,
    data: Vec<f64>,
}

impl HArray2F {
    /// OCCT NCollection_HArray2<double>(theLowerRow, theUpperRow,
    /// theLowerCol, theUpperCol).
    pub fn new(lower_row: usize, upper_row: usize, lower_col: usize, upper_col: usize) -> Self {
        HArray2F {
            lower_row,
            upper_row,
            lower_col,
            upper_col,
            data: vec![0.0; (upper_row - lower_row + 1) * (upper_col - lower_col + 1)],
        }
    }

    /// OCCT Value(i, j).
    pub fn value(&self, the_row: usize, the_col: usize) -> f64 {
        let r = the_row - self.lower_row;
        let c = the_col - self.lower_col;
        self.data[r * (self.upper_col - self.lower_col + 1) + c]
    }

    /// OCCT SetValue(i, j, V).
    pub fn set_value(&mut self, the_row: usize, the_col: usize, the_value: f64) {
        let r = the_row - self.lower_row;
        let c = the_col - self.lower_col;
        self.data[r * (self.upper_col - self.lower_col + 1) + c] = the_value;
    }
}

/// OCCT gp_GTrsf — the 3x4 general transformation (the matrix part + the
/// translation part), only the members the OCCT PositionTrsf body consumes.
struct GTrsf {
    /// OCCT myMatrix — Value(row, col), row 1..=3, col 1..=4 (the 4th column
    /// is the translation part).
    matrix: [[f64; 4]; 3],
}

impl GTrsf {
    /// OCCT gp_GTrsf() — the identity.
    fn new() -> Self {
        GTrsf {
            matrix: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
        }
    }

    /// OCCT SetValue(i, j, val) (1-based row/col).
    fn set_value(&mut self, the_row: usize, the_col: usize, the_val: f64) {
        self.matrix[the_row - 1][the_col - 1] = the_val;
    }

    /// OCCT Value(i, j) (1-based row/col).
    fn value(&self, the_row: usize, the_col: usize) -> f64 {
        self.matrix[the_row - 1][the_col - 1]
    }

    /// OCCT TranslationPart() — the (1, 4), (2, 4), (3, 4) column.
    fn translation_part(&self) -> DVec3 {
        DVec3::new(
            self.value(1, 4),
            self.value(2, 4),
            self.value(3, 4),
        )
    }
}

/// OCCT gp_Pln::Distance(theP) — the perpendicular distance from the point
/// to the plane (|Distance from Ax3|).
fn plane_distance(the_pln: &Plane, the_p: DVec3) -> f64 {
    (the_p - the_pln.origin).dot(the_pln.normal).abs()
}

// ---------------------------------------------------------------------------
// The class statics.
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_Geom::NearestPlane (cxx L27-96): finds the plane
/// best fitting the point cloud (the GProp_PEquation principal axes); when
/// the smallest extent is significantly smaller than the other two, sets
/// `the_pln` and computes in `the_max_dist` the maximal deviation of the
/// points from the plane (RealFirst() otherwise).
pub fn nearest_plane(
    the_pnts: &[DVec3],
    the_pln: &mut Plane,
    the_max_dist: &mut f64,
) -> bool {
    let an_eq = PEquation::new(the_pnts, 0.0);
    if an_eq.get_type() == EquationType::None {
        return false;
    }

    let a_dev1 = an_eq.extent(1);
    let a_dev2 = an_eq.extent(2);
    let a_dev3 = an_eq.extent(3);

    // Find the axis with the smallest extent - that's the plane normal candidate
    let mut an_axis = if a_dev1 < a_dev2 {
        if a_dev1 < a_dev3 {
            1
        } else {
            3
        }
    } else if a_dev2 < a_dev3 {
        2
    } else {
        3
    };

    // Check if the smallest extent is significantly smaller than the other two
    match an_axis {
        1 => {
            if (2.0 * a_dev1 > a_dev2) || (2.0 * a_dev1 > a_dev3) {
                an_axis = 0;
            } else {
                *the_pln = Plane::new(an_eq.barycentre(), an_eq.principal_axis(1));
            }
        }
        2 => {
            if (2.0 * a_dev2 > a_dev1) || (2.0 * a_dev2 > a_dev3) {
                an_axis = 0;
            } else {
                *the_pln = Plane::new(an_eq.barycentre(), an_eq.principal_axis(2));
            }
        }
        3 => {
            if (2.0 * a_dev3 > a_dev2) || (2.0 * a_dev3 > a_dev1) {
                an_axis = 0;
            } else {
                *the_pln = Plane::new(an_eq.barycentre(), an_eq.principal_axis(3));
            }
        }
        _ => {}
    }

    // OCCT L82: theMaxDist = RealFirst().
    *the_max_dist = -f64::MAX;
    if an_axis != 0 {
        // OCCT L85-93: for (i = thePnts.Lower(); i <= thePnts.Upper(); ++i).
        for &p in the_pnts.iter() {
            let a_d = plane_distance(the_pln, p);
            if *the_max_dist < a_d {
                *the_max_dist = a_d;
            }
        }
    }

    an_axis != 0
}

/// OCCT ShapeAnalysis_Geom::PositionTrsf (cxx L100-178): decodes the
/// placement transform from the 3x4 coefficient matrix; `unit` scales the
/// translation part, `prec` guards the scale/shear checks.  Returns false
/// when the matrix is not a rigid scaling+translation.
pub fn position_trsf(
    coefs: Option<&HArray2F>,
    trsf: &mut Trsf,
    unit: f64,
    prec: f64,
) -> bool {
    let result = true;

    // OCCT L107: trsf = gp_Trsf() (the identity).
    *trsf = Trsf::identity();

    let Some(coefs) = coefs else {
        return true;
    };

    let mut gtrsf = GTrsf::new();
    for i in 1..=3usize {
        for j in 1..=4usize {
            gtrsf.set_value(i, j, coefs.value(i, j));
        }
    }

    let mut v1 = DVec3::new(
        gtrsf.value(1, 1),
        gtrsf.value(2, 1),
        gtrsf.value(3, 1),
    );
    let mut v2 = DVec3::new(
        gtrsf.value(1, 2),
        gtrsf.value(2, 2),
        gtrsf.value(3, 2),
    );
    let mut v3 = DVec3::new(
        gtrsf.value(1, 3),
        gtrsf.value(2, 3),
        gtrsf.value(3, 3),
    );
    let m1 = v1.length();
    let m2 = v2.length();
    let m3 = v3.length();

    if m1 < prec || m2 < prec || m3 < prec {
        return false;
    }
    let mm = (m1 + m2 + m3) / 3.;
    let pmm = prec * mm;
    if (m1 - mm).abs() > pmm || (m2 - mm).abs() > pmm || (m3 - mm).abs() > pmm {
        return false;
    }
    v1 /= m1;
    v2 /= m2;
    v3 /= m3;
    if v1.dot(v2).abs() > prec || v2.dot(v3).abs() > prec || v3.dot(v1).abs() > prec {
        return false;
    }

    // OCCT L148-149: the exact component comparison against the identity.
    if v1.x != 1.
        || v1.y != 0.
        || v1.z != 0.
        || v2.x != 0.
        || v2.y != 1.
        || v2.z != 0.
        || v3.x != 0.
        || v3.y != 0.
        || v3.z != 1.
    {
        // OCCT L151-153: the normalized directions (already unit).
        let d1 = v1;
        let d2 = v2;
        let mut d3 = v3;
        // OCCT L154: gp_Ax3 axes(gp_Pnt(0, 0, 0), d3, d1).
        let mut axes = Ax3::from_pnt_n_vx(DVec3::ZERO, d3, d1);
        // OCCT L155: d3.Cross(d1) — gp_Dir::Cross is in place.
        d3 = d3.cross(d1);
        // OCCT L156-159: if (d3.Dot(d2) < 0) axes.YReverse().
        if d3.dot(d2) < 0. {
            // gp_Ax3::YReverse reverses the Y direction (the handedness
            // flip); the rcad Ax3 has no public YReverse — the same
            // transformation matrix results from the negated Y, because
            // Trsf::set_transformation_ax3 reads only the X/Y/main
            // direction vectors (bridge #4).
            axes.y_direction = -axes.y_direction;
        }
        *trsf = Trsf::set_transformation_ax3(&axes);
    }

    // OCCT L163-166: SetScale(gp_Pnt(0, 0, 0), mm) — the fixed point is the
    // origin, so the pure scale factor is the whole transformation change.
    if (mm - 1.).abs() > prec {
        trsf.set_scale_factor(mm);
    }
    let mut tp = gtrsf.translation_part();
    if unit != 1. {
        tp *= unit;
    }
    if tp.x != 0. || tp.y != 0. || tp.z != 0. {
        trsf.set_translation_part(tp);
    }

    result
}
