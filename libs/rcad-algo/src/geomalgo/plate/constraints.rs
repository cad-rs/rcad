//! OCCT Plate constraint classes (TKGeomAlgo/Plate) — 1:1 port.
//!
//! Ported files: Plate_GtoCConstraint.cxx, Plate_FreeGtoCConstraint.cxx,
//! Plate_LineConstraint.cxx, Plate_PlaneConstraint.cxx,
//! Plate_GlobalTranslationConstraint.cxx, Plate_SampledCurveConstraint.cxx,
//! Plate_LinearXYZConstraint.cxx, Plate_LinearScalarConstraint.cxx.
//!
//! Architecture mappings: gp_XY -> DVec2, gp_XYZ -> DVec3,
//! NCollection_Array1/HArray1 -> slices/Vec, NCollection_HArray2 -> CoefArray2,
//! Standard_DimensionMismatch -> panic, gp_Lin/gp_Pln flattened to their
//! location + direction invariants (rcad has no double-precision gp_Lin).

use glam::{DVec2, DVec3};

use super::d123::{PlateD1, PlateD2, PlateD3};
use super::pinpoint_constraint::PinpointConstraint;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{MatD, VecD};

/// OCCT static NORMIN (Plate_GtoCConstraint.cxx L26, Plate_FreeGtoCConstraint.cxx L32).
const NORMIN: f64 = 1.0e-10;
/// OCCT static COSMIN (Plate_GtoCConstraint.cxx L27, Plate_FreeGtoCConstraint.cxx L33).
const COSMIN: f64 = 1.0e-2;

/// OCCT NCollection_HArray2<T>(1..ColLen, 1..RowLen) stand-in with 1-based
/// (Row, Col) indexing.
#[derive(Debug, Clone, PartialEq)]
pub struct CoefArray2<T: Clone> {
    col_len: usize,
    row_len: usize,
    data: Vec<T>,
}

impl<T: Clone + Default> CoefArray2<T> {
    pub fn new(col_len: usize, row_len: usize) -> Self {
        CoefArray2 {
            col_len,
            row_len,
            data: vec![T::default(); col_len * row_len],
        }
    }

    /// OCCT HArray2::Init.
    pub fn init(&mut self, value: T) {
        self.data = vec![value; self.col_len * self.row_len];
    }

    pub fn col_len(&self) -> usize {
        self.col_len
    }

    pub fn row_len(&self) -> usize {
        self.row_len
    }

    /// OCCT Value(Row, Col) — 1-based.
    pub fn get(&self, row: usize, col: usize) -> T {
        self.data[(row - 1) * self.row_len + (col - 1)].clone()
    }

    /// OCCT ChangeValue(Row, Col) — 1-based.
    pub fn set(&mut self, row: usize, col: usize, value: T) {
        self.data[(row - 1) * self.row_len + (col - 1)] = value;
    }
}

// ---------------------------------------------------------------------------
// Plate_LinearXYZConstraint (Plate_LinearXYZConstraint.cxx, whole file)
// ---------------------------------------------------------------------------

/// OCCT Plate_LinearXYZConstraint.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearXYZConstraint {
    my_ppc: Vec<PinpointConstraint>,
    my_coef: CoefArray2<f64>,
}

impl Default for LinearXYZConstraint {
    /// OCCT Plate_LinearXYZConstraint() = default (.cxx L21).
    fn default() -> Self {
        LinearXYZConstraint {
            my_ppc: Vec::new(),
            my_coef: CoefArray2::new(0, 0),
        }
    }
}

impl LinearXYZConstraint {
    /// OCCT Plate_LinearXYZConstraint(thePPC, theCoeff) (.cxx L23-39).
    pub fn new(the_ppc: &[PinpointConstraint], the_coeff: &[f64]) -> Self {
        if the_coeff.len() != the_ppc.len() {
            panic!("Standard_DimensionMismatch");
        }
        let mut my_coef = CoefArray2::new(1, the_ppc.len());
        for (i, &c) in the_coeff.iter().enumerate() {
            // OCCT: myCoef->ChangeValue(1, i) = theCoeff(i + Lower - 1)
            my_coef.set(1, i + 1, c);
        }
        LinearXYZConstraint {
            my_ppc: the_ppc.to_vec(),
            my_coef,
        }
    }

    /// OCCT Plate_LinearXYZConstraint(thePPC, theCoeff2d) (.cxx L41-54).
    pub fn new_with_coeffs(the_ppc: &[PinpointConstraint], the_coeff: &CoefArray2<f64>) -> Self {
        if the_coeff.row_len() != the_ppc.len() {
            panic!("Standard_DimensionMismatch");
        }
        LinearXYZConstraint {
            my_ppc: the_ppc.to_vec(),
            my_coef: the_coeff.clone(),
        }
    }

    /// OCCT Plate_LinearXYZConstraint(ColLen, RowLen) (.cxx L56-61).
    pub fn with_size(col_len: usize, row_len: usize) -> Self {
        let mut my_coef = CoefArray2::new(col_len, row_len);
        my_coef.init(0.0);
        LinearXYZConstraint {
            my_ppc: vec![PinpointConstraint::default(); row_len],
            my_coef,
        }
    }

    /// OCCT SetPPC (.cxx L63-66).
    pub fn set_ppc(&mut self, index: usize, value: PinpointConstraint) {
        self.my_ppc[index - 1] = value;
    }

    /// OCCT SetCoeff (.cxx L68-71).
    pub fn set_coeff(&mut self, row: usize, col: usize, value: f64) {
        self.my_coef.set(row, col, value);
    }

    /// OCCT GetPPC() — the constraint array (hxx).
    pub fn get_ppc(&self) -> &[PinpointConstraint] {
        &self.my_ppc
    }

    /// OCCT Coeff() — the coefficient array (hxx).
    pub fn coeff(&self) -> &CoefArray2<f64> {
        &self.my_coef
    }
}

// ---------------------------------------------------------------------------
// Plate_LinearScalarConstraint (Plate_LinearScalarConstraint.cxx, whole file)
// ---------------------------------------------------------------------------

/// OCCT Plate_LinearScalarConstraint.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearScalarConstraint {
    my_ppc: Vec<PinpointConstraint>,
    my_coef: CoefArray2<DVec3>,
}

impl Default for LinearScalarConstraint {
    /// OCCT Plate_LinearScalarConstraint() = default (.cxx L93).
    fn default() -> Self {
        LinearScalarConstraint {
            my_ppc: Vec::new(),
            my_coef: CoefArray2::new(0, 0),
        }
    }
}

impl LinearScalarConstraint {
    /// OCCT Plate_LinearScalarConstraint(PPC1, coeff) (.cxx L95-104).
    pub fn new(ppc1: PinpointConstraint, coeff: DVec3) -> Self {
        let mut my_coef = CoefArray2::new(1, 1);
        my_coef.set(1, 1, coeff);
        LinearScalarConstraint {
            my_ppc: vec![ppc1],
            my_coef,
        }
    }

    /// OCCT Plate_LinearScalarConstraint(thePPC, theCoeff) (.cxx L106-122).
    pub fn new_from_coeffs(the_ppc: &[PinpointConstraint], the_coeff: &[DVec3]) -> Self {
        if the_coeff.len() != the_ppc.len() {
            panic!("Standard_DimensionMismatch");
        }
        let mut my_coef = CoefArray2::new(1, the_ppc.len());
        for (i, &c) in the_coeff.iter().enumerate() {
            my_coef.set(1, i + 1, c);
        }
        LinearScalarConstraint {
            my_ppc: the_ppc.to_vec(),
            my_coef,
        }
    }

    /// OCCT Plate_LinearScalarConstraint(thePPC, theCoeff2d) (.cxx L124-137).
    pub fn new_with_coeffs(the_ppc: &[PinpointConstraint], the_coeff: &CoefArray2<DVec3>) -> Self {
        if the_coeff.row_len() != the_ppc.len() {
            panic!("Standard_DimensionMismatch");
        }
        LinearScalarConstraint {
            my_ppc: the_ppc.to_vec(),
            my_coef: the_coeff.clone(),
        }
    }

    /// OCCT Plate_LinearScalarConstraint(ColLen, RowLen) (.cxx L139-144).
    pub fn with_size(col_len: usize, row_len: usize) -> Self {
        let mut my_coef = CoefArray2::new(col_len, row_len);
        my_coef.init(DVec3::ZERO);
        LinearScalarConstraint {
            my_ppc: vec![PinpointConstraint::default(); row_len],
            my_coef,
        }
    }

    /// OCCT SetPPC (.cxx L146-149).
    pub fn set_ppc(&mut self, index: usize, value: PinpointConstraint) {
        self.my_ppc[index - 1] = value;
    }

    /// OCCT SetCoeff (.cxx L151-154).
    pub fn set_coeff(&mut self, row: usize, col: usize, value: DVec3) {
        self.my_coef.set(row, col, value);
    }

    /// OCCT GetPPC() (hxx).
    pub fn get_ppc(&self) -> &[PinpointConstraint] {
        &self.my_ppc
    }

    /// OCCT Coeff() (hxx).
    pub fn coeff(&self) -> &CoefArray2<DVec3> {
        &self.my_coef
    }
}

// ---------------------------------------------------------------------------
// Plate_GlobalTranslationConstraint (whole file, .cxx L174-197)
// ---------------------------------------------------------------------------

/// OCCT Plate_GlobalTranslationConstraint — all points translate together.
#[derive(Debug, Clone, PartialEq)]
pub struct GlobalTranslationConstraint {
    my_lxyzc: LinearXYZConstraint,
}

impl GlobalTranslationConstraint {
    /// OCCT ctor (.cxx L174-197).
    pub fn new(sof_xy: &[DVec2]) -> Self {
        let mut my_lxyzc = LinearXYZConstraint::with_size(sof_xy.len() - 1, sof_xy.len());
        for (i, &xy) in sof_xy.iter().enumerate() {
            // OCCT: myLXYZC.SetPPC(i, Plate_PinpointConstraint(SOfXY(i), 0,0,0), 0, 0))
            my_lxyzc.set_ppc(
                i + 1,
                PinpointConstraint::new(xy, DVec3::ZERO, 0, 0),
            );
        }
        for i in 1..=sof_xy.len() - 1 {
            my_lxyzc.set_coeff(i, 1, -1.0);
            for j in 2..=sof_xy.len() {
                if j == i + 1 {
                    my_lxyzc.set_coeff(i, j, 1.0);
                } else {
                    my_lxyzc.set_coeff(i, j, 0.0);
                }
            }
        }
        GlobalTranslationConstraint { my_lxyzc }
    }

    /// OCCT LXYZC() (hxx).
    pub fn lxyzc(&self) -> &LinearXYZConstraint {
        &self.my_lxyzc
    }
}

// ---------------------------------------------------------------------------
// Plate_PlaneConstraint (whole file, .cxx L220-230)
// ---------------------------------------------------------------------------

/// OCCT Plate_PlaneConstraint — imposes the plate value on a plane.
#[derive(Debug, Clone, PartialEq)]
pub struct PlaneConstraint {
    my_lsc: LinearScalarConstraint,
}

impl PlaneConstraint {
    /// OCCT ctor (.cxx L220-230).  OCCT signature takes gp_Pln; the plane is
    /// flattened to its location + axis direction invariants (rcad has no
    /// double-precision gp_Pln).
    pub fn new(point2d: DVec2, pln_location: DVec3, pln_axis_direction: DVec3, iu: i32, iv: i32) -> Self {
        let mut my_lsc = LinearScalarConstraint::with_size(1, 1);
        let point = pln_location;
        my_lsc.set_ppc(1, PinpointConstraint::new(point2d, point, iu, iv));
        let dir = pln_axis_direction.normalize();
        my_lsc.set_coeff(1, 1, dir);
        PlaneConstraint { my_lsc }
    }

    /// OCCT LSC() (hxx).
    pub fn lsc(&self) -> &LinearScalarConstraint {
        &self.my_lsc
    }
}

// ---------------------------------------------------------------------------
// Plate_LineConstraint (whole file, .cxx L252-276)
// ---------------------------------------------------------------------------

/// OCCT Plate_LineConstraint — projects the plate value on a line.
#[derive(Debug, Clone, PartialEq)]
pub struct LineConstraint {
    my_lsc: LinearScalarConstraint,
}

impl LineConstraint {
    /// OCCT ctor (.cxx L252-276).  OCCT signature takes gp_Lin; the line is
    /// flattened to its location + direction invariants (rcad has no
    /// double-precision gp_Lin).
    pub fn new(point2d: DVec2, lin_location: DVec3, lin_direction: DVec3, iu: i32, iv: i32) -> Self {
        let mut my_lsc = LinearScalarConstraint::with_size(2, 1);
        let point = lin_location;
        my_lsc.set_ppc(1, PinpointConstraint::new(point2d, point, iu, iv));

        let dir = lin_direction;
        // one builds two directions orthogonal to dir
        let dx = DVec3::new(1.0, 0.0, 0.0);
        let dy = DVec3::new(0.0, 1.0, 0.0);

        let mut d1 = dx.cross(dir);
        let mut d2 = dy.cross(dir);
        if d2.length_squared() > d1.length_squared() {
            d1 = d2;
        }
        d1 = d1 / d1.length();
        d2 = dir.cross(d1);
        d2 = d2 / d2.length();
        my_lsc.set_coeff(1, 1, d1);
        my_lsc.set_coeff(2, 1, d2);
        LineConstraint { my_lsc }
    }

    /// OCCT LSC() (hxx).
    pub fn lsc(&self) -> &LinearScalarConstraint {
        &self.my_lsc
    }
}

// ---------------------------------------------------------------------------
// Plate_SampledCurveConstraint (whole file, .cxx L299-338)
// ---------------------------------------------------------------------------

/// OCCT static inline B0 (.cxx L299-312).
fn b0(t: f64) -> f64 {
    let mut s = t;
    if s < 0.0 {
        s = -s;
    }
    s = 1.0 - s;
    if s < 0.0 {
        s = 0.0;
    }
    s
}

/// OCCT Plate_SampledCurveConstraint — samples a curve constraint.
#[derive(Debug, Clone, PartialEq)]
pub struct SampledCurveConstraint {
    my_lxyzc: LinearXYZConstraint,
}

impl SampledCurveConstraint {
    /// OCCT ctor (.cxx L314-338).  OCCT signature takes
    /// NCollection_Sequence<Plate_PinpointConstraint>; a slice of the same
    /// elements (1-based order preserved) stands in.
    pub fn new(soppc: &[PinpointConstraint], n: usize) -> Self {
        let m = soppc.len();
        if n > m {
            panic!("Standard_DimensionMismatch");
        }
        let mut my_lxyzc = LinearXYZConstraint::with_size(n, m);
        for (index, ppc) in soppc.iter().enumerate() {
            my_lxyzc.set_ppc(index + 1, *ppc);
        }

        let ratio = (n + 1) as f64 / (m + 1) as f64;
        for i in 1..=n {
            for j in 1..=m {
                my_lxyzc.set_coeff(i, j, b0(ratio * j as f64 - i as f64));
            }
        }
        SampledCurveConstraint { my_lxyzc }
    }

    /// OCCT LXYZC() (hxx).
    pub fn lxyzc(&self) -> &LinearXYZConstraint {
        &self.my_lxyzc
    }
}

// ---------------------------------------------------------------------------
// Plate_GtoCConstraint (Plate_GtoCConstraint.cxx, whole file L26-560)
// ---------------------------------------------------------------------------

/// OCCT Plate_GtoCConstraint — G1/G2/G3 constraint between the plate and a
/// target surface.
#[derive(Debug, Clone, PartialEq)]
pub struct GtoCConstraint {
    my_ppc: Vec<PinpointConstraint>, // OCCT myPPC[9]
    my_d1_surf_init: PlateD1,
    pnt2d: DVec2,
    nb_pp_constraints: usize,
}

impl GtoCConstraint {
    /// OCCT G1 ctor (GtoCConstraint.cxx L40-79).
    pub fn new(point2d: DVec2, d1s: &PlateD1, d1t: &PlateD1) -> Self {
        let mut constraint = GtoCConstraint {
            my_ppc: Vec::new(),
            my_d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        };

        let mut normale = d1t.du().cross(d1t.dv());
        // alr le 12/11/96
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        // G1 Constraints

        let mut normale_s = d1s.du().cross(d1s.dv());
        // alr le 12/11/96
        if normale_s.length() < NORMIN {
            return constraint;
        }
        normale_s = normale_s / normale_s.length();
        let cos_normales = normale.dot(normale_s);
        if cos_normales.abs() < COSMIN {
            return constraint;
        }
        let invcos = 1.0 / cos_normales;

        let du = normale_s * (-(normale.dot(d1s.du()))) * invcos;
        let dv = normale_s * (-(normale.dot(d1s.dv()))) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dv, 0, 1));

        constraint.nb_pp_constraints = 2;
        constraint
    }

    /// OCCT G1 ctor with nP (GtoCConstraint.cxx L81-127).
    pub fn new_with_normal(point2d: DVec2, d1s: &PlateD1, d1t: &PlateD1, np: DVec3) -> Self {
        let mut constraint = GtoCConstraint {
            my_ppc: Vec::new(),
            my_d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        };

        let mut normale = d1t.du().cross(d1t.dv());
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        // G1 Constraints

        let mut normale_s = d1s.du().cross(d1s.dv());
        if normale_s.length() < NORMIN {
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        let mut nsp = normale_s - np * (np.dot(normale_s));
        if nsp.length() < NORMIN {
            return constraint;
        }
        nsp = nsp / nsp.length();

        let cos_normales = normale.dot(nsp);
        if cos_normales.abs() < COSMIN {
            return constraint;
        }
        let invcos = 1.0 / cos_normales;

        let du = nsp * (-(normale.dot(d1s.du()))) * invcos;
        let dv = nsp * (-(normale.dot(d1s.dv()))) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dv, 0, 1));

        constraint.nb_pp_constraints = 2;
        constraint
    }

    /// OCCT G1+G2 ctor (GtoCConstraint.cxx L129-214).
    pub fn new_g2(point2d: DVec2, d1s: &PlateD1, d1t: &PlateD1, d2s: &PlateD2, d2t: &PlateD2) -> Self {
        let mut constraint = GtoCConstraint {
            my_ppc: Vec::new(),
            my_d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        };

        let mut normale = d1t.du().cross(d1t.dv());
        // alr le 12/11/96
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        // G1 Constraints
        let mut normale_s = d1s.du().cross(d1s.dv());
        // alr le 12/11/96
        if normale_s.length() < NORMIN {
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        let cos_normales = normale.dot(normale_s);
        if cos_normales.abs() < COSMIN {
            return constraint;
        }
        let invcos = 1.0 / cos_normales;

        let du = normale_s * (-(normale.dot(d1s.du()))) * invcos;
        let dv = normale_s * (-(normale.dot(d1s.dv()))) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dv, 0, 1));
        constraint.nb_pp_constraints = 2;

        // G2 Constraints
        let su = d1s.du() + du;
        let sv = d1s.dv() + dv;

        let mut mat = MatD::new(2, 2);
        mat.set(1, 1, su.dot(d1t.du()));
        mat.set(1, 2, su.dot(d1t.dv()));
        mat.set(2, 1, sv.dot(d1t.du()));
        mat.set(2, 2, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return constraint;
        }

        let mut vec = VecD::new(2);
        vec.set(1, su.dot(su));
        vec.set(2, su.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let a = sol.get(1);
        let b = sol.get(2);

        vec.set(1, sv.dot(su));
        vec.set(2, sv.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let c = sol.get(1);
        let d = sol.get(2);

        let suu = d2t.duu() * (a * a) + d2t.duv() * (2.0 * a * b) + d2t.dvv() * (b * b);
        let suv = d2t.duu() * (a * c) + d2t.duv() * (a * d + b * c) + d2t.dvv() * (b * d);
        let svv = d2t.duu() * (c * c) + d2t.duv() * (2.0 * c * d) + d2t.dvv() * (d * d);

        let duu = normale_s * (normale.dot(suu - d2s.duu())) * invcos;
        let duv = normale_s * (normale.dot(suv - d2s.duv())) * invcos;
        let dvv = normale_s * (normale.dot(svv - d2s.dvv())) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, duu, 2, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duv, 1, 1));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dvv, 0, 2));
        constraint.nb_pp_constraints = 5;
        constraint
    }

    /// OCCT G1+G2 ctor with nP (GtoCConstraint.cxx L216-330).
    pub fn new_g2_with_normal(
        point2d: DVec2,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        np: DVec3,
    ) -> Self {
        let mut constraint = GtoCConstraint {
            my_ppc: Vec::new(),
            my_d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        };

        let mut normale = d1t.du().cross(d1t.dv());
        // alr le 12/11/96
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        let mut normale_s = d1s.du().cross(d1s.dv());
        if normale_s.length() < NORMIN {
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        let mut nsp = normale_s - np * (np.dot(normale_s));
        if nsp.length() < NORMIN {
            return constraint;
        }
        nsp = nsp / nsp.length();

        let cos_normales = normale.dot(nsp);
        if cos_normales.abs() < COSMIN {
            return constraint;
        }
        let invcos = 1.0 / cos_normales;

        let du = nsp * (-(normale.dot(d1s.du()))) * invcos;
        let dv = nsp * (-(normale.dot(d1s.dv()))) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dv, 0, 1));
        constraint.nb_pp_constraints = 2;

        let su = d1s.du() + du;
        let sv = d1s.dv() + dv;

        let mut mat = MatD::new(2, 2);
        mat.set(1, 1, su.dot(d1t.du()));
        mat.set(1, 2, su.dot(d1t.dv()));
        mat.set(2, 1, sv.dot(d1t.du()));
        mat.set(2, 2, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return constraint;
        }

        let mut vec = VecD::new(2);
        vec.set(1, su.dot(su));
        vec.set(2, su.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let a = sol.get(1);
        let b = sol.get(2);

        vec.set(1, sv.dot(su));
        vec.set(2, sv.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let c = sol.get(1);
        let d = sol.get(2);

        let suu = d2t.duu() * (a * a) + d2t.duv() * (2.0 * a * b) + d2t.dvv() * (b * b);
        let suv = d2t.duu() * (a * c) + d2t.duv() * (a * d + b * c) + d2t.dvv() * (b * d);
        let svv = d2t.duu() * (c * c) + d2t.duv() * (2.0 * c * d) + d2t.dvv() * (d * d);

        let duu = nsp * (normale.dot(suu - d2s.duu())) * invcos;
        let duv = nsp * (normale.dot(suv - d2s.duv())) * invcos;
        let dvv = nsp * (normale.dot(svv - d2s.dvv())) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, duu, 2, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duv, 1, 1));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dvv, 0, 2));
        constraint.nb_pp_constraints = 5;
        constraint
    }

    /// OCCT G1+G2+G3 ctor (GtoCConstraint.cxx L332-533) — flat translation.
    pub fn new_g3(
        point2d: DVec2,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        d3s: &PlateD3,
        d3t: &PlateD3,
    ) -> Self {
        let mut constraint = GtoCConstraint {
            my_ppc: Vec::new(),
            my_d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        };

        let mut normale = d1t.du().cross(d1t.dv());
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        let mut normale_s = d1s.du().cross(d1s.dv());
        if normale_s.length() < NORMIN {
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        let cos_normales = normale.dot(normale_s);
        if cos_normales.abs() < COSMIN {
            return constraint;
        }
        let invcos = 1.0 / cos_normales;

        let du = normale_s * (-(normale.dot(d1s.du()))) * invcos;
        let dv = normale_s * (-(normale.dot(d1s.dv()))) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dv, 0, 1));
        constraint.nb_pp_constraints = 2;

        let su = d1s.du() + du;
        let sv = d1s.dv() + dv;

        let mut mat = MatD::new(2, 2);
        mat.set(1, 1, su.dot(d1t.du()));
        mat.set(1, 2, su.dot(d1t.dv()));
        mat.set(2, 1, sv.dot(d1t.du()));
        mat.set(2, 2, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return constraint;
        }

        let mut vec = VecD::new(2);
        vec.set(1, su.dot(su));
        vec.set(2, su.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let a = sol.get(1);
        let b = sol.get(2);

        vec.set(1, sv.dot(su));
        vec.set(2, sv.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let c = sol.get(1);
        let d = sol.get(2);

        let suu = d2t.duu() * (a * a) + d2t.duv() * (2.0 * a * b) + d2t.dvv() * (b * b);
        let suv = d2t.duu() * (a * c) + d2t.duv() * (a * d + b * c) + d2t.dvv() * (b * d);
        let svv = d2t.duu() * (c * c) + d2t.duv() * (2.0 * c * d) + d2t.dvv() * (d * d);

        let duu = normale_s * (normale.dot(suu - d2s.duu())) * invcos;
        let duv = normale_s * (normale.dot(suv - d2s.duv())) * invcos;
        let dvv = normale_s * (normale.dot(svv - d2s.dvv())) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, duu, 2, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duv, 1, 1));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dvv, 0, 2));
        constraint.nb_pp_constraints = 5;

        // G3 Constraints

        vec.set(1, (d2s.duu() + duu - suu).dot(su));
        vec.set(2, (d2s.duu() + duu - suu).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0uu = sol.get(1);
        let b1uu = sol.get(2);

        vec.set(1, (d2s.duv() + duv - suv).dot(su));
        vec.set(2, (d2s.duv() + duv - suv).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0uv = sol.get(1);
        let b1uv = sol.get(2);

        vec.set(1, (d2s.dvv() + dvv - svv).dot(su));
        vec.set(2, (d2s.dvv() + dvv - svv).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0vv = sol.get(1);
        let b1vv = sol.get(2);

        let suuu = d3t.duuu() * (a * a * a)
            + d3t.duuv() * (3.0 * a * a * b)
            + d3t.duvv() * (3.0 * a * b * b)
            + d3t.dvvv() * (b * b * b);
        let suuv = d3t.duuu() * (a * a * c)
            + d3t.duuv() * (a * a * d + 2.0 * a * b * c)
            + d3t.duvv() * (b * b * c + 2.0 * a * b * d)
            + d3t.dvvv() * (b * b * d);
        let suvv = d3t.duuu() * (a * c * c)
            + d3t.duuv() * (b * c * c + 2.0 * a * c * d)
            + d3t.duvv() * (a * d * d + 2.0 * b * c * d)
            + d3t.dvvv() * (b * d * d);
        let svvv = d3t.duuu() * (c * c * c)
            + d3t.duuv() * (3.0 * c * c * d)
            + d3t.duvv() * (3.0 * c * d * d)
            + d3t.dvvv() * (d * d * d);

        // OCCT references: A0u = a, A1u = b, A0v = c, A1v = d.
        let (a0u, a1u, a0v, a1v) = (a, b, c, d);
        let mut suuu = suuu
            + d2t.duu() * (3.0 * a0u * b0uu)
            + d2t.duv() * (3.0 * (a0u * b1uu + a1u * b0uu))
            + d2t.dvv() * (3.0 * a1u * b1uu);
        let mut suuv = suuv
            + d2t.duu() * (2.0 * a0u * b0uv + a0v * b0uu)
            + d2t.duv() * (2.0 * (a0u * b1uv + a1u * b0uv) + a0v * b1uu + a1v * b0uu)
            + d2t.dvv() * (2.0 * a1u * b1uv + a1v * b1uu);
        let mut suvv = suvv
            + d2t.duu() * (a0u * b0vv + 2.0 * a0v * b0uv)
            + d2t.duv() * (2.0 * (a0v * b1uv + a1v * b0uv) + a0u * b1vv + a1u * b0vv)
            + d2t.dvv() * (2.0 * a1v * b1uv + a1u * b1vv);
        let mut svvv = svvv
            + d2t.duu() * (3.0 * a0v * b0vv)
            + d2t.duv() * (3.0 * (a0v * b1vv + a1v * b0vv))
            + d2t.dvv() * (3.0 * a1v * b1vv);
        let _ = (&mut suuu, &mut suuv, &mut suvv, &mut svvv);

        let duuu = normale_s * (normale.dot(suuu - d3s.duuu())) * invcos;
        let duuv = normale_s * (normale.dot(suuv - d3s.duuv())) * invcos;
        let duvv = normale_s * (normale.dot(suvv - d3s.duvv())) * invcos;
        let dvvv = normale_s * (normale.dot(svvv - d3s.dvvv())) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, duuu, 3, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duuv, 2, 1));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duvv, 1, 2));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dvvv, 0, 3));
        constraint.nb_pp_constraints = 9;
        constraint
    }

    /// OCCT G1+G2+G3 ctor with nP (GtoCConstraint.cxx L535-...) — flat
    /// translation; the G1 result routes through nSP instead of normaleS.
    pub fn new_g3_with_normal(
        point2d: DVec2,
        d1s: &PlateD1,
        d1t: &PlateD1,
        d2s: &PlateD2,
        d2t: &PlateD2,
        d3s: &PlateD3,
        d3t: &PlateD3,
        np: DVec3,
    ) -> Self {
        let mut constraint = GtoCConstraint {
            my_ppc: Vec::new(),
            my_d1_surf_init: *d1s,
            pnt2d: point2d,
            nb_pp_constraints: 0,
        };

        let mut normale = d1t.du().cross(d1t.dv());
        // alr le 12/11/96
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        let mut normale_s = d1s.du().cross(d1s.dv());
        if normale_s.length() < NORMIN {
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        let mut nsp = normale_s - np * (np.dot(normale_s));
        if nsp.length() < NORMIN {
            return constraint;
        }
        nsp = nsp / nsp.length();

        let cos_normales = normale.dot(nsp);
        if cos_normales.abs() < COSMIN {
            return constraint;
        }
        let invcos = 1.0 / cos_normales;

        let du = nsp * (-(normale.dot(d1s.du()))) * invcos;
        let dv = nsp * (-(normale.dot(d1s.dv()))) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dv, 0, 1));
        constraint.nb_pp_constraints = 2;

        let su = d1s.du() + du;
        let sv = d1s.dv() + dv;

        let mut mat = MatD::new(2, 2);
        mat.set(1, 1, su.dot(d1t.du()));
        mat.set(1, 2, su.dot(d1t.dv()));
        mat.set(2, 1, sv.dot(d1t.du()));
        mat.set(2, 2, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return constraint;
        }

        let mut vec = VecD::new(2);
        vec.set(1, su.dot(su));
        vec.set(2, su.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let a = sol.get(1);
        let b = sol.get(2);

        vec.set(1, sv.dot(su));
        vec.set(2, sv.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let c = sol.get(1);
        let d = sol.get(2);

        let suu = d2t.duu() * (a * a) + d2t.duv() * (2.0 * a * b) + d2t.dvv() * (b * b);
        let suv = d2t.duu() * (a * c) + d2t.duv() * (a * d + b * c) + d2t.dvv() * (b * d);
        let svv = d2t.duu() * (c * c) + d2t.duv() * (2.0 * c * d) + d2t.dvv() * (d * d);

        let duu = nsp * (normale.dot(suu - d2s.duu())) * invcos;
        let duv = nsp * (normale.dot(suv - d2s.duv())) * invcos;
        let dvv = nsp * (normale.dot(svv - d2s.dvv())) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, duu, 2, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duv, 1, 1));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dvv, 0, 2));
        constraint.nb_pp_constraints = 5;

        // G3 Constraints

        vec.set(1, (d2s.duu() + duu - suu).dot(su));
        vec.set(2, (d2s.duu() + duu - suu).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0uu = sol.get(1);
        let b1uu = sol.get(2);

        vec.set(1, (d2s.duv() + duv - suv).dot(su));
        vec.set(2, (d2s.duv() + duv - suv).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0uv = sol.get(1);
        let b1uv = sol.get(2);

        vec.set(1, (d2s.dvv() + dvv - svv).dot(su));
        vec.set(2, (d2s.dvv() + dvv - svv).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0vv = sol.get(1);
        let b1vv = sol.get(2);

        let suuu = d3t.duuu() * (a * a * a)
            + d3t.duuv() * (3.0 * a * a * b)
            + d3t.duvv() * (3.0 * a * b * b)
            + d3t.dvvv() * (b * b * b);
        let suuv = d3t.duuu() * (a * a * c)
            + d3t.duuv() * (a * a * d + 2.0 * a * b * c)
            + d3t.duvv() * (b * b * c + 2.0 * a * b * d)
            + d3t.dvvv() * (b * b * d);
        let suvv = d3t.duuu() * (a * c * c)
            + d3t.duuv() * (b * c * c + 2.0 * a * c * d)
            + d3t.duvv() * (a * d * d + 2.0 * b * c * d)
            + d3t.dvvv() * (b * d * d);
        let svvv = d3t.duuu() * (c * c * c)
            + d3t.duuv() * (3.0 * c * c * d)
            + d3t.duvv() * (3.0 * c * d * d)
            + d3t.dvvv() * (d * d * d);

        let (a0u, a1u, a0v, a1v) = (a, b, c, d);
        let mut suuu = suuu
            + d2t.duu() * (3.0 * a0u * b0uu)
            + d2t.duv() * (3.0 * (a0u * b1uu + a1u * b0uu))
            + d2t.dvv() * (3.0 * a1u * b1uu);
        let mut suuv = suuv
            + d2t.duu() * (2.0 * a0u * b0uv + a0v * b0uu)
            + d2t.duv() * (2.0 * (a0u * b1uv + a1u * b0uv) + a0v * b1uu + a1v * b0uu)
            + d2t.dvv() * (2.0 * a1u * b1uv + a1v * b1uu);
        let mut suvv = suvv
            + d2t.duu() * (a0u * b0vv + 2.0 * a0v * b0uv)
            + d2t.duv() * (2.0 * (a0v * b1uv + a1v * b0uv) + a0u * b1vv + a1u * b0vv)
            + d2t.dvv() * (2.0 * a1v * b1uv + a1u * b1vv);
        let mut svvv = svvv
            + d2t.duu() * (3.0 * a0v * b0vv)
            + d2t.duv() * (3.0 * (a0v * b1vv + a1v * b0vv))
            + d2t.dvv() * (3.0 * a1v * b1vv);
        let _ = (&mut suuu, &mut suuv, &mut suvv, &mut svvv);

        let duuu = nsp * (normale.dot(suuu - d3s.duuu())) * invcos;
        let duuv = nsp * (normale.dot(suuv - d3s.duuv())) * invcos;
        let duvv = nsp * (normale.dot(suvv - d3s.duvv())) * invcos;
        let dvvv = nsp * (normale.dot(svvv - d3s.dvvv())) * invcos;

        constraint.my_ppc.push(PinpointConstraint::new(point2d, duuu, 3, 0));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duuv, 2, 1));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, duvv, 1, 2));
        constraint.my_ppc.push(PinpointConstraint::new(point2d, dvvv, 0, 3));
        constraint.nb_pp_constraints = 9;
        constraint
    }

    /// OCCT nb_PPC() (lxx L110-113).
    pub fn nb_ppc(&self) -> usize {
        self.nb_pp_constraints
    }

    /// OCCT GetPPC(Index) (lxx L115-118) — 0-based Index per OCCT.
    pub fn get_ppc(&self, index: usize) -> PinpointConstraint {
        self.my_ppc[index]
    }

    /// OCCT D1SurfInit() (lxx L120-123).
    pub fn d1_surf_init(&self) -> PlateD1 {
        self.my_d1_surf_init
    }
}

// ---------------------------------------------------------------------------
// Plate_FreeGtoCConstraint (Plate_FreeGtoCConstraint.cxx, whole file)
// ---------------------------------------------------------------------------

/// OCCT Plate_FreeGtoCConstraint — G1/G2/G3 constraint on the plate using a
/// weaker formulation than GtoCConstraint.
#[derive(Debug, Clone, PartialEq)]
pub struct FreeGtoCConstraint {
    pnt2d: DVec2,
    nb_pp_constraints: usize,
    nb_ls_constraints: usize,
    my_ppc: Vec<PinpointConstraint>, // OCCT myPPC[5]
    my_lsc: Vec<LinearScalarConstraint>, // OCCT myLSC[4]
}

impl FreeGtoCConstraint {
    /// OCCT G1 ctor (Plate_FreeGtoCConstraint.cxx L37-101).
    pub fn new(
        point2d: DVec2,
        d1s: &PlateD1,
        d1t: &PlateD1,
        incremental_load: f64,
        orientation: i32,
    ) -> Self {
        let mut constraint = FreeGtoCConstraint {
            pnt2d: point2d,
            nb_pp_constraints: 0,
            nb_ls_constraints: 0,
            my_ppc: Vec::new(),
            my_lsc: Vec::new(),
        };

        let mut normale = d1t.du().cross(d1t.dv());
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        if incremental_load != 1.0 {
            let mut n0 = d1s.du().cross(d1s.dv());
            if n0.length() < NORMIN {
                return constraint;
            }
            n0 = n0 / n0.length();
            let mut n1 = normale;
            if orientation != 0 {
                n1 *= orientation as f64;
            }
            let mut c = n0.dot(n1);
            if orientation == 0 && c < 0.0 {
                c *= -1.0;
                n1 *= -1.0;
            }

            let s = n0.cross(n1).length();
            if s < 1.0e-2 && c < 0.0 {
                return constraint;
            }
            let angle = f64::atan2(c, s);

            let mut d = n0.cross(n1);
            d = d / d.length();
            // OCCT: rota.SetRotation(gp_Ax1(gp_Pnt(0,0,0), dir),
            //                        angle * (IncrementalLoad - 1.));
            let rota_ang = angle * (incremental_load - 1.0);
            normale = rotate_about_origin(d, rota_ang, normale);
        }

        let du = d1s.du() * -1.0;
        let dv = d1s.dv() * -1.0;

        constraint.my_lsc.push(LinearScalarConstraint::new(
            PinpointConstraint::new(point2d, du, 1, 0),
            normale,
        ));
        constraint.my_lsc.push(LinearScalarConstraint::new(
            PinpointConstraint::new(point2d, dv, 0, 1),
            normale,
        ));
        constraint.nb_ls_constraints = 2;
        constraint
    }

    /// OCCT G1+G2 ctor (Plate_FreeGtoCConstraint.cxx L105-256).
    pub fn new_g2(
        point2d: DVec2,
        d1s: &PlateD1,
        d1t0: &PlateD1,
        d2s: &PlateD2,
        d2t0: &PlateD2,
        incremental_load: f64,
        orientation: i32,
    ) -> Self {
        let mut constraint = FreeGtoCConstraint {
            pnt2d: point2d,
            nb_pp_constraints: 0,
            nb_ls_constraints: 0,
            my_ppc: Vec::new(),
            my_lsc: Vec::new(),
        };
        let mut d1t = *d1t0;
        let mut d2t = *d2t0;

        let mut normale = d1t.du().cross(d1t.dv());
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        // G1 Constraints
        let mut normale_s = d1s.du().cross(d1s.dv());
        if normale_s.length() < NORMIN {
            if incremental_load != 1.0 {
                return constraint;
            }
            let du = d1s.du() * -1.0;
            let dv = d1s.dv() * -1.0;

            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, du, 1, 0),
                    normale,
                ));
            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, dv, 0, 1),
                    normale,
                ));
            constraint.nb_ls_constraints = 2;
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        if incremental_load != 1.0 {
            let n0 = normale_s;
            let mut n1 = normale;
            if orientation != 0 {
                n1 *= orientation as f64;
            }
            let mut c = n0.dot(n1);
            if orientation == 0 && c < 0.0 {
                c *= -1.0;
                n1 *= -1.0;
            }

            let s = n0.cross(n1).length();
            if s < 1.0e-2 && c < 0.0 {
                return constraint;
            }
            let angle = f64::atan2(c, s);

            let mut d = n0.cross(n1);
            d = d / d.length();
            let rota_ang = angle * (incremental_load - 1.0);
            normale = rotate_about_origin(d, rota_ang, normale);
            d1t = PlateD1::new(
                rotate_about_origin(d, rota_ang, d1t.du()),
                rotate_about_origin(d, rota_ang, d1t.dv()),
            );
            d2t = PlateD2::new(
                rotate_about_origin(d, rota_ang, d2t.duu()),
                rotate_about_origin(d, rota_ang, d2t.duv()),
                rotate_about_origin(d, rota_ang, d2t.dvv()),
            );
        }

        let cos_normales = normale.dot(normale_s);
        if cos_normales.abs() < COSMIN {
            let du = d1s.du() * -1.0;
            let dv = d1s.dv() * -1.0;

            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, du, 1, 0),
                    normale,
                ));
            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, dv, 0, 1),
                    normale,
                ));
            constraint.nb_ls_constraints = 2;
            return constraint;
        }

        let invcos = 1.0 / cos_normales;

        let du = normale_s * (-(normale.dot(d1s.du()))) * invcos;
        let dv = normale_s * (-(normale.dot(d1s.dv()))) * invcos;

        constraint
            .my_ppc
            .push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint
            .my_ppc
            .push(PinpointConstraint::new(point2d, dv, 0, 1));
        constraint.nb_pp_constraints = 2;

        // G2 Constraints
        let su = d1s.du() + du;
        let sv = d1s.dv() + dv;

        let mut mat = MatD::new(2, 2);
        mat.set(1, 1, su.dot(d1t.du()));
        mat.set(1, 2, su.dot(d1t.dv()));
        mat.set(2, 1, sv.dot(d1t.du()));
        mat.set(2, 2, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return constraint;
        }

        let mut vec = VecD::new(2);
        vec.set(1, su.dot(su));
        vec.set(2, su.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let a = sol.get(1);
        let b = sol.get(2);

        vec.set(1, sv.dot(su));
        vec.set(2, sv.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let c = sol.get(1);
        let d = sol.get(2);

        let suu = d2t.duu() * (a * a) + d2t.duv() * (2.0 * a * b) + d2t.dvv() * (b * b);
        let suv = d2t.duu() * (a * c) + d2t.duv() * (a * d + b * c) + d2t.dvv() * (b * d);
        let svv = d2t.duu() * (c * c) + d2t.duv() * (2.0 * c * d) + d2t.dvv() * (d * d);

        let mut duu = suu - d2s.duu();
        let mut duv = suv - d2s.duv();
        let mut dvv = svv - d2s.dvv();
        duu *= incremental_load;
        duv *= incremental_load;
        dvv *= incremental_load;

        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duu, 2, 0),
                normale,
            ));
        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duv, 1, 1),
                normale,
            ));
        constraint.nb_ls_constraints = 2;
        constraint
    }

    /// OCCT G1+G2+G3 ctor (Plate_FreeGtoCConstraint.cxx L255-463).
    pub fn new_g3(
        point2d: DVec2,
        d1s: &PlateD1,
        d1t0: &PlateD1,
        d2s: &PlateD2,
        d2t0: &PlateD2,
        d3s: &PlateD3,
        d3t0: &PlateD3,
        incremental_load: f64,
        orientation: i32,
    ) -> Self {
        let mut constraint = FreeGtoCConstraint {
            pnt2d: point2d,
            nb_pp_constraints: 0,
            nb_ls_constraints: 0,
            my_ppc: Vec::new(),
            my_lsc: Vec::new(),
        };
        let mut d1t = *d1t0;
        let mut d2t = *d2t0;
        let mut d3t = *d3t0;

        let mut normale = d1t.du().cross(d1t.dv());
        if normale.length() < NORMIN {
            return constraint;
        }
        normale = normale / normale.length();

        // G1 Constraints
        let mut normale_s = d1s.du().cross(d1s.dv());
        if normale_s.length() < NORMIN {
            if incremental_load != 1.0 {
                return constraint;
            }
            let du = d1s.du() * -1.0;
            let dv = d1s.dv() * -1.0;

            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, du, 1, 0),
                    normale,
                ));
            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, dv, 0, 1),
                    normale,
                ));
            constraint.nb_ls_constraints = 2;
            return constraint;
        }
        normale_s = normale_s / normale_s.length();

        if incremental_load != 1.0 {
            let n0 = normale_s;
            let mut n1 = normale;
            if orientation != 0 {
                n1 *= orientation as f64;
            }
            let mut c = n0.dot(n1);
            if orientation == 0 && c < 0.0 {
                c *= -1.0;
                n1 *= -1.0;
            }
            let s = n0.cross(n1).length();
            if s < 1.0e-2 && c < 0.0 {
                return constraint;
            }
            let angle = f64::atan2(c, s);

            let mut d = n0.cross(n1);
            d = d / d.length();
            let rota_ang = angle * (incremental_load - 1.0);
            normale = rotate_about_origin(d, rota_ang, normale);
            d1t = PlateD1::new(
                rotate_about_origin(d, rota_ang, d1t.du()),
                rotate_about_origin(d, rota_ang, d1t.dv()),
            );
            d2t = PlateD2::new(
                rotate_about_origin(d, rota_ang, d2t.duu()),
                rotate_about_origin(d, rota_ang, d2t.duv()),
                rotate_about_origin(d, rota_ang, d2t.dvv()),
            );
            d3t = PlateD3::new(
                rotate_about_origin(d, rota_ang, d3t.duuu()),
                rotate_about_origin(d, rota_ang, d3t.duuv()),
                rotate_about_origin(d, rota_ang, d3t.duvv()),
                rotate_about_origin(d, rota_ang, d3t.dvvv()),
            );
        }

        let cos_normales = normale.dot(normale_s);
        if cos_normales.abs() < COSMIN {
            let du = d1s.du() * -1.0;
            let dv = d1s.dv() * -1.0;

            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, du, 1, 0),
                    normale,
                ));
            constraint
                .my_lsc
                .push(LinearScalarConstraint::new(
                    PinpointConstraint::new(point2d, dv, 0, 1),
                    normale,
                ));
            constraint.nb_ls_constraints = 2;
            return constraint;
        }

        let invcos = 1.0 / cos_normales;

        let du = normale_s * (-(normale.dot(d1s.du()))) * invcos;
        let dv = normale_s * (-(normale.dot(d1s.dv()))) * invcos;

        constraint
            .my_ppc
            .push(PinpointConstraint::new(point2d, du, 1, 0));
        constraint
            .my_ppc
            .push(PinpointConstraint::new(point2d, dv, 0, 1));
        constraint.nb_pp_constraints = 2;

        // G2 Constraints
        let su = d1s.du() + du;
        let sv = d1s.dv() + dv;

        let mut mat = MatD::new(2, 2);
        mat.set(1, 1, su.dot(d1t.du()));
        mat.set(1, 2, su.dot(d1t.dv()));
        mat.set(2, 1, sv.dot(d1t.du()));
        mat.set(2, 2, sv.dot(d1t.dv()));
        let gauss = MathGauss::new(&mat);
        if !gauss.is_done() {
            return constraint;
        }

        let mut vec = VecD::new(2);
        vec.set(1, su.dot(su));
        vec.set(2, su.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let a = sol.get(1);
        let b = sol.get(2);

        vec.set(1, sv.dot(su));
        vec.set(2, sv.dot(sv));

        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let c = sol.get(1);
        let d = sol.get(2);

        let suu = d2t.duu() * (a * a) + d2t.duv() * (2.0 * a * b) + d2t.dvv() * (b * b);
        let suv = d2t.duu() * (a * c) + d2t.duv() * (a * d + b * c) + d2t.dvv() * (b * d);
        let svv = d2t.duu() * (c * c) + d2t.duv() * (2.0 * c * d) + d2t.dvv() * (d * d);

        let mut duu = suu - d2s.duu();
        let mut duv = suv - d2s.duv();
        let mut dvv = svv - d2s.dvv();
        duu *= incremental_load;
        duv *= incremental_load;
        dvv *= incremental_load;

        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duu, 2, 0),
                normale,
            ));
        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duv, 1, 1),
                normale,
            ));

        // G3 Constraints

        vec.set(1, (d2s.duu() + duu - suu).dot(su));
        vec.set(2, (d2s.duu() + duu - suu).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0uu = sol.get(1);
        let b1uu = sol.get(2);

        vec.set(1, (d2s.duv() + duv - suv).dot(su));
        vec.set(2, (d2s.duv() + duv - suv).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0uv = sol.get(1);
        let b1uv = sol.get(2);

        vec.set(1, (d2s.dvv() + dvv - svv).dot(su));
        vec.set(2, (d2s.dvv() + dvv - svv).dot(sv));
        let mut sol = vec.clone();
        gauss.solve(&mut sol);
        let b0vv = sol.get(1);
        let b1vv = sol.get(2);

        let mut suuu = d3t.duuu() * (a * a * a)
            + d3t.duuv() * (3.0 * a * a * b)
            + d3t.duvv() * (3.0 * a * b * b)
            + d3t.dvvv() * (b * b * b);
        let mut suuv = d3t.duuu() * (a * a * c)
            + d3t.duuv() * (a * a * d + 2.0 * a * b * c)
            + d3t.duvv() * (b * b * c + 2.0 * a * b * d)
            + d3t.dvvv() * (b * b * d);
        let mut suvv = d3t.duuu() * (a * c * c)
            + d3t.duuv() * (b * c * c + 2.0 * a * c * d)
            + d3t.duvv() * (a * d * d + 2.0 * b * c * d)
            + d3t.dvvv() * (b * d * d);
        let mut svvv = d3t.duuu() * (c * c * c)
            + d3t.duuv() * (3.0 * c * c * d)
            + d3t.duvv() * (3.0 * c * d * d)
            + d3t.dvvv() * (d * d * d);

        // OCCT references: A0u = a, A1u = b, A0v = c, A1v = d.
        let (a0u, a1u, a0v, a1v) = (a, b, c, d);
        suuu += d2t.duu() * (3.0 * a0u * b0uu)
            + d2t.duv() * (3.0 * (a0u * b1uu + a1u * b0uu))
            + d2t.dvv() * (3.0 * a1u * b1uu);
        suuv += d2t.duu() * (2.0 * a0u * b0uv + a0v * b0uu)
            + d2t.duv() * (2.0 * (a0u * b1uv + a1u * b0uv) + a0v * b1uu + a1v * b0uu)
            + d2t.dvv() * (2.0 * a1u * b1uv + a1v * b1uu);
        suvv += d2t.duu() * (a0u * b0vv + 2.0 * a0v * b0uv)
            + d2t.duv() * (2.0 * (a0v * b1uv + a1v * b0uv) + a0u * b1vv + a1u * b0vv)
            + d2t.dvv() * (2.0 * a1v * b1uv + a1u * b1vv);
        svvv += d2t.duu() * (3.0 * a0v * b0vv)
            + d2t.duv() * (3.0 * (a0v * b1vv + a1v * b0vv))
            + d2t.dvv() * (3.0 * a1v * b1vv);

        let duuu = suuu - d3s.duuu();
        let duuv = suuv - d3s.duuv();
        let duvv = suvv - d3s.duvv();
        let dvvv = svvv - d3s.dvvv();
        let duuu = duuu * incremental_load;
        let duuv = duuv * incremental_load;
        let duvv = duvv * incremental_load;
        let dvvv = dvvv * incremental_load;

        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duuu, 3, 0),
                normale,
            ));
        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duuv, 2, 1),
                normale,
            ));
        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, duvv, 1, 2),
                normale,
            ));
        constraint
            .my_lsc
            .push(LinearScalarConstraint::new(
                PinpointConstraint::new(point2d, dvvv, 0, 3),
                normale,
            ));
        constraint.nb_ls_constraints = 4;
        constraint
    }

    /// OCCT nb_PPC() (lxx L17-20).
    pub fn nb_ppc(&self) -> usize {
        self.nb_pp_constraints
    }

    /// OCCT GetPPC(Index) (lxx L22-25) — 0-based Index per OCCT.
    pub fn get_ppc(&self, index: usize) -> PinpointConstraint {
        self.my_ppc[index]
    }

    /// OCCT nb_LSC() (lxx L27-30).
    pub fn nb_lsc(&self) -> usize {
        self.nb_ls_constraints
    }

    /// OCCT LSC(Index) (lxx L32-35) — 0-based Index per OCCT.
    pub fn lsc(&self, index: usize) -> LinearScalarConstraint {
        self.my_lsc[index].clone()
    }
}

/// OCCT gp_Trsf::SetRotation (gp_Trsf.cxx L90-99) with the axis through the
/// origin, applied to a vector — gp_Mat::SetRotation (gp_Mat.cxx L122-158).
fn rotate_about_origin(axis_dir: DVec3, ang: f64, v: DVec3) -> DVec3 {
    let a_v = axis_dir / axis_dir.length();
    let a = a_v.x;
    let b = a_v.y;
    let c = a_v.z;

    let a_cos = ang.cos();
    let a_sin = ang.sin();
    let a_om_cos = 1.0 - a_cos;

    let a2 = a * a;
    let b2 = b * b;
    let c2 = c * c;
    let ab = a * b;
    let ac = a * c;
    let bc = b * c;

    // gp_Mat rows (L148-157).
    let m00 = 1.0 + a_om_cos * (-(b2 + c2));
    let m01 = a_om_cos * ab - a_sin * c;
    let m02 = a_om_cos * ac + a_sin * b;
    let m10 = a_om_cos * ab + a_sin * c;
    let m11 = 1.0 + a_om_cos * (-(a2 + c2));
    let m12 = a_om_cos * bc - a_sin * a;
    let m20 = a_om_cos * ac - a_sin * b;
    let m21 = a_om_cos * bc + a_sin * a;
    let m22 = 1.0 + a_om_cos * (-(a2 + b2));

    DVec3::new(
        m00 * v.x + m01 * v.y + m02 * v.z,
        m10 * v.x + m11 * v.y + m12 * v.z,
        m20 * v.x + m21 * v.y + m22 * v.z,
    )
}
