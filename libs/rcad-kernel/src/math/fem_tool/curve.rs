//! OCCT FEmTool_Curve (ModelingData/TKGeomBase/FEmTool).
//!
//! 1:1 translation of `FEmTool_Curve.hxx` (L17-100) and `FEmTool_Curve.cxx`
//! (L1-580).  Also hosts `PLib::EvalLength` (PLib.cxx L2244-2301), the PLib
//! free function consumed by `FEmTool_Curve::Length`, because the `math`
//! module registration is closed and the legacy `math/plib.rs` is not part
//! of this batch's edit scope.
//!
//! Container mapping: the OCCT 1-based flat storages (`myKnots`,
//! `myDegree`, `myCoeff`, `myPoly`, `myDeri`, `myDsecn`, `HasPoly`,
//! `HasDeri`, `HasSecn`, `myLength`) become `Vec`s indexed with explicit
//! `i - 1`; `NCollection_Array2<double>` element coefficients use the
//! bound-aware `math_matrix::Matrix`; output arrays (`Pnt`, `Vec`) are
//! 0-based slices (Lower() == 0 convention).

use super::super::gauss_points::{gauss_points, gauss_weights};
use super::super::math_matrix::Matrix;
use super::super::plib::hermit_jacobi::HermitJacobi;
use super::super::plib::no_derivative_eval_polynomial_flat;

/// OCCT PLib::EvalLength (PLib.cxx L2244-2301) - length of a polynomial
/// curve segment computed by Gauss integration.  `polynomial_coeff` is the
/// flat coefficient array addressed like the OCCT `double&` pointer.
pub fn eval_length(
    degree: i32,
    dimension: i32,
    polynomial_coeff: &[f64],
    u1: f64,
    u2: f64,
    length: &mut f64,
) {
    let nb_gauss_points = 4 * ((degree / 4) + 1).min(10);

    let mut gauss_points_v = super::super::math_matrix::Vector::new(1, nb_gauss_points);
    gauss_points(nb_gauss_points as usize, &mut gauss_points_v.data);

    let mut gauss_weights_v = super::super::math_matrix::Vector::new(1, nb_gauss_points);
    gauss_weights(nb_gauss_points as usize, &mut gauss_weights_v.data);

    let c1 = (u2 + u1) / 2.0;
    let c2 = (u2 - u1) / 2.0;

    //-----------------------------------------------------------
    //****** Integration - Boucle sur les intervalles de GAUSS **
    //-----------------------------------------------------------

    let mut sum = 0.0;

    for j in 1..=nb_gauss_points / 2 {
        // Integration en tenant compte de la symetrie
        let tran = c2 * gauss_points_v.get(j);
        let x1 = c1 + tran;
        let x2 = c1 - tran;

        //****** Derivation sur la dimension de l'espace **

        let degdim = (degree * dimension) as usize;
        let mut der1 = 0.0;
        let mut der2 = 0.0;
        for idim in 0..dimension {
            let idimu = idim as usize;
            let mut d1 = degree as f64 * polynomial_coeff[idimu + degdim];
            let mut d2 = d1;
            for i in (1..degree).rev() {
                let dd = i as f64 * polynomial_coeff[idimu + (i * dimension) as usize];
                d1 = d1 * x1 + dd;
                d2 = d2 * x2 + dd;
            }
            der1 += d1 * d1;
            der2 += d2 * d2;
        }

        //****** Integration **

        sum += gauss_weights_v.get(j) * c2 * (der1.sqrt() + der2.sqrt());

        //****** Fin de boucle dur les intervalles de GAUSS **
    }
    *length = sum;
}

/// OCCT FEmTool_Curve - curve defined by polynomial elements
/// (FEmTool_Curve.hxx L32-98).
#[derive(Debug, Clone)]
pub struct Curve {
    my_nb_elements: i32,
    my_dimension: i32,
    my_base: HermitJacobi,
    /// OCCT myKnots HArray1<double>(1, myNbElements + 1).
    my_knots: Vec<f64>,
    /// OCCT myDegree NCollection_Array1<int>(1, myNbElements).
    my_degree: Vec<i32>,
    /// OCCT myCoeff NCollection_Array1<double>(1, dim*nb*(wd+1)).
    my_coeff: Vec<f64>,
    /// OCCT myPoly NCollection_Array1<double>(1, dim*nb*(wd+1)).
    my_poly: Vec<f64>,
    /// OCCT myDeri NCollection_Array1<double>(1, dim*nb*wd).
    my_deri: Vec<f64>,
    /// OCCT myDsecn NCollection_Array1<double>(1, dim*nb*(wd-1)).
    my_dsecn: Vec<f64>,
    /// OCCT HasPoly NCollection_Array1<int>(1, myNbElements).
    has_poly: Vec<i32>,
    /// OCCT HasDeri NCollection_Array1<int>(1, myNbElements).
    has_deri: Vec<i32>,
    /// OCCT HasSecn NCollection_Array1<int>(1, myNbElements).
    has_secn: Vec<i32>,
    /// OCCT myLength NCollection_Array1<double>(1, myNbElements).
    my_length: Vec<f64>,
    uf: f64,
    ul: f64,
    denom: f64,
    usum: f64,
    my_index: i32,
    my_ptr: i32,
}

impl Curve {
    /// OCCT FEmTool_Curve::FEmTool_Curve (FEmTool_Curve.cxx L30-54).  OCCT
    /// copies `TheBase` and ignores the (unnamed) Tolerance parameter.
    pub fn new(dimension: i32, nb_elements: i32, the_base: &HermitJacobi, _tolerance: f64) -> Self {
        let my_nb_elements = nb_elements;
        let my_dimension = dimension;
        let wd = the_base.work_degree();
        let deg_base = (wd + 1) as usize;
        let my_base = the_base.clone();

        let my_knots = vec![0.0f64; (my_nb_elements + 1) as usize];
        // myDegree.Init(myBase.WorkDegree());
        let my_degree = vec![wd; my_nb_elements as usize];
        let has_poly = vec![0i32; my_nb_elements as usize];
        let has_deri = vec![0i32; my_nb_elements as usize];
        let has_secn = vec![0i32; my_nb_elements as usize];
        let my_length = vec![-1.0f64; my_nb_elements as usize];

        Curve {
            my_nb_elements,
            my_dimension,
            my_base,
            my_knots,
            my_degree,
            my_coeff: vec![0.0; my_dimension as usize * my_nb_elements as usize * deg_base],
            my_poly: vec![0.0; my_dimension as usize * my_nb_elements as usize * deg_base],
            my_deri: vec![0.0; my_dimension as usize * my_nb_elements as usize * wd as usize],
            my_dsecn: vec![0.0; my_dimension as usize * my_nb_elements as usize * (wd - 1) as usize],
            has_poly,
            has_deri,
            has_secn,
            my_length,
            uf: 0.0,
            ul: 0.0,
            denom: 0.0,
            usum: 0.0,
            my_index: 0,
            my_ptr: 0,
        }
    }

    /// OCCT FEmTool_Curve::Knots (cxx L56-59) - OCCT returns
    /// `NCollection_Array1<double>&` so callers fill the knots; the slice is
    /// 0-based (OCCT Knots()(i) is `knots()[i - 1]`).
    pub fn knots(&mut self) -> &mut [f64] {
        &mut self.my_knots
    }

    /// OCCT FEmTool_Curve::NbElements (cxx L470-473).
    pub fn nb_elements(&self) -> i32 {
        self.my_nb_elements
    }

    /// OCCT FEmTool_Curve::Dimension (cxx L475-478).
    pub fn dimension(&self) -> i32 {
        self.my_dimension
    }

    /// OCCT FEmTool_Curve::Base (cxx L480-483).
    pub fn base(&self) -> &HermitJacobi {
        &self.my_base
    }

    /// OCCT FEmTool_Curve::Degree (cxx L485-488).
    pub fn degree(&self, index_of_element: i32) -> i32 {
        self.my_degree[(index_of_element - 1) as usize]
    }

    /// OCCT FEmTool_Curve::SetElement (cxx L63-102).
    pub fn set_element(&mut self, index_of_element: i32, coeffs: &Matrix) {
        assert!(
            index_of_element <= self.my_nb_elements && index_of_element >= 1,
            "Standard_OutOfRange: FEmTool_Curve::SetElement"
        );
        let deg_base = self.my_base.work_degree();
        let deg = self.my_degree[(index_of_element - 1) as usize];
        let i_base = ((index_of_element - 1) * (deg_base + 1) * self.my_dimension) as usize;
        let mut i1 = i_base as i32 - self.my_dimension;
        // OCCT: i2 = Coeffs.LowerRow() - 1, j1 = Coeffs.LowerCol() - 1; the
        // Matrix bound-aware access makes the explicit offsets unnecessary.
        let mut i2 = coeffs.lower_row() - 1;
        let j1 = coeffs.lower_col() - 1;
        for _i in 1..=deg + 1 {
            i1 += self.my_dimension;
            i2 += 1;
            for j in 1..=self.my_dimension {
                self.my_coeff[(i1 + j) as usize - 1] = coeffs.get(i2, j1 + j);
            }
        }

        let stenor = (self.my_knots[index_of_element as usize]
            - self.my_knots[(index_of_element - 1) as usize])
            / 2.0;
        let mut i1 = i_base as i32;
        let mut i2 = i1 + (self.my_base.niv_constr() + 1) * self.my_dimension;
        for i in 1..=self.my_base.niv_constr() {
            i1 += self.my_dimension;
            i2 += self.my_dimension;
            let mfact = stenor.powi(i);
            for j in 1..=self.my_dimension {
                self.my_coeff[(i1 + j) as usize - 1] *= mfact;
                self.my_coeff[(i2 + j) as usize - 1] *= mfact;
            }
        }

        self.has_poly[(index_of_element - 1) as usize] = 0;
        self.has_deri[(index_of_element - 1) as usize] = 0;
        self.has_secn[(index_of_element - 1) as usize] = 0;
        self.my_length[(index_of_element - 1) as usize] = -1.0;
    }

    /// OCCT FEmTool_Curve::GetElement (cxx L106-141).
    pub fn get_element(&mut self, index_of_element: i32, coeffs: &mut Matrix) {
        assert!(
            index_of_element <= self.my_nb_elements && index_of_element >= 1,
            "Standard_OutOfRange: FEmTool_Curve::GetElement"
        );
        let deg_base = self.my_base.work_degree();
        let deg = self.my_degree[(index_of_element - 1) as usize];
        let i_base = ((index_of_element - 1) * (deg_base + 1) * self.my_dimension) as usize;
        let mut i1 = i_base as i32 - self.my_dimension;
        let mut i2 = coeffs.lower_row() - 1;
        let j1 = coeffs.lower_col() - 1;
        for _i in 1..=deg + 1 {
            i1 += self.my_dimension;
            i2 += 1;
            for j in 1..=self.my_dimension {
                coeffs.set(i2, j1 + j, self.my_coeff[(i1 + j) as usize - 1]);
            }
        }

        let stenor = 2.0
            / (self.my_knots[index_of_element as usize]
                - self.my_knots[(index_of_element - 1) as usize]);

        let i2 = coeffs.lower_row();
        let i3 = i2 + self.my_base.niv_constr() + 1;

        for i in 1..=self.my_base.niv_constr() {
            let mfact = stenor.powi(i);
            for j in (j1 + 1)..=self.my_dimension {
                let v = coeffs.get(i2 + i, j) * mfact;
                coeffs.set(i2 + i, j, v);
                let v = coeffs.get(i3 + i, j) * mfact;
                coeffs.set(i3 + i, j, v);
            }
        }
    }

    /// OCCT FEmTool_Curve::GetPolynom (cxx L145-162) - returns coefficients
    /// of all elements in canonical base (0-based slice convention).
    pub fn get_polynom(&mut self, coeffs: &mut [f64]) {
        for index_of_element in 1..=self.my_nb_elements {
            if self.has_poly[(index_of_element - 1) as usize] == 0 {
                self.update(index_of_element, 0);
            }
        }

        // OCCT: Coeffs(di + i) = myPoly(i), di = Coeffs.Lower() - myPoly.Lower();
        // both storages normalized to 0-based, so the indices coincide.
        let n = self.my_poly.len();
        coeffs[..n].copy_from_slice(&self.my_poly);
    }

    /// Search of the span (common prologue of D0 / D1 / D2,
    /// FEmTool_Curve.cxx L171-198 / L223-250 / L281-308).
    fn find_span(&mut self, u: f64) {
        if self.my_index == 0
            || u < self.uf
            || u > self.ul
            || self.my_knots[(self.my_index - 1) as usize] != self.uf
            || self.my_knots[self.my_index as usize] != self.ul
        {
            // Search the span
            if u <= self.my_knots[1] {
                self.my_index = 1;
            } else {
                self.my_index = 2;
                while self.my_index <= self.my_nb_elements {
                    if u >= self.my_knots[(self.my_index - 1) as usize]
                        && u <= self.my_knots[self.my_index as usize]
                    {
                        break;
                    }
                    self.my_index += 1;
                }
                if self.my_index > self.my_nb_elements {
                    self.my_index = self.my_nb_elements;
                }
            }
            self.uf = self.my_knots[(self.my_index - 1) as usize];
            self.ul = self.my_knots[self.my_index as usize];
            self.denom = 1.0 / (self.ul - self.uf);
            self.usum = self.uf + self.ul;
            self.my_ptr =
                (self.my_index - 1) * (self.my_base.work_degree() + 1) * self.my_dimension + 1;
        }
    }

    /// OCCT FEmTool_Curve::D0 (cxx L166-214).
    pub fn d0(&mut self, u: f64, pnt: &mut [f64]) {
        self.find_span(u);

        let deg = self.my_degree[(self.my_index - 1) as usize];
        if self.has_poly[(self.my_index - 1) as usize] == 0 {
            self.update(self.my_index, 0);
        }

        // Parameter normalization: S [-1, 1]
        let s = (2.0 * u - self.usum) * self.denom;
        no_derivative_eval_polynomial_flat(
            s,
            deg,
            self.my_dimension,
            deg * self.my_dimension,
            &self.my_poly[(self.my_ptr - 1) as usize..],
            pnt,
        );
    }

    /// OCCT FEmTool_Curve::D1 (cxx L218-272).
    pub fn d1(&mut self, u: f64, vec: &mut [f64]) {
        self.find_span(u);

        let deg = self.my_degree[(self.my_index - 1) as usize];
        if self.has_deri[(self.my_index - 1) as usize] == 0 {
            self.update(self.my_index, 1);
        }

        // Parameter normalization: S [-1, 1]
        let s = (2.0 * u - self.usum) * self.denom;
        let offset = (self.my_index - 1) * self.my_base.work_degree() * self.my_dimension;
        no_derivative_eval_polynomial_flat(
            s,
            deg - 1,
            self.my_dimension,
            (deg - 1) * self.my_dimension,
            &self.my_deri[offset as usize..],
            vec,
        );

        let s = 2.0 * self.denom;
        for v in vec.iter_mut() {
            *v *= s;
        }
    }

    /// OCCT FEmTool_Curve::D2 (cxx L276-331).
    pub fn d2(&mut self, u: f64, vec: &mut [f64]) {
        self.find_span(u);

        let deg = self.my_degree[(self.my_index - 1) as usize];
        if self.has_secn[(self.my_index - 1) as usize] == 0 {
            self.update(self.my_index, 2);
        }

        // Parameter normalization: S [-1, 1]
        let s = (2.0 * u - self.usum) * self.denom;
        let offset = (self.my_index - 1) * (self.my_base.work_degree() - 1) * self.my_dimension;
        no_derivative_eval_polynomial_flat(
            s,
            deg - 2,
            self.my_dimension,
            (deg - 2) * self.my_dimension,
            &self.my_dsecn[offset as usize..],
            vec,
        );

        let s = 4.0 * self.denom * self.denom;
        for v in vec.iter_mut() {
            *v *= s;
        }
    }

    /// OCCT FEmTool_Curve::Length (cxx L335-468).
    pub fn length(&mut self, first_u: f64, last_u: f64, length: &mut f64) {
        assert!(first_u <= last_u, "Standard_OutOfRange: FEmTool_Curve::Length");

        let low: i32;
        if self.my_knots[0] > first_u {
            low = 1;
        } else {
            let mut l = 1;
            while l <= self.my_nb_elements {
                if first_u >= self.my_knots[(l - 1) as usize]
                    && first_u <= self.my_knots[l as usize]
                {
                    break;
                }
                l += 1;
            }
            low = l;
        }
        let low = if low > self.my_nb_elements { self.my_nb_elements } else { low };

        let high: i32;
        if self.my_knots[0] > last_u {
            high = 1;
        } else {
            let mut h = low;
            while h <= self.my_nb_elements {
                if last_u >= self.my_knots[(h - 1) as usize] && last_u <= self.my_knots[h as usize]
                {
                    break;
                }
                h += 1;
            }
            high = h;
        }
        let mut high = high;
        if self.my_knots[self.my_nb_elements as usize] < last_u {
            high = self.my_nb_elements;
        }

        let deg_base = self.my_base.work_degree();
        *length = 0.0;

        let first_s = (2.0 * first_u - self.my_knots[(low - 1) as usize] - self.my_knots[low as usize])
            / (self.my_knots[low as usize] - self.my_knots[(low - 1) as usize]);
        let last_s = (2.0 * last_u - self.my_knots[(high - 1) as usize] - self.my_knots[high as usize])
            / (self.my_knots[high as usize] - self.my_knots[(high - 1) as usize]);

        if low == high {
            let ptr = ((low - 1) * (deg_base + 1) * self.my_dimension + 1) as usize;
            let deg = self.my_degree[(low - 1) as usize];

            if self.has_poly[(low - 1) as usize] == 0 {
                self.update(low, 0);
            }
            let mut li = 0.0;
            eval_length(
                deg,
                self.my_dimension,
                &self.my_poly[ptr - 1..],
                first_s,
                last_s,
                &mut li,
            );
            *length = li;
            return;
        }

        let mut deg = self.my_degree[(low - 1) as usize];
        let mut ptr = ((low - 1) * (deg_base + 1) * self.my_dimension + 1) as usize;

        if self.has_poly[(low - 1) as usize] == 0 {
            self.update(low, 0);
        }
        let mut li = 0.0;
        if first_s < -1.0 {
            eval_length(
                deg,
                self.my_dimension,
                &self.my_poly[ptr - 1..],
                first_s,
                -1.0,
                &mut li,
            );
            *length += li;
            if self.my_length[(low - 1) as usize] < 0.0 {
                eval_length(
                    deg,
                    self.my_dimension,
                    &self.my_poly[ptr - 1..],
                    -1.0,
                    1.0,
                    &mut li,
                );
                self.my_length[(low - 1) as usize] = li;
            }
            *length += self.my_length[(low - 1) as usize];
        } else {
            eval_length(
                deg,
                self.my_dimension,
                &self.my_poly[ptr - 1..],
                first_s,
                1.0,
                &mut li,
            );
            *length += li;
        }

        deg = self.my_degree[(high - 1) as usize];
        ptr = ((high - 1) * (deg_base + 1) * self.my_dimension + 1) as usize;

        if self.has_poly[(high - 1) as usize] == 0 {
            self.update(high, 0);
        }
        if last_s > 1.0 {
            eval_length(
                deg,
                self.my_dimension,
                &self.my_poly[ptr - 1..],
                1.0,
                last_s,
                &mut li,
            );
            *length += li;
            if self.my_length[(high - 1) as usize] < 0.0 {
                eval_length(
                    deg,
                    self.my_dimension,
                    &self.my_poly[ptr - 1..],
                    -1.0,
                    1.0,
                    &mut li,
                );
                self.my_length[(high - 1) as usize] = li;
            }
            *length += self.my_length[(high - 1) as usize];
        } else {
            eval_length(
                deg,
                self.my_dimension,
                &self.my_poly[ptr - 1..],
                -1.0,
                last_s,
                &mut li,
            );
            *length += li;
        }

        for i in (low + 1)..high {
            if self.my_length[(i - 1) as usize] < 0.0 {
                let ptr = ((i - 1) * (deg_base + 1) * self.my_dimension + 1) as usize;
                let deg = self.my_degree[(i - 1) as usize];
                if self.has_poly[(i - 1) as usize] == 0 {
                    self.update(i, 0);
                }
                let mut li = 0.0;
                eval_length(
                    deg,
                    self.my_dimension,
                    &self.my_poly[ptr - 1..],
                    -1.0,
                    1.0,
                    &mut li,
                );
                self.my_length[(i - 1) as usize] = li;
            }
            *length += self.my_length[(i - 1) as usize];
        }
    }

    /// OCCT FEmTool_Curve::SetDegree (cxx L492-504).
    pub fn set_degree(&mut self, index_of_element: i32, degree: i32) {
        if degree <= self.my_base.work_degree() {
            self.my_degree[(index_of_element - 1) as usize] = degree;
            self.has_poly[(index_of_element - 1) as usize] = 0;
            self.has_deri[(index_of_element - 1) as usize] = 0;
            self.has_secn[(index_of_element - 1) as usize] = 0;
            self.my_length[(index_of_element - 1) as usize] = -1.0;
        } else {
            panic!("Standard_OutOfRange: FEmTool_Curve::SetDegree");
        }
    }

    /// OCCT FEmTool_Curve::ReduceDegree (cxx L508-527).
    pub fn reduce_degree(
        &mut self,
        index_of_element: i32,
        tol: f64,
        new_degree: &mut i32,
        max_error: &mut f64,
    ) {
        let deg = self.my_degree[(index_of_element - 1) as usize];

        let ptr = ((index_of_element - 1) * (self.my_base.work_degree() + 1) * self.my_dimension + 1)
            as usize;

        self.my_base.reduce_degree(
            self.my_dimension,
            deg,
            tol,
            &self.my_coeff[ptr - 1..],
            new_degree,
            max_error,
        );

        *new_degree = (*new_degree).max(2 * self.my_base.niv_constr() + 1);

        if *new_degree < deg {
            self.my_degree[(index_of_element - 1) as usize] = *new_degree;
            self.has_poly[(index_of_element - 1) as usize] = 0;
            self.has_deri[(index_of_element - 1) as usize] = 0;
            self.has_secn[(index_of_element - 1) as usize] = 0;
            self.my_length[(index_of_element - 1) as usize] = -1.0;
        }
    }

    /// OCCT FEmTool_Curve::Update (cxx L531-580).
    fn update(&mut self, index: i32, order: i32) {
        let deg_base = self.my_base.work_degree();
        let deg = self.my_degree[(index - 1) as usize];

        if self.has_poly[(index - 1) as usize] == 0 {
            let ptr = ((index - 1) * (deg_base + 1) * self.my_dimension + 1) as usize;
            let n = (self.my_dimension * (deg + 1)) as usize;

            let coeff = &mut self.my_poly[ptr - 1..ptr - 1 + n];
            let base_coeff = &self.my_coeff[ptr - 1..ptr - 1 + n];

            self.my_base.to_coefficients(self.my_dimension, deg, base_coeff, &mut *coeff);
            self.has_poly[(index - 1) as usize] = 1;
        }

        if order >= 1 {
            let mut i1 = (index - 1) * deg_base * self.my_dimension - self.my_dimension;
            let mut i2 = (index - 1) * (deg_base + 1) * self.my_dimension;
            if self.has_deri[(index - 1) as usize] == 0 {
                for i in 1..=deg {
                    i1 += self.my_dimension;
                    i2 += self.my_dimension;
                    for j in 1..=self.my_dimension {
                        self.my_deri[(i1 + j) as usize - 1] =
                            i as f64 * self.my_poly[(i2 + j) as usize - 1];
                    }
                }
                self.has_deri[(index - 1) as usize] = 1;
            }
            if order >= 2 && self.has_secn[(index - 1) as usize] == 0 {
                i1 = (index - 1) * (deg_base - 1) * self.my_dimension - self.my_dimension;
                i2 = (index - 1) * deg_base * self.my_dimension;
                for i in 1..deg {
                    i1 += self.my_dimension;
                    i2 += self.my_dimension;
                    for j in 1..=self.my_dimension {
                        self.my_dsecn[(i1 + j) as usize - 1] =
                            i as f64 * self.my_deri[(i2 + j) as usize - 1];
                    }
                }
                self.has_secn[(index - 1) as usize] = 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::GeomAbsShape;
    use super::super::super::math_matrix::Matrix;
    use super::super::super::plib::hermit_jacobi::HermitJacobi;
    use super::*;

    /// FEmTool_Curve D0/D1/D2 on a single HermitJacobi element, hand-derived.
    /// Base: PLib_HermitJacobi(2, C0) (NivConstr = 0, DegreeH = 1, Jacobi
    /// degree 0, W(t) = 1 - t^2, J0 = sqrt(15/16)); element on knots
    /// [0, 2]; coefficients (H00, H01, J0) = (1, 2, 3).
    /// With s = (2u - 2)/2 = u - 1:
    ///   P(s) = H00 + 2*H01 + 3*TN0*(1 - s^2),  TN0 = sqrt(15/16)
    ///   dP/ds = -1/2 + 1 + 3*TN0*(-2s) = 1/2 - 6*TN0*s
    ///   d2P/ds2 = -6*TN0
    /// At u = 0.5 (s = -0.5), with ds/du = 1:
    ///   D0 = 3/4 + 2/4 + 3*TN0*3/4 = 5/4 + (9/4)*TN0
    ///   D1 = 1/2 + 3*TN0
    ///   D2 = -6*TN0
    #[test]
    fn fem_tool_curve_d0_d1_d2_single_element_hand() {
        let base = HermitJacobi::new(2, GeomAbsShape::C0);
        let mut curve = Curve::new(1, 1, &base, 0.0);
        curve.knots()[0] = 0.0;
        curve.knots()[1] = 2.0;

        // Element coefficients: rows 1..3 = (H00, H01, J0) coefficients.
        let mut coeffs = Matrix::new(1, 3, 1, 1);
        coeffs.set(1, 1, 1.0);
        coeffs.set(2, 1, 2.0);
        coeffs.set(3, 1, 3.0);
        curve.set_element(1, &coeffs);

        let tn0 = (15.0f64 / 16.0).sqrt();

        let mut pnt = [0.0f64; 1];
        curve.d0(0.5, &mut pnt);
        let expected0 = 1.25 + 2.25 * tn0;
        assert!((pnt[0] - expected0).abs() < 1e-13, "D0 = {}", pnt[0]);

        let mut der = [0.0f64; 1];
        curve.d1(0.5, &mut der);
        let expected1 = 0.5 + 3.0 * tn0;
        assert!((der[0] - expected1).abs() < 1e-13, "D1 = {}", der[0]);

        let mut sec = [0.0f64; 1];
        curve.d2(0.5, &mut sec);
        let expected2 = -6.0 * tn0;
        assert!((sec[0] - expected2).abs() < 1e-13, "D2 = {}", sec[0]);
    }

    /// FEmTool_Curve::Length hand-check on the same element: with
    /// coefficients (H00, H01, J0) = (1, 2, 0) the canonical polynomial is
    /// P(s) = 3/2 + s/2, so |P'| = 1/2 over the whole span and the length of
    /// [0, 2] is exactly 1 (and the sub-segment [0.5, 2] is exactly 3/4).
    #[test]
    fn fem_tool_curve_length_hand() {
        let base = HermitJacobi::new(2, GeomAbsShape::C0);
        let mut curve = Curve::new(1, 1, &base, 0.0);
        curve.knots()[0] = 0.0;
        curve.knots()[1] = 2.0;

        let mut coeffs = Matrix::new(1, 3, 1, 1);
        coeffs.set(1, 1, 1.0);
        coeffs.set(2, 1, 2.0);
        coeffs.set(3, 1, 0.0);
        curve.set_element(1, &coeffs);

        let mut len = 0.0;
        curve.length(0.0, 2.0, &mut len);
        assert!((len - 1.0).abs() < 1e-13, "length = {}", len);

        curve.length(0.5, 2.0, &mut len);
        assert!((len - 0.75).abs() < 1e-13, "length = {}", len);
    }
}
