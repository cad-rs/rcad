//! OCCT Convert_GridPolynomialToPoles (TKMath/Convert,
//! Convert_GridPolynomialToPoles.cxx L28-422 + .hxx) — the conversion of a
//! grid of 2-var polynomial coefficient patches into a BSpline surface's
//! pole grid, consumed by AdvApp2Var (Patch::Poles, Patch::ConvertBS /
//! AdvApp2Var_ApproxAFunc2Var::ConvertBS).
//!
//! Architecture differences (Rust <-> C++):
//! - The handle-based delegating ctors (the HArray1/HArray2 overloads,
//!   cxx L28-43 and L45-71) collapse: rcad callers pass slices.
//! - `NCollection_Array1/Array2` become `Vec` with explicit 1-based
//!   indexing helpers (the `ati`-style offset convention of bspl_lib.rs).
//! - `BSplSLib::Interpolate` (the non-homogeneous surface form,
//!   BSplSLib.cxx L3955-4015) is translated here as `bspl_slib_interpolate`
//!   over the kernel `bspl_lib::interpolate` (BSplCLib.cxx L3353-3395)
//!   two-pass structure.

use glam::DVec3;

use super::bspl_lib::{build_schoenberg_points, interpolate as bspl_clib_interpolate, knot_sequence};
use super::plib::eval_poly2_var;

/// 1-based read helper (NCollection Array1::Value semantics).
fn at1(v: &[f64], i: i32) -> f64 {
    v[(i - 1) as usize]
}
fn at1i(v: &[i32], i: i32) -> i32 {
    v[(i - 1) as usize]
}

/// OCCT Convert_GridPolynomialToPoles (the .hxx member set).
#[derive(Debug, Clone, Default)]
pub struct ConvertGridPolynomialToPoles {
    /// OCCT: int myUDegree / myVDegree.
    my_u_degree: i32,
    my_v_degree: i32,
    /// OCCT: bool myDone.
    my_done: bool,
    /// OCCT: NCollection_Array1<double> myUKnots / myVKnots.
    my_u_knots: Vec<f64>,
    my_v_knots: Vec<f64>,
    /// OCCT: NCollection_Array1<double> myUFlatKnots / myVFlatKnots.
    my_u_flat_knots: Vec<f64>,
    my_v_flat_knots: Vec<f64>,
    /// OCCT: NCollection_Array1<int> myUMults / myVMults.
    my_u_mults: Vec<i32>,
    my_v_mults: Vec<i32>,
    /// OCCT: NCollection_Array2<gp_Pnt> myPoles (row = U pole index,
    /// column = V pole index).
    my_poles: Vec<Vec<DVec3>>,
}

impl ConvertGridPolynomialToPoles {
    /// OCCT ctor #3 (cxx L73-120): the single-patch form (continuity 0).
    #[allow(clippy::too_many_arguments)]
    pub fn from_single_patch(
        the_max_u_degree: i32,
        the_max_v_degree: i32,
        the_num_coeff: &[i32],
        the_coefficients: &[f64],
        the_polynomial_u_intervals: &[f64],
        the_polynomial_v_intervals: &[f64],
    ) -> Self {
        let mut conv = ConvertGridPolynomialToPoles {
            my_u_degree: 0,
            my_v_degree: 0,
            my_done: false,
            my_u_knots: Vec::new(),
            my_v_knots: Vec::new(),
            my_u_flat_knots: Vec::new(),
            my_v_flat_knots: Vec::new(),
            my_u_mults: Vec::new(),
            my_v_mults: Vec::new(),
            my_poles: Vec::new(),
        };

        // L84-92: the domain checks.
        if the_num_coeff.len() != 2 {
            panic!("Convert : Wrong Coefficients");
        }
        if the_coefficients.len() != (3 * (the_max_u_degree + 1) * (the_max_v_degree + 1)) as usize
        {
            panic!("Convert : Wrong Coefficients");
        }

        // L94-104: the actual degrees.
        conv.my_u_degree = at1i(the_num_coeff, 1) - 1;
        conv.my_v_degree = at1i(the_num_coeff, 2) - 1;
        if conv.my_u_degree > the_max_u_degree {
            panic!("Convert : Incoherence between NumCoeffPerSurface and MaxUDegree");
        }
        if conv.my_v_degree > the_max_v_degree {
            panic!("Convert : Incoherence between NumCoeffPerSurface and MaxVDegree");
        }

        // L106-108: NCollection_Array2<int> aNumCoeff2D(1, 1, 1, 2).
        let a_num_coeff_2d = vec![at1i(the_num_coeff, 1), at1i(the_num_coeff, 2)];

        // L110-119: Perform(0, 0, ...) with the polynomial intervals as the
        // true intervals.
        conv.perform(
            0,
            0,
            the_max_u_degree,
            the_max_v_degree,
            &a_num_coeff_2d,
            the_coefficients,
            the_polynomial_u_intervals,
            the_polynomial_v_intervals,
            the_polynomial_u_intervals,
            the_polynomial_v_intervals,
        );
        conv
    }

    /// OCCT ctor #4 (cxx L122-188): the 12-argument grid form.
    #[allow(clippy::too_many_arguments)]
    pub fn from_grid(
        the_nb_u_surfaces: i32,
        the_nb_v_surfaces: i32,
        the_u_continuity: i32,
        the_v_continuity: i32,
        the_max_u_degree: i32,
        the_max_v_degree: i32,
        the_num_coeff_per_surface: &[i32],
        the_coefficients: &[f64],
        the_polynomial_u_intervals: &[f64],
        the_polynomial_v_intervals: &[f64],
        the_true_u_intervals: &[f64],
        the_true_v_intervals: &[f64],
    ) -> Self {
        let mut conv = ConvertGridPolynomialToPoles {
            my_u_degree: 0,
            my_v_degree: 0,
            my_done: false,
            my_u_knots: Vec::new(),
            my_v_knots: Vec::new(),
            my_u_flat_knots: Vec::new(),
            my_v_flat_knots: Vec::new(),
            my_u_mults: Vec::new(),
            my_v_mults: Vec::new(),
            my_poles: Vec::new(),
        };

        // L139-140.
        let a_real_u_degree = the_max_u_degree.max(2 * the_u_continuity + 1);
        let a_real_v_degree = the_max_v_degree.max(2 * the_v_continuity + 1);

        // L142-147: the NumCoeffPerSurface bounds check (row count
        // NbUSurfaces * NbVSurfaces, 2 columns).
        if the_num_coeff_per_surface.len() != (the_nb_u_surfaces * the_nb_v_surfaces * 2) as usize
        {
            panic!("Convert : Wrong NumCoeffPerSurface");
        }

        // L149-154.
        if the_coefficients.len()
            != (3 * the_nb_u_surfaces * the_nb_v_surfaces
                * (a_real_u_degree + 1)
                * (a_real_v_degree + 1)) as usize
        {
            panic!("Convert : Wrong Coefficients");
        }

        // L156-167: compute the actual degrees from the coefficient counts.
        for ii in 1..=(the_nb_u_surfaces * the_nb_v_surfaces) {
            // Value(ii, 1).
            let nc_u = the_num_coeff_per_surface[((ii - 1) * 2) as usize];
            // Value(ii, 2).
            let nc_v = the_num_coeff_per_surface[((ii - 1) * 2 + 1) as usize];
            if nc_u > conv.my_u_degree + 1 {
                conv.my_u_degree = nc_u - 1;
            }
            if nc_v > conv.my_v_degree + 1 {
                conv.my_v_degree = nc_v - 1;
            }
        }

        // L169-176.
        if conv.my_u_degree > a_real_u_degree {
            panic!("Convert : Incoherence between NumCoeffPerSurface and MaxUDegree");
        }
        if conv.my_v_degree > a_real_v_degree {
            panic!("Convert : Incoherence between NumCoeffPerSurface and MaxVDegree");
        }

        // L178-187.
        conv.perform(
            the_u_continuity,
            the_v_continuity,
            a_real_u_degree,
            a_real_v_degree,
            the_num_coeff_per_surface,
            the_coefficients,
            the_polynomial_u_intervals,
            the_polynomial_v_intervals,
            the_true_u_intervals,
            the_true_v_intervals,
        );
        conv
    }

    /// OCCT Convert_GridPolynomialToPoles::Perform (cxx L190-293).
    #[allow(clippy::too_many_arguments)]
    fn perform(
        &mut self,
        the_u_continuity: i32,
        the_v_continuity: i32,
        the_max_u_degree: i32,
        the_max_v_degree: i32,
        the_num_coeff_per_surface: &[i32],
        the_coefficients: &[f64],
        the_polynomial_u_intervals: &[f64],
        the_polynomial_v_intervals: &[f64],
        the_true_u_intervals: &[f64],
        the_true_v_intervals: &[f64],
    ) {
        // (1) Build the monodimensional knot/multiplicity/parameter tables
        // (L202-209).
        self.my_u_knots = the_true_u_intervals.to_vec();
        self.my_v_knots = the_true_v_intervals.to_vec();

        let a_u_parameters =
            Self::build_array(self.my_u_degree, &self.my_u_knots, the_u_continuity, &mut self.my_u_flat_knots, &mut self.my_u_mults);
        let a_v_parameters =
            Self::build_array(self.my_v_degree, &self.my_v_knots, the_v_continuity, &mut self.my_v_flat_knots, &mut self.my_v_mults);

        // (2) Digitalisation (L211-279).
        let a_siz_patch = 3 * (the_max_u_degree + 1) * (the_max_v_degree + 1);
        let u_len = a_u_parameters.len();
        let v_len = a_v_parameters.len();
        self.my_poles = vec![vec![DVec3::ZERO; v_len]; u_len];

        let mut a_patch = vec![0.0f64; ((self.my_u_degree + 1) * 3 * (self.my_v_degree + 1)) as usize];
        let mut a_point = [0.0f64; 3];

        let mut a_patch_indice: i32 = 0;
        let u_knot_count = self.my_u_knots.len() as i32;
        let v_knot_count = self.my_v_knots.len() as i32;
        for ii in 1..=(u_len as i32) {
            let u_par = at1(&a_u_parameters, ii);
            let mut a_u_index: i32 = 1;
            while u_par > at1(the_true_u_intervals, a_u_index + 1) && a_u_index < u_knot_count - 1 {
                a_u_index += 1;
            }
            let a_n_value = (u_par - at1(the_true_u_intervals, a_u_index))
                / (at1(the_true_u_intervals, a_u_index + 1) - at1(the_true_u_intervals, a_u_index));
            let a_u_value = (1.0 - a_n_value) * at1(the_polynomial_u_intervals, 1)
                + a_n_value * at1(the_polynomial_u_intervals, 2);

            for jj in 1..=(v_len as i32) {
                let v_par = at1(&a_v_parameters, jj);
                let mut a_v_index: i32 = 1;
                while v_par > at1(the_true_v_intervals, a_v_index + 1)
                    && a_v_index < v_knot_count - 1
                {
                    a_v_index += 1;
                }
                let a_n_value = (v_par - at1(the_true_v_intervals, a_v_index))
                    / (at1(the_true_v_intervals, a_v_index + 1)
                        - at1(the_true_v_intervals, a_v_index));
                let a_v_value = (1.0 - a_n_value) * at1(the_polynomial_v_intervals, 1)
                    + a_n_value * at1(the_polynomial_v_intervals, 2);

                // (2.1) Extract the active patch coefficients (L248-264).
                if a_patch_indice != a_u_index + (u_knot_count - 1) * (a_v_index - 1) {
                    a_patch_indice = a_u_index + (u_knot_count - 1) * (a_v_index - 1);
                    let mut ll = 1usize;
                    // Value(aPatchIndice, 1) / Value(aPatchIndice, 2).
                    let nc_u = the_num_coeff_per_surface[((a_patch_indice - 1) * 2) as usize];
                    let nc_v =
                        the_num_coeff_per_surface[((a_patch_indice - 1) * 2 + 1) as usize];
                    for k1 in 1..=nc_u {
                        let mut pos = (a_siz_patch * (a_patch_indice - 1)
                            + 3 * (the_max_v_degree + 1) * (k1 - 1)
                            + 1) as usize;
                        for _k2 in 1..=nc_v {
                            a_patch[ll - 1] = the_coefficients[pos - 1];
                            a_patch[ll] = the_coefficients[pos];
                            a_patch[ll + 1] = the_coefficients[pos + 1];
                            ll += 3;
                            pos += 3;
                        }
                    }
                }

                // (2.2) Evaluate at (aUValue, aVValue) (L266-275):
                // PLib::EvalPoly2Var with the active patch's degrees.
                eval_poly2_var(
                    a_u_value,
                    a_v_value,
                    0,
                    0,
                    (the_num_coeff_per_surface[((a_patch_indice - 1) * 2) as usize] - 1) as usize,
                    (the_num_coeff_per_surface[((a_patch_indice - 1) * 2 + 1) as usize] - 1)
                        as usize,
                    3,
                    &a_patch,
                    &mut a_point,
                );

                self.my_poles[(ii - 1) as usize][(jj - 1) as usize] =
                    DVec3::new(a_point[0], a_point[1], a_point[2]);
            }
        }

        // (3) Interpolation (L281-292): BSplSLib::Interpolate.
        let mut inversion_problem = 0i32;
        bspl_slib_interpolate(
            self.my_u_degree,
            self.my_v_degree,
            &self.my_u_flat_knots,
            &self.my_v_flat_knots,
            &a_u_parameters,
            &a_v_parameters,
            &mut self.my_poles,
            &mut inversion_problem,
        );
        self.my_done = inversion_problem == 0;
    }

    /// OCCT Convert_GridPolynomialToPoles::BuildArray (cxx L295-327) — the
    /// multiplicities, flat knot sequence and Schoenberg interpolation
    /// parameters for one direction; returns the parameters.
    fn build_array(
        degree: i32,
        knots: &[f64],
        continuity: i32,
        flat_knots: &mut Vec<f64>,
        mults: &mut Vec<i32>,
    ) -> Vec<f64> {
        let num_curves = knots.len() as i32 - 1;

        // Calcul des Multiplicites (L304-313).
        let multiplicities = degree - continuity;
        *mults = vec![0i32; knots.len()];
        for ii in 2..(knots.len() as i32) {
            mults[(ii - 1) as usize] = multiplicities;
        }
        mults[0] = degree + 1;
        mults[(num_curves + 1 - 1) as usize] = degree + 1;

        // Calcul des Noeuds Plats (L315-319): BSplCLib::KnotSequence(Knots,
        // Mults, Degree, false, FlatKnots).
        let num_flat_knots = (multiplicities * (num_curves - 1) + 2 * degree + 2) as usize;
        let mut fk = vec![0.0f64; num_flat_knots];
        knot_sequence(knots, mults, degree as usize, false, &mut fk);
        *flat_knots = fk;

        // Calcul du nombre de Poles (L321-322).
        let num_poles = num_flat_knots as i32 - degree - 1;

        // Calcul des parametres d'interpolation (L324-326):
        // BSplCLib::BuildSchoenbergPoints.
        let mut parameters = vec![0.0f64; num_poles as usize];
        build_schoenberg_points(degree as usize, flat_knots, &mut parameters);
        parameters
    }

    /// OCCT IsDone() (cxx L419-422).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT UDegree() (cxx L355-359).
    pub fn u_degree(&self) -> i32 {
        self.my_u_degree
    }

    /// OCCT VDegree() (cxx L363-367).
    pub fn v_degree(&self) -> i32 {
        self.my_v_degree
    }

    /// OCCT NbUPoles() (cxx L331-335).
    pub fn nb_u_poles(&self) -> usize {
        self.my_poles.len()
    }

    /// OCCT NbVPoles() (cxx L339-343).
    pub fn nb_v_poles(&self) -> usize {
        self.my_poles.first().map(|r| r.len()).unwrap_or(0)
    }

    /// OCCT Poles() (cxx L347-351).
    pub fn poles(&self) -> &Vec<Vec<DVec3>> {
        &self.my_poles
    }

    /// OCCT NbUKnots() (cxx L371-375).
    pub fn nb_u_knots(&self) -> usize {
        self.my_u_knots.len()
    }

    /// OCCT NbVKnots() (cxx L379-383).
    pub fn nb_v_knots(&self) -> usize {
        self.my_v_knots.len()
    }

    /// OCCT UKnots() (cxx L387-391).
    pub fn u_knots(&self) -> &[f64] {
        &self.my_u_knots
    }

    /// OCCT VKnots() (cxx L395-399).
    pub fn v_knots(&self) -> &[f64] {
        &self.my_v_knots
    }

    /// OCCT UMultiplicities() (cxx L403-407).
    pub fn u_multiplicities(&self) -> &[i32] {
        &self.my_u_mults
    }

    /// OCCT VMultiplicities() (cxx L411-415).
    pub fn v_multiplicities(&self) -> &[i32] {
        &self.my_v_mults
    }
}

/// OCCT BSplSLib::Interpolate (the non-homogeneous surface form,
/// BSplSLib.cxx L3955-4015): two passes of BSplCLib::Interpolate — first
/// the iso-u extraction interpolated along V, then the iso-v extraction
/// interpolated along U.
#[allow(clippy::too_many_arguments)]
fn bspl_slib_interpolate(
    u_degree: i32,
    v_degree: i32,
    u_flat_knots: &[f64],
    v_flat_knots: &[f64],
    u_parameters: &[f64],
    v_parameters: &[f64],
    poles: &mut [Vec<DVec3>],
    inversion_problem: &mut i32,
) {
    let u_length = u_parameters.len() as i32;
    let v_length = v_parameters.len() as i32;

    // extraction of iso u (L3969-3982): Points(1, VLength, 1, 3*ULength);
    // each V row carries the U pole row flattened (x, y, z).
    let dimension = 3 * u_length;
    let mut points = vec![0.0f64; (v_length * dimension) as usize];
    for ii in 1..=(v_length as i32) {
        for (jj, ll) in (1i32..).zip((0i32..).step_by(3)).take(u_length as usize) {
            let p = poles[(jj - 1) as usize][(ii - 1) as usize];
            let base = ((ii - 1) * dimension + ll) as usize;
            points[base] = p.x;
            points[base + 1] = p.y;
            points[base + 2] = p.z;
        }
    }

    // interpolation of iso u (L3984-3990): BSplCLib::Interpolate along the
    // V direction with contact order 0.
    let contact_order = vec![0i32; v_length as usize];
    let error_code = bspl_clib_interpolate(
        v_degree as usize,
        v_flat_knots,
        v_parameters,
        &contact_order,
        dimension as usize,
        &mut points,
    );
    if error_code != 0 {
        *inversion_problem = error_code;
        return;
    }

    // extraction of iso v (L3993-4005): IsoPoles(1, ULength, 1,
    // 3*VLength); the transpose of the interpolated Points.
    let dimension = 3 * v_length;
    let mut iso_poles = vec![0.0f64; (u_length * dimension) as usize];
    for ii in 1..=(u_length as i32) {
        for (jj, ll) in (1i32..).zip((0i32..).step_by(3)).take(v_length as usize) {
            let base = ((ii - 1) * dimension + ll) as usize;
            let src = ((jj - 1) * (3 * u_length) + (ii - 1) * 3) as usize;
            iso_poles[base] = points[src];
            iso_poles[base + 1] = points[src + 1];
            iso_poles[base + 2] = points[src + 2];
        }
    }

    // interpolation of iso v (L4007-4013): BSplCLib::Interpolate along the
    // U direction.
    let contact_order = vec![0i32; u_length as usize];
    let error_code = bspl_clib_interpolate(
        u_degree as usize,
        u_flat_knots,
        u_parameters,
        &contact_order,
        dimension as usize,
        &mut iso_poles,
    );
    if error_code != 0 {
        *inversion_problem = error_code;
        return;
    }

    // return results (L4016-4027).
    for (ii, row) in (1i32..).zip(poles.iter_mut()).take(u_length as usize) {
        for (jj, p) in (1i32..).zip(row.iter_mut()).take(v_length as usize) {
            let base = ((ii - 1) * dimension + (jj - 1) * 3) as usize;
            *p = DVec3::new(iso_poles[base], iso_poles[base + 1], iso_poles[base + 2]);
        }
    }
}

// The banded-solver ladder (BuildBSpMatrix -> FactorBandedMatrix ->
// SolveBandedSystem) runs inside bspl_lib::interpolate per the OCCT
// BSplCLib::Interpolate body (BSplCLib.cxx L3353-3395).
