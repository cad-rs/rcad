//! OCCT Plate_Plate (TKGeomAlgo/Plate) — 1:1 port of Plate_Plate.cxx
//! (whole file, L34-2001) + Plate_Plate.lxx inline bodies.
//!
//! Architecture mappings: gp_XY -> DVec2, gp_XYZ -> DVec3,
//! NCollection_Sequence -> Vec, the void* member buffers (solution/points/
//! deru/derv) -> Option<Vec<...>>, math_Matrix/math_Vector -> MatD/VecD
//! (Plate indexes them 0-based; the 1-based MatD/VecD accessors take +1),
//! math_Gauss -> MathGauss, Message_ProgressRange/Scope -> dropped (no
//! progress reporting; UserBreak() is never true).
//!
//! The OCCT `mutable` scratch members (Uold/Vold/U2/R/L) mutated through a
//! const_cast in SolEm become `Cell<f64>` so that Evaluate/SolEm keep `&self`.

use core::cell::Cell;

use glam::{DVec2, DVec3};

use super::constraints::{LinearScalarConstraint, LinearXYZConstraint};
use super::pinpoint_constraint::PinpointConstraint;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{MatD, VecD};

/// OCCT math_Matrix x math_Vector product (row-major dot per row).
fn mat_vec_mul(mat: &MatD, v: &VecD) -> VecD {
    let n = mat.n_rows();
    let m = mat.n_cols();
    let mut out = VecD::new(n);
    for i in 1..=n {
        let mut sum = 0.0f64;
        for j in 1..=m {
            sum += mat.get(i, j) * v.get(j);
        }
        out.set(i, sum);
    }
    out
}

/// OCCT math_Vector addition (sol += sol1).
fn vec_add(a: &mut VecD, b: &VecD) {
    for i in 1..=a.len() {
        a.set(i, a.get(i) + b.get(i));
    }
}

/// OCCT math_Vector subtraction (sec_member1 = sec_member - mat * sol).
fn vec_sub(a: &VecD, b: &VecD) -> VecD {
    let mut out = VecD::new(a.len());
    for i in 1..=a.len() {
        out.set(i, a.get(i) - b.get(i));
    }
    out
}

/// OCCT gp_XYZ::SetCoord(icoor, v) on a DVec3.
fn xyz_set_coord(v: &mut DVec3, icoor: i32, value: f64) {
    match icoor {
        1 => v.x = value,
        2 => v.y = value,
        _ => v.z = value,
    }
}

/// OCCT gp_XYZ::Coord(icoor).
fn xyz_coord(v: DVec3, icoor: i32) -> f64 {
    match icoor {
        1 => v.x,
        2 => v.y,
        _ => v.z,
    }
}

/// OCCT Plate_Plate.
#[derive(Clone)]
pub struct Plate {
    order: i32,
    n_el: i32,
    n_dim: i32,
    solution: Option<Vec<DVec3>>,
    points: Option<Vec<DVec2>>,
    deru: Option<Vec<i32>>,
    derv: Option<Vec<i32>>,
    ok: bool,
    my_constraints: Vec<PinpointConstraint>,
    my_lxyz_constraints: Vec<LinearXYZConstraint>,
    my_lscalar_constraints: Vec<LinearScalarConstraint>,
    ddu: [f64; 10],
    ddv: [f64; 10],
    max_constraint_order: i32,
    polynomial_part_only: bool,
    // OCCT mutable scratch (mutated through const_cast in SolEm).
    u_old: Cell<f64>,
    v_old: Cell<f64>,
    u2: Cell<f64>,
    r: Cell<f64>,
    l: Cell<f64>,
}

impl Default for Plate {
    fn default() -> Self {
        Self::new()
    }
}

impl Plate {
    /// OCCT Plate_Plate() (Plate_Plate.cxx L34-53).
    pub fn new() -> Self {
        Plate {
            order: 0,
            n_el: 0,
            n_dim: 0,
            solution: None,
            points: None,
            deru: None,
            derv: None,
            ok: false,
            my_constraints: Vec::new(),
            my_lxyz_constraints: Vec::new(),
            my_lscalar_constraints: Vec::new(),
            ddu: [0.0; 10],
            ddv: [0.0; 10],
            max_constraint_order: 0,
            polynomial_part_only: false,
            u_old: Cell::new(1.0e20),
            v_old: Cell::new(1.0e20),
            u2: Cell::new(0.0),
            r: Cell::new(0.0),
            l: Cell::new(0.0),
        }
    }

    /// OCCT Load(PinpointConstraint) (Plate_Plate.cxx L195-205).
    pub fn load_pinpoint(&mut self, p_const: PinpointConstraint) {
        self.ok = false;
        self.n_el += 1;
        self.my_constraints.push(p_const);
        let ordre_const = p_const.idu() + p_const.idv();
        if self.max_constraint_order < ordre_const {
            self.max_constraint_order = ordre_const;
        }
    }

    /// OCCT Load(LinearXYZConstraint) (Plate_Plate.cxx L207-221).
    pub fn load_linear_xyz(&mut self, lxyz_const: LinearXYZConstraint) {
        self.ok = false;
        self.n_el += lxyz_const.coeff().row_len() as i32;

        self.my_lxyz_constraints.push(lxyz_const.clone());
        for j in 1..=lxyz_const.get_ppc().len() {
            let ordre_const = lxyz_const.get_ppc()[j - 1].idu()
                + lxyz_const.get_ppc()[j - 1].idv();
            if self.max_constraint_order < ordre_const {
                self.max_constraint_order = ordre_const;
            }
        }
    }

    /// OCCT Load(LinearScalarConstraint) (Plate_Plate.cxx L223-236).
    pub fn load_linear_scalar(&mut self, lscalar_const: LinearScalarConstraint) {
        self.ok = false;
        self.n_el += lscalar_const.coeff().row_len() as i32;
        self.my_lscalar_constraints.push(lscalar_const.clone());
        for j in 1..=lscalar_const.get_ppc().len() {
            let ordre_const = lscalar_const.get_ppc()[j - 1].idu()
                + lscalar_const.get_ppc()[j - 1].idv();
            if self.max_constraint_order < ordre_const {
                self.max_constraint_order = ordre_const;
            }
        }
    }

    /// OCCT Load(GtoCConstraint) (Plate_Plate.cxx L253-259).
    pub fn load_gto_c(&mut self, gto_c_const: &super::constraints::GtoCConstraint) {
        for i in 0..gto_c_const.nb_ppc() {
            self.load_pinpoint(gto_c_const.get_ppc(i));
        }
    }

    /// OCCT Load(FreeGtoCConstraint) (Plate_Plate.cxx L261-272).
    pub fn load_free_gto_c(&mut self, f_gto_c_const: &super::constraints::FreeGtoCConstraint) {
        for i in 0..f_gto_c_const.nb_ppc() {
            self.load_pinpoint(f_gto_c_const.get_ppc(i));
        }
        for i in 0..f_gto_c_const.nb_lsc() {
            self.load_linear_scalar(f_gto_c_const.lsc(i));
        }
    }

    /// OCCT Load(GlobalTranslationConstraint) (Plate_Plate.cxx L274-277).
    pub fn load_global_translation(
        &mut self,
        gt_const: &super::constraints::GlobalTranslationConstraint,
    ) {
        let lxyzc = gt_const.lxyzc().clone();
        self.load_linear_xyz(lxyzc);
    }

    /// OCCT Load(LineConstraint) (Plate_Plate.cxx L238-241).
    pub fn load_line(&mut self, l_const: &super::constraints::LineConstraint) {
        let lsc = l_const.lsc().clone();
        self.load_linear_scalar(lsc);
    }

    /// OCCT Load(PlaneConstraint) (Plate_Plate.cxx L243-246).
    pub fn load_plane(&mut self, p_const: &super::constraints::PlaneConstraint) {
        let lsc = p_const.lsc().clone();
        self.load_linear_scalar(lsc);
    }

    /// OCCT Load(SampledCurveConstraint) (Plate_Plate.cxx L248-251).
    pub fn load_sampled_curve(&mut self, sc_const: &super::constraints::SampledCurveConstraint) {
        let lxyzc = sc_const.lxyzc().clone();
        self.load_linear_xyz(lxyzc);
    }

    /// OCCT SolveTI (Plate_Plate.cxx L284-362).  Message_ProgressRange is
    /// dropped (no progress reporting in rcad).
    pub fn solve_ti(&mut self, ord: i32, anisotropie: f64) {
        let iteration_number = 0;
        self.ok = false;
        self.order = ord;
        if ord <= 1 {
            return;
        }
        if ord > 9 {
            return;
        }
        if self.n_el < 1 {
            return;
        }
        if anisotropie < 1.0e-6 {
            return;
        }
        if anisotropie > 1.0e6 {
            return;
        }

        // computation of the bounding box of the 2d PPconstraints
        let (mut xmin, mut xmax, mut ymin, mut ymax) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        self.uv_box(&mut xmin, &mut xmax, &mut ymin, &mut ymax);

        let mut du = 0.5 * (xmax - xmin);
        if anisotropie > 1.0 {
            du *= anisotropie;
        }
        if du < 1.0e-10 {
            return;
        }
        self.ddu[0] = 1.0;
        for i in 1..=9 {
            self.ddu[i] = self.ddu[i - 1] / du;
        }

        let mut dv = 0.5 * (ymax - ymin);
        if anisotropie < 1.0 {
            dv /= anisotropie;
        }
        if dv < 1.0e-10 {
            return;
        }
        self.ddv[0] = 1.0;
        for i in 1..=9 {
            self.ddv[i] = self.ddv[i - 1] / dv;
        }

        if self.my_lscalar_constraints.is_empty() {
            if self.my_lxyz_constraints.is_empty() {
                self.solve_ti1(iteration_number);
            } else {
                self.solve_ti2(iteration_number);
            }
        } else {
            self.solve_ti3(iteration_number);
        }
    }

    /// OCCT SolveTI1 (Plate_Plate.cxx L370-505).
    fn solve_ti1(&mut self, iteration_number: i32) {
        // computation of square matrix members

        self.n_dim = self.n_el + self.order * (self.order + 1) / 2;
        let n_dim = self.n_dim as usize;
        let mut mat = MatD::new(n_dim, n_dim);

        self.points = Some(vec![DVec2::ZERO; self.n_el as usize]);
        let n_el = self.n_el as usize;
        {
            let points = self.points.as_mut().expect("points");
            for i in 0..n_el {
                points[i] = self.my_constraints[i].pnt2d();
            }
        }

        self.deru = Some(vec![0i32; n_el]);
        {
            let deru = self.deru.as_mut().expect("deru");
            for i in 0..n_el {
                deru[i] = self.my_constraints[i].idu();
            }
        }

        self.derv = Some(vec![0i32; n_el]);
        {
            let derv = self.derv.as_mut().expect("derv");
            for i in 0..n_el {
                derv[i] = self.my_constraints[i].idv();
            }
        }

        {
            let deru = self.deru.as_ref().expect("deru");
            let derv = self.derv.as_ref().expect("derv");
            let points = self.points.as_ref().expect("points");
            for i in 0..n_el {
                for j in 0..i {
                    let signe: f64 = if (deru[j] + derv[j]) % 2 == 1 {
                        -1.0
                    } else {
                        1.0
                    };
                    let iu = deru[i] + deru[j];
                    let iv = derv[i] + derv[j];
                    mat.set(
                        (i + 1) as usize,
                        (j + 1) as usize,
                        signe * self.sol_em(points[i] - points[j], iu, iv),
                    );
                }
            }
        }

        let mut i = n_el;
        for iu in 0..self.order {
            let mut iv = 0i32;
            while iu + iv < self.order {
                {
                    let deru = self.deru.as_ref().expect("deru");
                    let derv = self.derv.as_ref().expect("derv");
                    let points = self.points.as_ref().expect("points");
                    for j in 0..n_el {
                        let idu = deru[j];
                        let idv = derv[j];
                        mat.set(
                            (i + 1) as usize,
                            (j + 1) as usize,
                            self.polm(points[j], iu, iv, idu, idv),
                        );
                    }
                }
                i += 1;
                iv += 1;
            }
        }

        for i in 0..n_dim {
            for j in (i + 1)..n_dim {
                let v = mat.get((j + 1) as usize, (i + 1) as usize);
                mat.set((i + 1) as usize, (j + 1) as usize, v);
            }
        }

        // initialisation of the Gauss algorithm
        let mut pivot_max = 1.0e-12;
        self.ok = true;

        // OCCT: math_Gauss algo_gauss(mat, pivot_max, aScope.Next(7));
        let mut algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);

        if !algo_gauss.is_done() {
            let nbm = self.order * (self.order + 1) / 2;
            for i in self.n_el..(self.n_el + nbm) {
                mat.set((i + 1) as usize, (i + 1) as usize, 1.0e-8);
            }
            pivot_max = 1.0e-18;

            let the_algo = MathGauss::with_min_pivot(&mat, pivot_max);
            algo_gauss = the_algo;
            self.ok = algo_gauss.is_done();
        }

        if self.ok {
            //   computation of the linear system solution for the X, Y and Z
            //   coordinates
            let mut sec_member = VecD::new(n_dim);
            let mut sol = VecD::new(n_dim);

            self.solution = Some(vec![DVec3::ZERO; n_dim]);

            for icoor in 1..=3 {
                for i in 0..n_el {
                    sec_member.set(
                        (i + 1) as usize,
                        xyz_coord(self.my_constraints[i].value(), icoor),
                    );
                }
                let mut sol_ = sec_member.clone();
                algo_gauss.solve(&mut sol_);
                sol = sol_;
                // alr iteration pour affiner la solution
                {
                    let mut sol1 = VecD::new(n_dim);
                    let mut sec_member1 = VecD::new(n_dim);
                    for _ in 1..=iteration_number {
                        sec_member1 = vec_sub(&sec_member, &mat_vec_mul(&mat, &sol));
                        sol1 = sec_member1.clone();
                        algo_gauss.solve(&mut sol1);
                        vec_add(&mut sol, &sol1);
                    }
                }
                // finalr

                for i in 0..n_dim {
                    let v = sol.get((i + 1) as usize);
                    let solution = self.solution.as_mut().expect("solution");
                    xyz_set_coord(&mut solution[i], icoor, v);
                }
            }
        }
    }

    /// OCCT SolveTI2 (Plate_Plate.cxx L513-663).
    fn solve_ti2(&mut self, iteration_number: i32) {
        // computation of square matrix members

        let n_cc1 = self.my_constraints.len();
        let mut n_cc2 = 0usize;
        for lxyz in &self.my_lxyz_constraints {
            n_cc2 += lxyz.coeff().col_len();
        }

        let n_dimat = n_cc1 + n_cc2 + (self.order * (self.order + 1) / 2) as usize;

        self.points = Some(vec![DVec2::ZERO; self.n_el as usize]);
        self.deru = Some(vec![0i32; self.n_el as usize]);
        self.derv = Some(vec![0i32; self.n_el as usize]);

        {
            let points = self.points.as_mut().expect("points");
            let deru = self.deru.as_mut().expect("deru");
            let derv = self.derv.as_mut().expect("derv");
            for i in 0..n_cc1 {
                points[i] = self.my_constraints[i].pnt2d();
                deru[i] = self.my_constraints[i].idu();
                derv[i] = self.my_constraints[i].idv();
            }
        }

        let mut k = n_cc1;
        for lxyz in &self.my_lxyz_constraints {
            for ppc in lxyz.get_ppc() {
                let points = self.points.as_mut().expect("points");
                let deru = self.deru.as_mut().expect("deru");
                let derv = self.derv.as_mut().expect("derv");
                points[k] = ppc.pnt2d();
                deru[k] = ppc.idu();
                derv[k] = ppc.idv();
                k += 1;
            }
        }

        let mut mat = MatD::new(n_dimat, n_dimat);

        self.fill_xyz_matrix(&mut mat, 0, 0, n_cc1, n_cc2);

        // initialisation of the Gauss algorithm
        let mut pivot_max = 1.0e-12;
        self.ok = true; // ************ JHH

        // OCCT: math_Gauss algo_gauss(mat, pivot_max, aScope.Next(7));
        let mut algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);

        if !algo_gauss.is_done() {
            for i in (n_cc1 + n_cc2)..n_dimat {
                mat.set(i + 1, i + 1, 1.0e-8);
            }
            pivot_max = 1.0e-18;

            let the_algo1 = MathGauss::with_min_pivot(&mat, pivot_max);
            algo_gauss = the_algo1;
            self.ok = algo_gauss.is_done();
        }

        if self.ok {
            //   computation of the linear system solution for the X, Y and Z
            //   coordinates
            let mut sec_member = VecD::new(n_dimat);
            let mut sol = VecD::new(n_dimat);

            self.n_dim = self.n_el + self.order * (self.order + 1) / 2;
            self.solution = Some(vec![DVec3::ZERO; self.n_dim as usize]);

            for icoor in 1..=3 {
                for i in 0..n_cc1 {
                    sec_member.set(
                        i + 1,
                        xyz_coord(self.my_constraints[i].value(), icoor),
                    );
                }

                let mut k = n_cc1;
                for lxyz in &self.my_lxyz_constraints {
                    for irow in 1..=lxyz.coeff().col_len() {
                        for icol in 1..=lxyz.coeff().row_len() {
                            let v = sec_member.get(k + 1)
                                + lxyz.coeff().get(irow, icol)
                                    * xyz_coord(lxyz.get_ppc()[icol - 1].value(), icoor);
                            sec_member.set(k + 1, v);
                        }
                        k += 1;
                    }
                }

                let mut sol_ = sec_member.clone();
                algo_gauss.solve(&mut sol_);
                sol = sol_;
                // alr iteration pour affiner la solution
                {
                    let mut sol1 = VecD::new(n_dimat);
                    let mut sec_member1 = VecD::new(n_dimat);
                    for _ in 1..=iteration_number {
                        sec_member1 = vec_sub(&sec_member, &mat_vec_mul(&mat, &sol));
                        sol1 = sec_member1.clone();
                        algo_gauss.solve(&mut sol1);
                        vec_add(&mut sol, &sol1);
                    }
                }
                // finalr

                for i in 0..n_cc1 {
                    let v = sol.get(i + 1);
                    let solution = self.solution.as_mut().expect("solution");
                    xyz_set_coord(&mut solution[i], icoor, v);
                }

                let mut k_solution = n_cc1;
                let mut ksol = n_cc1;

                for lxyz in &self.my_lxyz_constraints {
                    for icol in 1..=lxyz.coeff().row_len() {
                        let mut vsol = 0.0f64;
                        for irow in 1..=lxyz.coeff().col_len() {
                            vsol += lxyz.coeff().get(irow, icol)
                                * sol.get(ksol + irow);
                        }
                        let solution = self.solution.as_mut().expect("solution");
                        xyz_set_coord(&mut solution[k_solution], icoor, vsol);
                        k_solution += 1;
                    }
                    ksol += lxyz.coeff().col_len();
                }

                for i in 0..(self.order * (self.order + 1) / 2) as usize {
                    let v = sol.get(ksol + i + 1);
                    let solution = self.solution.as_mut().expect("solution");
                    xyz_set_coord(&mut solution[self.n_el as usize + i], icoor, v);
                }
            }
        }
    }

    /// OCCT SolveTI3 (Plate_Plate.cxx L670-1056).
    fn solve_ti3(&mut self, iteration_number: i32) {
        // computation of square matrix members

        let n_cc1 = self.my_constraints.len();

        let mut n_cc2 = 0usize;
        for lxyz in &self.my_lxyz_constraints {
            n_cc2 += lxyz.coeff().col_len();
        }

        let mut n_cc3 = 0usize;
        for lscalar in &self.my_lscalar_constraints {
            n_cc3 += lscalar.coeff().col_len();
        }

        let nbm = (self.order * (self.order + 1) / 2) as usize;
        let n_dimsousmat = n_cc1 + n_cc2 + nbm;
        let n_dimat = 3 * n_dimsousmat + n_cc3;

        self.points = Some(vec![DVec2::ZERO; self.n_el as usize]);
        self.deru = Some(vec![0i32; self.n_el as usize]);
        self.derv = Some(vec![0i32; self.n_el as usize]);

        {
            let points = self.points.as_mut().expect("points");
            let deru = self.deru.as_mut().expect("deru");
            let derv = self.derv.as_mut().expect("derv");
            for i in 0..n_cc1 {
                points[i] = self.my_constraints[i].pnt2d();
                deru[i] = self.my_constraints[i].idu();
                derv[i] = self.my_constraints[i].idv();
            }
        }

        let mut k = n_cc1;
        for lxyz in &self.my_lxyz_constraints {
            for ppc in lxyz.get_ppc() {
                let points = self.points.as_mut().expect("points");
                let deru = self.deru.as_mut().expect("deru");
                let derv = self.derv.as_mut().expect("derv");
                points[k] = ppc.pnt2d();
                deru[k] = ppc.idu();
                derv[k] = ppc.idv();
                k += 1;
            }
        }
        let n_ppc2 = k;
        for lscalar in &self.my_lscalar_constraints {
            for ppc in lscalar.get_ppc() {
                let points = self.points.as_mut().expect("points");
                let deru = self.deru.as_mut().expect("deru");
                let derv = self.derv.as_mut().expect("derv");
                points[k] = ppc.pnt2d();
                deru[k] = ppc.idu();
                derv[k] = ppc.idv();
                k += 1;
            }
        }

        let mut mat = MatD::new(n_dimat, n_dimat);

        self.fill_xyz_matrix(&mut mat, 0, 0, n_cc1, n_cc2);
        self.fill_xyz_matrix(&mut mat, n_dimsousmat as i32, n_dimsousmat as i32, n_cc1, n_cc2);
        self.fill_xyz_matrix(
            &mut mat,
            2 * n_dimsousmat as i32,
            2 * n_dimsousmat as i32,
            n_cc1,
            n_cc2,
        );

        let mut k = 3 * n_dimsousmat;
        let mut kppc = n_ppc2;
        for i in 1..=self.my_lscalar_constraints.len() {
            let lscalar = &self.my_lscalar_constraints[i - 1];
            for j in 0..n_cc1 {
                let mut vmat = VecD::new(lscalar.get_ppc().len());

                for ippc in 1..=lscalar.get_ppc().len() {
                    let deru = self.deru.as_ref().expect("deru");
                    let derv = self.derv.as_ref().expect("derv");
                    let points = self.points.as_ref().expect("points");
                    let mut signe = 1.0f64;
                    if (deru[j] + derv[j]) % 2 == 1 {
                        signe = -1.0;
                    }
                    let iu = deru[kppc + ippc - 1] + deru[j];
                    let iv = derv[kppc + ippc - 1] + derv[j];
                    vmat.set(ippc as usize, signe * self.sol_em(points[kppc + ippc - 1] - points[j], iu, iv));
                }

                for irow in 1..=lscalar.coeff().col_len() {
                    for icol in 1..=lscalar.coeff().row_len() {
                        let coeff = lscalar.coeff().get(irow, icol);
                        let v = mat.get(k + irow, j + 1) + coeff.x * vmat.get(icol as usize);
                        mat.set(k + irow, j + 1, v);
                        let v = mat.get(k + irow, n_dimsousmat + j + 1)
                            + coeff.y * vmat.get(icol as usize);
                        mat.set(k + irow, n_dimsousmat + j + 1, v);
                        let v = mat.get(k + irow, 2 * n_dimsousmat + j + 1)
                            + coeff.z * vmat.get(icol as usize);
                        mat.set(k + irow, 2 * n_dimsousmat + j + 1, v);
                    }
                }
            }

            let mut k2 = n_cc1;
            let mut kppc2 = n_cc1;
            for i2 in 1..=i {
                let lxyz2 = &self.my_lxyz_constraints[i2 - 1];

                let mut tmpmat = MatD::new(lscalar.get_ppc().len(), lxyz2.get_ppc().len());

                for ippc in 1..=lscalar.get_ppc().len() {
                    for ippc2 in 1..=lxyz2.get_ppc().len() {
                        let deru = self.deru.as_ref().expect("deru");
                        let derv = self.derv.as_ref().expect("derv");
                        let points = self.points.as_ref().expect("points");
                        let mut signe = 1.0f64;
                        if (deru[kppc2 + ippc2 - 1] + derv[kppc2 + ippc2 - 1]) % 2 == 1 {
                            signe = -1.0;
                        }
                        let iu = deru[kppc + ippc - 1] + deru[kppc2 + ippc2 - 1];
                        let iv = derv[kppc + ippc - 1] + derv[kppc2 + ippc2 - 1];
                        tmpmat.set(
                            ippc as usize,
                            ippc2 as usize,
                            signe
                                * self.sol_em(
                                    points[kppc + ippc - 1] - points[kppc2 + ippc2 - 1],
                                    iu,
                                    iv,
                                ),
                        );
                    }
                }

                for irow in 1..=lscalar.coeff().col_len() {
                    for irow2 in 1..=lxyz2.coeff().col_len() {
                        for icol in 1..=lscalar.coeff().row_len() {
                            for icol2 in 1..=lxyz2.coeff().row_len() {
                                let v = mat.get(k + irow, k2 + irow2)
                                    + lscalar.coeff().get(irow, icol).x
                                        * lxyz2.coeff().get(irow2, icol2)
                                        * tmpmat.get(icol as usize, icol2 as usize);
                                mat.set(k + irow, k2 + irow2, v);
                                let v = mat.get(k + irow, n_dimsousmat + k2 + irow2)
                                    + lscalar.coeff().get(irow, icol).y
                                        * lxyz2.coeff().get(irow2, icol2)
                                        * tmpmat.get(icol as usize, icol2 as usize);
                                mat.set(k + irow, n_dimsousmat + k2 + irow2, v);
                                let v = mat.get(k + irow, 2 * n_dimsousmat + k2 + irow2)
                                    + lscalar.coeff().get(irow, icol).z
                                        * lxyz2.coeff().get(irow2, icol2)
                                        * tmpmat.get(icol as usize, icol2 as usize);
                                mat.set(k + irow, 2 * n_dimsousmat + k2 + irow2, v);
                            }
                        }
                    }
                }

                k2 += lxyz2.coeff().col_len();
                kppc2 += lxyz2.coeff().row_len();
            }

            let mut j = n_cc1 + n_cc2;
            for iu in 0..self.order {
                let mut iv = 0i32;
                while iu + iv < self.order {
                    let lscalar = &self.my_lscalar_constraints[i - 1];
                    let mut vmat = VecD::new(lscalar.get_ppc().len());
                    {
                        let deru = self.deru.as_ref().expect("deru");
                        let derv = self.derv.as_ref().expect("derv");
                        let points = self.points.as_ref().expect("points");
                        for ippc in 1..=lscalar.get_ppc().len() {
                            let idu = deru[kppc + ippc - 1];
                            let idv = derv[kppc + ippc - 1];
                            vmat.set(
                                ippc as usize,
                                self.polm(points[kppc + ippc - 1], iu, iv, idu, idv),
                            );
                        }
                    }

                    for irow in 1..=lscalar.coeff().col_len() {
                        for icol in 1..=lscalar.coeff().row_len() {
                            let coeff = lscalar.coeff().get(irow, icol);
                            let v = mat.get(k + irow, j + 1) + coeff.x * vmat.get(icol as usize);
                            mat.set(k + irow, j + 1, v);
                            let v = mat.get(k + irow, n_dimsousmat + j + 1)
                                + coeff.y * vmat.get(icol as usize);
                            mat.set(k + irow, n_dimsousmat + j + 1, v);
                            let v = mat.get(k + irow, 2 * n_dimsousmat + j + 1)
                                + coeff.z * vmat.get(icol as usize);
                            mat.set(k + irow, 2 * n_dimsousmat + j + 1, v);
                        }
                    }

                    j += 1;
                    iv += 1;
                }
            }

            let mut k2 = 3 * n_dimsousmat;
            let mut kppc2 = n_ppc2;
            for i2 in 1..=i {
                let lscalar2 = &self.my_lscalar_constraints[i2 - 1];

                let mut tmpmat = MatD::new(lscalar.get_ppc().len(), lscalar2.get_ppc().len());

                for ippc in 1..=lscalar.get_ppc().len() {
                    for ippc2 in 1..=lscalar2.get_ppc().len() {
                        let deru = self.deru.as_ref().expect("deru");
                        let derv = self.derv.as_ref().expect("derv");
                        let points = self.points.as_ref().expect("points");
                        let mut signe = 1.0f64;
                        if (deru[kppc2 + ippc2 - 1] + derv[kppc2 + ippc2 - 1]) % 2 == 1 {
                            signe = -1.0;
                        }
                        let a_iu = deru[kppc + ippc - 1] + deru[kppc2 + ippc2 - 1];
                        let iv = derv[kppc + ippc - 1] + derv[kppc2 + ippc2 - 1];
                        tmpmat.set(
                            ippc as usize,
                            ippc2 as usize,
                            signe
                                * self.sol_em(
                                    points[kppc + ippc - 1] - points[kppc2 + ippc2 - 1],
                                    a_iu,
                                    iv,
                                ),
                        );
                    }
                }

                for irow in 1..=lscalar.coeff().col_len() {
                    for irow2 in 1..=lscalar2.coeff().col_len() {
                        for icol in 1..=lscalar.coeff().row_len() {
                            for icol2 in 1..=lscalar2.coeff().row_len() {
                                // OCCT: Coeff()(irow, icol) * Coeff2()(irow2, icol2)
                                // is the gp_XYZ dot product.
                                let dot = lscalar.coeff().get(irow, icol).x
                                    * lscalar2.coeff().get(irow2, icol2).x
                                    + lscalar.coeff().get(irow, icol).y
                                        * lscalar2.coeff().get(irow2, icol2).y
                                    + lscalar.coeff().get(irow, icol).z
                                        * lscalar2.coeff().get(irow2, icol2).z;
                                let v = mat.get(k + irow, k2 + irow2)
                                    + dot * tmpmat.get(icol as usize, icol2 as usize);
                                mat.set(k + irow, k2 + irow2, v);
                            }
                        }
                    }
                }

                k2 += lscalar2.coeff().col_len();
                kppc2 += lscalar2.coeff().row_len();
            }

            k += lscalar.coeff().col_len();
            kppc += lscalar.coeff().row_len();
        }

        for j in (3 * n_dimsousmat)..n_dimat {
            for i in 0..j {
                let v = mat.get(j + 1, i + 1);
                mat.set(i + 1, j + 1, v);
            }
        }

        // initialisation of the Gauss algorithm
        let mut pivot_max = 1.0e-12;
        self.ok = true; // ************ JHH

        // OCCT: math_Gauss algo_gauss(mat, pivot_max, aScope.Next(7));
        let mut algo_gauss = MathGauss::with_min_pivot(&mat, pivot_max);

        if !algo_gauss.is_done() {
            for i in (n_cc1 + n_cc2)..(n_cc1 + n_cc2 + nbm) {
                mat.set(i + 1, i + 1, 1.0e-8);
                mat.set(n_dimsousmat + i + 1, n_dimsousmat + i + 1, 1.0e-8);
                mat.set(2 * n_dimsousmat + i + 1, 2 * n_dimsousmat + i + 1, 1.0e-8);
            }
            pivot_max = 1.0e-18;

            let the_algo2 = MathGauss::with_min_pivot(&mat, pivot_max);
            algo_gauss = the_algo2;
            self.ok = algo_gauss.is_done();
        }

        if self.ok {
            //   computation of the linear system solution for the X, Y and Z
            //   coordinates
            let mut sec_member = VecD::new(n_dimat);
            let mut sol = VecD::new(n_dimat);

            self.n_dim = self.n_el + self.order * (self.order + 1) / 2;
            self.solution = Some(vec![DVec3::ZERO; self.n_dim as usize]);

            for icoor in 1..=3 {
                for i in 0..n_cc1 {
                    sec_member.set(
                        ((icoor - 1) * n_dimsousmat as i32 + i as i32 + 1) as usize,
                        xyz_coord(self.my_constraints[i].value(), icoor),
                    );
                }

                let mut k = n_cc1;
                for lxyz in &self.my_lxyz_constraints {
                    for irow in 1..=lxyz.coeff().col_len() {
                        for icol in 1..=lxyz.coeff().row_len() {
                            let v = sec_member.get((icoor as usize - 1) * n_dimsousmat + k + 1)
                                + lxyz.coeff().get(irow, icol)
                                    * xyz_coord(lxyz.get_ppc()[icol - 1].value(), icoor);
                            sec_member.set(
                                (icoor as usize - 1) * n_dimsousmat + k + 1,
                                v,
                            );
                        }
                        k += 1;
                    }
                }
            }
            let mut k = 3 * n_dimsousmat;
            for lscalar in &self.my_lscalar_constraints {
                for irow in 1..=lscalar.coeff().col_len() {
                    for icol in 1..=lscalar.coeff().row_len() {
                        // OCCT: Coeff()(irow, icol) * GetPPC()(icol).Value() is
                        // the gp_XYZ dot product.
                        let coeff = lscalar.coeff().get(irow, icol);
                        let value = lscalar.get_ppc()[icol - 1].value();
                        let v = sec_member.get(k + 1)
                            + coeff.x * value.x
                            + coeff.y * value.y
                            + coeff.z * value.z;
                        sec_member.set(k + 1, v);
                    }
                    k += 1;
                }
            }

            let mut sol_ = sec_member.clone();
            algo_gauss.solve(&mut sol_);
            sol = sol_;
            // iteration to refine the solution
            {
                let mut sol1 = VecD::new(n_dimat);
                let mut sec_member1 = VecD::new(n_dimat);
                for _ in 1..=iteration_number {
                    sec_member1 = vec_sub(&sec_member, &mat_vec_mul(&mat, &sol));
                    sol1 = sec_member1.clone();
                    algo_gauss.solve(&mut sol1);
                    vec_add(&mut sol, &sol1);
                }
            }

            for icoor in 1..=3 {
                for i in 0..n_cc1 {
                    let v = sol.get(((icoor - 1) * n_dimsousmat as i32 + i as i32 + 1) as usize);
                    let solution = self.solution.as_mut().expect("solution");
                    xyz_set_coord(&mut solution[i], icoor, v);
                }

                let mut k_solution = n_cc1;
                let mut ksol = n_cc1;

                for lxyz in &self.my_lxyz_constraints {
                    for icol in 1..=lxyz.coeff().row_len() {
                        let mut vsol = 0.0f64;
                        for irow in 1..=lxyz.coeff().col_len() {
                            vsol += lxyz.coeff().get(irow, icol)
                                * sol.get(
                                    ((icoor - 1) * n_dimsousmat as i32
                                        + ksol as i32
                                        + irow as i32) as usize,
                                );
                        }
                        let solution = self.solution.as_mut().expect("solution");
                        xyz_set_coord(&mut solution[k_solution], icoor, vsol);
                        k_solution += 1;
                    }
                    ksol += lxyz.coeff().col_len();
                }

                ksol = n_cc1 + n_cc2;
                for i in 0..(self.order * (self.order + 1) / 2) as usize {
                    let v = sol.get(((icoor - 1) * n_dimsousmat as i32 + ksol as i32 + i as i32 + 1)
                        as usize);
                    let solution = self.solution.as_mut().expect("solution");
                    xyz_set_coord(&mut solution[self.n_el as usize + i], icoor, v);
                }
            }

            let mut ksol = 3 * n_dimsousmat;
            let mut k_solution = n_ppc2;
            for lscalar in &self.my_lscalar_constraints {
                for icol in 1..=lscalar.coeff().row_len() {
                    let mut vsol = DVec3::ZERO;
                    for irow in 1..=lscalar.coeff().col_len() {
                        vsol += lscalar.coeff().get(irow, icol)
                            * sol.get(ksol + irow);
                    }
                    let solution = self.solution.as_mut().expect("solution");
                    solution[k_solution] = vsol;
                    k_solution += 1;
                }
                ksol += lscalar.coeff().col_len();
            }
        }
    }

    /// OCCT fillXYZmatrix (Plate_Plate.cxx L1060-1213).
    fn fill_xyz_matrix(&self, mat: &mut MatD, i0: i32, j0: i32, ncc1: usize, ncc2: usize) {
        let deru = self.deru.as_ref().expect("deru");
        let derv = self.derv.as_ref().expect("derv");
        let points = self.points.as_ref().expect("points");

        let i0u = i0 as usize;
        let j0u = j0 as usize;

        for i in 0..ncc1 {
            for j in 0..i {
                let mut signe = 1.0f64;
                if (deru[j] + derv[j]) % 2 == 1 {
                    signe = -1.0;
                }
                let iu = deru[i] + deru[j];
                let iv = derv[i] + derv[j];
                mat.set(
                    i0u + i + 1,
                    j0u + j + 1,
                    signe * self.sol_em(points[i] - points[j], iu, iv),
                );
            }
        }

        let mut k = ncc1;
        let mut kppc = ncc1;
        for (lxyz_idx, lxyz) in self.my_lxyz_constraints.iter().enumerate() {
            for a_j in 0..ncc1 {
                let mut vmat = VecD::new(lxyz.get_ppc().len());

                for ippc in 1..=lxyz.get_ppc().len() {
                    let mut signe = 1.0f64;
                    if (deru[a_j] + derv[a_j]) % 2 == 1 {
                        signe = -1.0;
                    }
                    let iu = deru[kppc + ippc - 1] + deru[a_j];
                    let iv = derv[kppc + ippc - 1] + derv[a_j];
                    vmat.set(
                        ippc as usize,
                        signe * self.sol_em(points[kppc + ippc - 1] - points[a_j], iu, iv),
                    );
                }

                for irow in 1..=lxyz.coeff().col_len() {
                    for icol in 1..=lxyz.coeff().row_len() {
                        let v = mat.get(i0u + k + irow, j0u + a_j + 1)
                            + lxyz.coeff().get(irow, icol) * vmat.get(icol as usize);
                        mat.set(i0u + k + irow, j0u + a_j + 1, v);
                    }
                }
            }

            let mut k2 = ncc1;
            let mut kppc2 = ncc1;
            for i2 in 1..=lxyz_idx + 1 {
                let lxyz2 = &self.my_lxyz_constraints[i2 - 1];

                let mut tmpmat = MatD::new(lxyz.get_ppc().len(), lxyz2.get_ppc().len());

                for ippc in 1..=lxyz.get_ppc().len() {
                    for ippc2 in 1..=lxyz2.get_ppc().len() {
                        let mut signe = 1.0f64;
                        if (deru[kppc2 + ippc2 - 1] + derv[kppc2 + ippc2 - 1]) % 2 == 1 {
                            signe = -1.0;
                        }
                        let iu = deru[kppc + ippc - 1] + deru[kppc2 + ippc2 - 1];
                        let iv = derv[kppc + ippc - 1] + derv[kppc2 + ippc2 - 1];
                        tmpmat.set(
                            ippc as usize,
                            ippc2 as usize,
                            signe
                                * self.sol_em(
                                    points[kppc + ippc - 1] - points[kppc2 + ippc2 - 1],
                                    iu,
                                    iv,
                                ),
                        );
                    }
                }

                for irow in 1..=lxyz.coeff().col_len() {
                    for irow2 in 1..=lxyz2.coeff().col_len() {
                        for icol in 1..=lxyz.coeff().row_len() {
                            for icol2 in 1..=lxyz2.coeff().row_len() {
                                let v = mat.get(i0u + k + irow, j0u + k2 + irow2)
                                    + lxyz.coeff().get(irow, icol)
                                        * lxyz2.coeff().get(irow2, icol2)
                                        * tmpmat.get(icol as usize, icol2 as usize);
                                mat.set(i0u + k + irow, j0u + k2 + irow2, v);
                            }
                        }
                    }
                }

                k2 += lxyz2.coeff().col_len();
                kppc2 += lxyz2.coeff().row_len();
            }

            k += lxyz.coeff().col_len();
            kppc += lxyz.coeff().row_len();
        }

        let mut i = ncc1 + ncc2;
        for iu in 0..self.order {
            let mut iv = 0i32;
            while iu + iv < self.order {
                for a_j in 0..ncc1 {
                    let idu = deru[a_j];
                    let idv = derv[a_j];
                    mat.set(
                        i0u + i + 1,
                        j0u + a_j + 1,
                        self.polm(points[a_j], iu, iv, idu, idv),
                    );
                }

                let mut k2 = ncc1;
                let mut kppc2 = ncc1;
                for lxyz2 in &self.my_lxyz_constraints {
                    let mut vmat = VecD::new(lxyz2.get_ppc().len());
                    for ippc2 in 1..=lxyz2.get_ppc().len() {
                        let idu = deru[kppc2 + ippc2 - 1];
                        let idv = derv[kppc2 + ippc2 - 1];
                        vmat.set(
                            ippc2 as usize,
                            self.polm(points[kppc2 + ippc2 - 1], iu, iv, idu, idv),
                        );
                    }

                    for irow2 in 1..=lxyz2.coeff().col_len() {
                        for icol2 in 1..=lxyz2.coeff().row_len() {
                            let v = mat.get(i0u + i + 1, j0u + k2 + irow2)
                                + lxyz2.coeff().get(irow2, icol2) * vmat.get(icol2 as usize);
                            mat.set(i0u + i + 1, j0u + k2 + irow2, v);
                        }
                    }

                    k2 += lxyz2.coeff().col_len();
                    kppc2 += lxyz2.coeff().row_len();
                }

                i += 1;
                iv += 1;
            }
        }

        let n_dimat = ncc1 + ncc2 + (self.order * (self.order + 1) / 2) as usize;

        for i in 0..n_dimat {
            for a_j in (i + 1)..n_dimat {
                let v = mat.get(i0u + a_j + 1, j0u + i + 1);
                mat.set(i0u + i + 1, j0u + a_j + 1, v);
            }
        }
    }

    /// OCCT IsDone (Plate_Plate.cxx L1217-1220).
    pub fn is_done(&self) -> bool {
        self.ok
    }

    /// OCCT destroy (Plate_Plate.cxx L1224-1227) — delegates to Init.
    pub fn destroy(&mut self) {
        self.init();
    }

    /// OCCT Init (Plate_Plate.cxx L1231-1254).  Note OCCT sets OK = **true**
    /// here (finding #30 — kept as-is on purpose).
    pub fn init(&mut self) {
        self.my_constraints.clear();
        self.my_lxyz_constraints.clear();
        self.my_lscalar_constraints.clear();

        self.solution = None;
        self.points = None;
        self.deru = None;
        self.derv = None;

        self.order = 0;
        self.n_el = 0;
        self.n_dim = 0;
        self.ok = true;
        self.max_constraint_order = 0;
    }

    /// OCCT Evaluate (Plate_Plate.cxx L1258-1293).
    pub fn evaluate(&self, point2d: DVec2) -> DVec3 {
        let solution = match &self.solution {
            Some(s) => s,
            None => return DVec3::ZERO,
        };
        if !self.ok {
            return DVec3::ZERO;
        }

        let mut valeur = DVec3::ZERO;

        if !self.polynomial_part_only {
            let deru = self.deru.as_ref().expect("deru");
            let derv = self.derv.as_ref().expect("derv");
            let points = self.points.as_ref().expect("points");
            for i in 0..self.n_el as usize {
                let mut signe = 1.0f64;
                if (deru[i] + derv[i]) % 2 == 1 {
                    signe = -1.0;
                }
                valeur += solution[i]
                    * (signe * self.sol_em(point2d - points[i], deru[i], derv[i]));
            }
        }
        let mut i = self.n_el as usize;
        for idu in 0..self.order {
            let mut idv = 0i32;
            while idu + idv < self.order {
                valeur += solution[i] * self.polm(point2d, idu, idv, 0, 0);
                i += 1;
                idv += 1;
            }
        }
        valeur
    }

    /// OCCT EvaluateDerivative (Plate_Plate.cxx L1297-1331).
    pub fn evaluate_derivative(&self, point2d: DVec2, iu: i32, iv: i32) -> DVec3 {
        let solution = match &self.solution {
            Some(s) => s,
            None => return DVec3::ZERO,
        };
        if !self.ok {
            return DVec3::ZERO;
        }

        let mut valeur = DVec3::ZERO;
        if !self.polynomial_part_only {
            let deru = self.deru.as_ref().expect("deru");
            let derv = self.derv.as_ref().expect("derv");
            let points = self.points.as_ref().expect("points");
            for i in 0..self.n_el as usize {
                let mut signe = 1.0f64;
                if (deru[i] + derv[i]) % 2 == 1 {
                    signe = -1.0;
                }
                valeur += solution[i]
                    * (signe * self.sol_em(point2d - points[i], deru[i] + iu, derv[i] + iv));
            }
        }
        let mut i = self.n_el as usize;
        for idu in 0..self.order {
            let mut idv = 0i32;
            while idu + idv < self.order {
                valeur += solution[i] * self.polm(point2d, idu, idv, iu, iv);
                i += 1;
                idv += 1;
            }
        }
        valeur
    }

    /// OCCT CoefPol (Plate_Plate.cxx L1339-1355) — power-basis coefficients of
    /// the polynomial part; 0-based [iu][iv] grid of size order x order
    /// (OCCT HArray2(0, order-1, 0, order-1)).
    pub fn coef_pol(&self) -> Vec<Vec<DVec3>> {
        let order = self.order as usize;
        let mut coefs = vec![vec![DVec3::ZERO; order]; order];
        let mut i = self.n_el as usize;
        for iu in 0..self.order as usize {
            for iv in 0..(self.order as usize - iu) {
                coefs[iu][iv] = self.solution.as_ref().expect("solution")[i]
                    * self.ddu[iu]
                    * self.ddv[iv];
                i += 1;
            }
        }
        coefs
    }

    /// OCCT Continuity (Plate_Plate.cxx L1367-1370).
    pub fn continuity(&self) -> i32 {
        2 * self.order - 3 - self.max_constraint_order
    }

    /// OCCT SolEm (Plate_Plate.cxx L1378-1886) — (iu,iv)th derivative of the
    /// fundamental solution of the Laplacian at the power order.
    fn sol_em(&self, point2d: DVec2, iu: i32, iv: i32) -> f64 {
        let u;
        let v;
        let i_u;
        let i_v;

        if iv > iu {
            // SolEm is symmetric in (u<->v) : we swap u and v if iv>iu
            // to avoid some code
            i_u = iv;
            i_v = iu;
            u = point2d.y * self.ddv[1];
            v = point2d.x * self.ddu[1];
        } else {
            i_u = iu;
            i_v = iv;
            u = point2d.x * self.ddu[1];
            v = point2d.y * self.ddv[1];
        }

        if (u == self.u_old.get()) && (v == self.v_old.get()) {
            if self.r.get() < 1.0e-20 {
                return 0.0;
            }
        } else {
            self.u_old.set(u);
            self.v_old.set(v);
            self.u2.set(u * u);
            self.r.set(self.u2.get() + v * v);
            if self.r.get() < 1.0e-20 {
                return 0.0;
            }
            self.l.set(self.r.get().ln());
        }
        let mut d_uv = 0.0f64;

        let m = self.order;
        let mm1 = m - 1;
        let r = self.r.get();
        let u2 = self.u2.get();
        let l = self.l.get();

        // pr = pow(R, mm1 - IU - IV) with a small integer exponent.
        let expo = mm1 - i_u - i_v;
        let pr;
        if expo < 0 {
            let mut prr = r;
            for _ in 1..(-expo) {
                prr *= r;
            }
            pr = 1.0 / prr;
        } else if expo > 0 {
            let mut prr = r;
            for _ in 1..expo {
                prr *= r;
            }
            pr = prr;
        } else {
            pr = 1.0;
        }

        match i_u {
            0 => match i_v {
                0 => {
                    d_uv = pr * l;
                }
                _ => {}
            },
            1 => match i_v {
                0 => {
                    d_uv = 2.0 * pr * u * (1.0 + l * mm1 as f64);
                }
                1 => {
                    let m2 = (m * m) as f64;
                    // DUV = 4*pr*U*V*(-3+2*L+2*m-3*L*m+L*m2);
                    d_uv = 4.0 * pr * u * v * ((2.0 * m as f64 - 3.0) + (m2 - 3.0 * m as f64 + 2.0) * l);
                }
                _ => {}
            },
            2 => match i_v {
                0 => {
                    let m2 = (m * m) as f64;
                    d_uv = 2.0 * pr
                        * (r - l * r
                            + l * m as f64 * r
                            - 6.0 * u2
                            + 4.0 * l * u2
                            + 4.0 * m as f64 * u2
                            - 6.0 * l * m as f64 * u2
                            + 2.0 * l * m2 * u2);
                }
                1 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    d_uv = -3.0 * r + 2.0 * l * r + 2.0 * m as f64 * r - 3.0 * l * m as f64 * r
                        + l * m2 * r
                        + 22.0 * u2
                        - 12.0 * l * u2
                        - 24.0 * m as f64 * u2
                        + 22.0 * l * m as f64 * u2
                        + 6.0 * m2 * u2
                        - 12.0 * l * m2 * u2
                        + 2.0 * l * m3 * u2;
                    d_uv = d_uv * 4.0 * pr * v;
                }
                2 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    let v2 = v * v;
                    let r2 = r * r;
                    d_uv = -3.0 * r2 + 2.0 * l * r2 + 2.0 * m as f64 * r2 - 3.0 * l * m as f64 * r2
                        + l * m2 * r2
                        + 22.0 * r * u2
                        - 12.0 * l * r * u2
                        - 24.0 * m as f64 * r * u2
                        + 22.0 * l * m as f64 * r * u2
                        + 6.0 * m2 * r * u2
                        - 12.0 * l * m2 * r * u2;
                    d_uv += 2.0 * l * m3 * r * u2
                        + 22.0 * r * v2
                        - 12.0 * l * r * v2
                        - 24.0 * m as f64 * r * v2
                        + 22.0 * l * m as f64 * r * v2
                        + 6.0 * m2 * r * v2
                        - 12.0 * l * m2 * r * v2
                        + 2.0 * l * m3 * r * v2
                        - 200.0 * u2 * v2
                        + 96.0 * l * u2 * v2;
                    d_uv += 280.0 * m as f64 * u2 * v2
                        - 200.0 * l * m as f64 * u2 * v2
                        - 120.0 * m2 * u2 * v2
                        + 140.0 * l * m2 * u2 * v2
                        + 16.0 * m3 * u2 * v2
                        - 40.0 * l * m3 * u2 * v2
                        + 4.0 * l * m4 * u2 * v2;
                    d_uv = 4.0 * pr * d_uv;
                }
                _ => {}
            },
            3 => match i_v {
                0 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    d_uv = -9.0 * r + 6.0 * l * r + 6.0 * m as f64 * r - 9.0 * l * m as f64 * r
                        + 3.0 * l * m2 * r
                        + 22.0 * u2
                        - 12.0 * l * u2
                        - 24.0 * m as f64 * u2
                        + 22.0 * l * m as f64 * u2
                        + 6.0 * m2 * u2
                        - 12.0 * l * m2 * u2
                        + 2.0 * l * m3 * u2;
                    d_uv = d_uv * 4.0 * pr * u;
                }
                1 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    d_uv = 33.0 * r - 18.0 * l * r - 36.0 * m as f64 * r + 33.0 * l * m as f64 * r
                        + 9.0 * m2 * r
                        - 18.0 * l * m2 * r
                        + 3.0 * l * m3 * r
                        - 100.0 * u2
                        + 48.0 * l * u2
                        + 140.0 * m as f64 * u2
                        - 100.0 * l * m as f64 * u2
                        - 60.0 * m2 * u2
                        + 70.0 * l * m2 * u2;
                    d_uv += 8.0 * m3 * u2 - 20.0 * l * m3 * u2 + 2.0 * l * m4 * u2;
                    d_uv = 8.0 * pr * u * v * d_uv;
                }
                2 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    let m5 = m4 * m as f64;
                    let ru2 = r * u2;
                    let v2 = v * v;
                    let rv2 = r * v2;
                    let u2v2 = v2 * u2;
                    let r2 = r * r;

                    // copy-paste the mathematics
                    d_uv = -100.0 * ru2 + 48.0 * l * ru2 + 140.0 * m as f64 * ru2
                        - 100.0 * l * m as f64 * ru2
                        - 60.0 * m2 * ru2
                        + 70.0 * l * m2 * ru2
                        + 8.0 * m3 * ru2
                        - 20.0 * l * m3 * ru2
                        + 2.0 * l * m4 * ru2
                        - 300.0 * rv2
                        + 144.0 * l * rv2
                        + 420.0 * m as f64 * rv2
                        - 300.0 * l * m as f64 * rv2
                        - 180.0 * m2 * rv2
                        + 210.0 * l * m2 * rv2
                        + 24.0 * m3 * rv2
                        - 60.0 * l * m3 * rv2
                        + 6.0 * l * m4 * rv2
                        + 33.0 * r2
                        - 18.0 * l * r2
                        - 36.0 * m as f64 * r2
                        + 33.0 * l * m as f64 * r2
                        + 9.0 * m2 * r2
                        - 18.0 * l * m2 * r2
                        + 3.0 * l * m3 * r2
                        + 1096.0 * u2v2
                        - 480.0 * l * u2v2
                        - 1800.0 * m as f64 * u2v2
                        + 1096.0 * l * m as f64 * u2v2
                        + 1020.0 * m2 * u2v2
                        - 900.0 * l * m2 * u2v2
                        - 240.0 * m3 * u2v2
                        + 340.0 * l * m3 * u2v2
                        + 20.0 * m4 * u2v2
                        - 60.0 * l * m4 * u2v2
                        + 4.0 * l * m5 * u2v2;

                    d_uv = 8.0 * pr * u * d_uv;
                }
                3 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    let m5 = m3 * m as f64;
                    let m6 = m3 * m as f64;
                    let ru2 = r * u2;
                    let v2 = v * v;
                    let rv2 = r * v2;
                    let u2v2 = v2 * u2;
                    let r2 = r * r;

                    d_uv = 1644.0 * ru2 - 720.0 * l * ru2 - 2700.0 * m as f64 * ru2
                        + 1644.0 * l * m as f64 * ru2
                        + 1530.0 * m2 * ru2
                        - 1350.0 * l * m2 * ru2
                        - 360.0 * m3 * ru2
                        + 510.0 * l * m3 * ru2
                        + 30.0 * m4 * ru2
                        - 90.0 * l * m4 * ru2
                        + 6.0 * l * m5 * ru2
                        + 1644.0 * rv2
                        - 720.0 * l * rv2
                        - 2700.0 * m as f64 * rv2
                        + 1644.0 * l * m as f64 * rv2
                        + 1530.0 * m2 * rv2
                        - 1350.0 * l * m2 * rv2
                        - 360.0 * m3 * rv2
                        + 510.0 * l * m3 * rv2
                        + 30.0 * m4 * rv2
                        - 90.0 * l * m4 * rv2
                        + 6.0 * l * m5 * rv2
                        - 450.0 * r2
                        + 216.0 * l * r2
                        + 630.0 * m as f64 * r2
                        - 450.0 * l * m as f64 * r2
                        - 270.0 * m2 * r2
                        + 315.0 * l * m2 * r2
                        + 36.0 * m3 * r2
                        - 90.0 * l * m3 * r2
                        + 9.0 * l * m4 * r2
                        - 7056.0 * u2v2
                        + 2880.0 * l * u2v2
                        + 12992.0 * m as f64 * u2v2
                        - 7056.0 * l * m as f64 * u2v2
                        - 8820.0 * m2 * u2v2
                        + 6496.0 * l * m2 * u2v2
                        + 2800.0 * m3 * u2v2
                        - 2940.0 * l * m3 * u2v2
                        - 420.0 * m4 * u2v2
                        + 700.0 * l * m4 * u2v2
                        + 24.0 * m5 * u2v2
                        - 84.0 * l * m5 * u2v2
                        + 4.0 * l * m6 * u2v2;

                    d_uv = 16.0 * pr * u * v * d_uv;
                }
                _ => {}
            },
            4 => match i_v {
                0 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    let u4 = u2 * u2;
                    let r2 = r * r;
                    d_uv = -9.0 * r2 + 6.0 * l * r2 + 6.0 * m as f64 * r2 - 9.0 * l * m as f64 * r2
                        + 3.0 * l * m2 * r2
                        + 132.0 * r * u2
                        - 72.0 * l * r * u2
                        - 144.0 * m as f64 * r * u2
                        + 132.0 * l * m as f64 * r * u2
                        + 36.0 * m2 * r * u2
                        - 72.0 * l * m2 * r * u2;
                    d_uv += 12.0 * l * m3 * r * u2 - 200.0 * u4 + 96.0 * l * u4
                        + 280.0 * m as f64 * u4
                        - 200.0 * l * m as f64 * u4
                        - 120.0 * m2 * u4
                        + 140.0 * l * m2 * u4
                        + 16.0 * m3 * u4
                        - 40.0 * l * m3 * u4
                        + 4.0 * l * m4 * u4;
                    d_uv = 4.0 * pr * d_uv;
                }
                1 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    let m5 = m2 * m3;
                    let u4 = u2 * u2;
                    let ru2 = r * u2;
                    let r2 = r * r;

                    d_uv = -600.0 * ru2 + 288.0 * l * ru2 + 840.0 * m as f64 * ru2
                        - 600.0 * l * m as f64 * ru2
                        - 360.0 * m2 * ru2
                        + 420.0 * l * m2 * ru2
                        + 48.0 * m3 * ru2
                        - 120.0 * l * m3 * ru2
                        + 12.0 * l * m4 * ru2
                        + 33.0 * r2
                        - 18.0 * l * r2
                        - 36.0 * m as f64 * r2
                        + 33.0 * l * m as f64 * r2
                        + 9.0 * m2 * r2
                        - 18.0 * l * m2 * r2
                        + 3.0 * l * m3 * r2
                        + 1096.0 * u4
                        - 480.0 * l * u4
                        - 1800.0 * m as f64 * u4
                        + 1096.0 * l * m as f64 * u4
                        + 1020.0 * m2 * u4
                        - 900.0 * l * m2 * u4
                        - 240.0 * m3 * u4
                        + 340.0 * l * m3 * u4
                        + 20.0 * m4 * u4
                        - 60.0 * l * m4 * u4
                        + 4.0 * l * m5 * u4;

                    d_uv = 8.0 * pr * v * d_uv;
                }
                2 => {
                    let m2 = (m * m) as f64;
                    let m3 = m2 * m as f64;
                    let m4 = m2 * m2;
                    let m5 = m2 * m3;
                    let m6 = m3 * m as f64;
                    let u4 = u2 * u2;
                    let r2 = r * r;
                    let r3 = r2 * r;
                    let v2 = v * v;
                    let u2v2 = v2 * u2;
                    let ru2v2 = r * u2v2;
                    let u4v2 = u4 * v2;
                    let r2u2 = r2 * u2;
                    let ru4 = r * u4;
                    let r2v2 = r2 * v2;

                    d_uv = 6576.0 * ru2v2 - 2880.0 * l * ru2v2 - 10800.0 * m as f64 * ru2v2
                        + 6576.0 * l * m as f64 * ru2v2
                        + 6120.0 * m2 * ru2v2
                        - 5400.0 * l * m2 * ru2v2
                        - 1440.0 * m3 * ru2v2
                        + 2040.0 * l * m3 * ru2v2
                        + 120.0 * m4 * ru2v2
                        - 360.0 * l * m4 * ru2v2
                        + 24.0 * l * m5 * ru2v2
                        + 1096.0 * ru4
                        - 480.0 * l * ru4
                        - 1800.0 * m as f64 * ru4
                        + 1096.0 * l * m as f64 * ru4
                        + 1020.0 * m2 * ru4
                        - 900.0 * l * m2 * ru4
                        - 240.0 * m3 * ru4
                        + 340.0 * l * m3 * ru4
                        + 20.0 * m4 * ru4
                        - 60.0 * l * m4 * ru4
                        + 4.0 * l * m5 * ru4
                        - 600.0 * r2u2
                        + 288.0 * l * r2u2
                        + 840.0 * m as f64 * r2u2
                        - 600.0 * l * m as f64 * r2u2
                        - 360.0 * m2 * r2u2
                        + 420.0 * l * m2 * r2u2
                        + 48.0 * m3 * r2u2
                        - 120.0 * l * m3 * r2u2
                        + 12.0 * l * m4 * r2u2
                        - 300.0 * r2v2
                        + 144.0 * l * r2v2
                        + 420.0 * m as f64 * r2v2
                        - 300.0 * l * m as f64 * r2v2
                        - 180.0 * m2 * r2v2
                        + 210.0 * l * m2 * r2v2
                        + 24.0 * m3 * r2v2
                        - 60.0 * l * m3 * r2v2
                        + 6.0 * l * m4 * r2v2
                        + 33.0 * r3
                        - 18.0 * l * r3
                        - 36.0 * m as f64 * r3
                        + 33.0 * l * m as f64 * r3
                        + 9.0 * m2 * r3
                        - 18.0 * l * m2 * r3
                        + 3.0 * l * m3 * r3
                        - 14112.0 * u4v2
                        + 5760.0 * l * u4v2
                        + 25984.0 * m as f64 * u4v2
                        - 14112.0 * l * m as f64 * u4v2
                        - 17640.0 * m2 * u4v2
                        + 12992.0 * l * m2 * u4v2
                        + 5600.0 * m3 * u4v2
                        - 5880.0 * l * m3 * u4v2
                        - 840.0 * m4 * u4v2
                        + 11760.0 * l * m4 * u4v2
                        + 408.0 * m5 * u4v2
                        - 1680.0 * l * m5 * u4v2
                        + 120.0 * m6 * u4v2;
                    d_uv = 4.0 * pr * u * d_uv;
                }
                _ => {}
            },
            _ => {}
        }

        d_uv
    }

    /// OCCT Polm (Plate_Plate.lxx L25-56).
    fn polm(&self, point2d: DVec2, iu: i32, iv: i32, idu: i32, idv: i32) -> f64 {
        if idu > iu {
            return 0.0;
        }
        if idv > iv {
            return 0.0;
        }
        let u = point2d.x;
        let v = point2d.y;

        let mut value = 1.0f64;

        let degu = iu - idu;
        for _i in 0..degu {
            value *= u;
        }
        for i in (degu + 1)..=iu {
            value *= i as f64;
        }

        let degv = iv - idv;
        for _i in 0..degv {
            value *= v;
        }
        for i in (degv + 1)..=iv {
            value *= i as f64;
        }

        // le produit par ddu[iu]*ddv[iv] n'est pas indispensable !! (il change
        // les valeurs calculees pour la partie coef polynomiaux de Sol
        // de telle facon que les methodes Evaluate et EvaluateDerivative donnent
        // en theorie les memes valeurs. Toutefois, il nous semble que ce produit
        // ameliore le conditionnement de la matrice
        value * self.ddu[iu as usize] * self.ddv[iv as usize]
    }

    /// OCCT UVBox (Plate_Plate.cxx L1888-1953).
    pub fn uv_box(&self, u_min: &mut f64, u_max: &mut f64, v_min: &mut f64, v_max: &mut f64) {
        const BMIN: f64 = 1.0e-3;
        *u_min = self.my_constraints[0].pnt2d().x;
        *u_max = *u_min;
        *v_min = self.my_constraints[0].pnt2d().y;
        *v_max = *v_min;

        for c in &self.my_constraints {
            let x = c.pnt2d().x;
            if x < *u_min {
                *u_min = x;
            }
            if x > *u_max {
                *u_max = x;
            }
            let y = c.pnt2d().y;
            if y < *v_min {
                *v_min = y;
            }
            if y > *v_max {
                *v_max = y;
            }
        }

        for lxyz in &self.my_lxyz_constraints {
            for ppc in lxyz.get_ppc() {
                let x = ppc.pnt2d().x;
                if x < *u_min {
                    *u_min = x;
                }
                if x > *u_max {
                    *u_max = x;
                }
                let y = ppc.pnt2d().y;
                if y < *v_min {
                    *v_min = y;
                }
                if y > *v_max {
                    *v_max = y;
                }
            }
        }

        for lscalar in &self.my_lscalar_constraints {
            for ppc in lscalar.get_ppc() {
                let x = ppc.pnt2d().x;
                if x < *u_min {
                    *u_min = x;
                }
                if x > *u_max {
                    *u_max = x;
                }
                let y = ppc.pnt2d().y;
                if y < *v_min {
                    *v_min = y;
                }
                if y > *v_max {
                    *v_max = y;
                }
            }
        }

        if *u_max - *u_min < BMIN {
            let um = 0.5 * (*u_min + *u_max);
            *u_min = um - 0.5 * BMIN;
            *u_max = um + 0.5 * BMIN;
        }
        if *v_max - *v_min < BMIN {
            let vm = 0.5 * (*v_min + *v_max);
            *v_min = vm - 0.5 * BMIN;
            *v_max = vm + 0.5 * BMIN;
        }
    }

    /// OCCT UVConstraints (Plate_Plate.cxx L1985-1996) — appends the 2D points
    /// of the G0 pinpoint constraints.
    pub fn uv_constraints(&self, seq: &mut Vec<DVec2>) {
        for c in &self.my_constraints {
            if c.idu() == 0 && c.idv() == 0 {
                seq.push(c.pnt2d());
            }
        }
    }

    /// OCCT SetPolynomialPartOnly (Plate_Plate.cxx L1998-2001).
    pub fn set_polynomial_part_only(&mut self, pp_only: bool) {
        self.polynomial_part_only = pp_only;
    }
}
