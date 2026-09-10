// The myTolerance member mirrors the OCCT field (read by OCCT callers).
#![allow(dead_code)]
//! OCCT Law_Interpolate (TKGeomAlgo/Law) — 1:1 port of Law_Interpolate.cxx
//! (whole file L26-601) + .hxx members.  Dependencies: BSplCLib::Interpolate
//! (rcad bspl::interpolate) and PLib::EvalLagrange (rcad plib::eval_lagrange).

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::bspl_lib::interpolate as bspl_interpolate;
use rcad_kernel::math::plib::eval_lagrange;

use super::law_bspline::LawBSpline;

// OCCT RealSmall().
const REAL_SMALL: f64 = 2.2250738585072014e-308;

/// OCCT static CheckParameters (Law_Interpolate.cxx L28-41).
fn check_parameters(parameters: &[f64]) -> bool {
    let mut result = true;
    let mut ii = 0usize;
    while result && ii < parameters.len() - 1 {
        let distance = parameters[ii + 1] - parameters[ii];
        result = distance >= REAL_SMALL;
        ii += 1;
    }
    result
}

/// OCCT static BuildParameters (L43-66).
fn build_parameters(periodic_flag: bool, points_array: &[f64]) -> Vec<f64> {
    let mut index = 1usize; // OCCT index starts at 2 (1-based) == slot 2.
    let mut num_parameters = points_array.len();
    if periodic_flag {
        num_parameters += 1;
    }
    let mut parameters = vec![0.0f64; num_parameters];
    parameters[0] = 0.0;
    for ii in 0..points_array.len() - 1 {
        let distance = (points_array[ii] - points_array[ii + 1]).abs();
        parameters[index] = parameters[index - 1] + distance;
        index += 1;
    }
    if periodic_flag {
        let distance = (points_array[points_array.len() - 1] - points_array[0]).abs();
        parameters[index] = parameters[index - 1] + distance;
    }
    parameters
}

/// OCCT static BuildPeriodicTangent (L68-108).
fn build_periodic_tangent(
    points_array: &[f64],
    tangents_array: &mut [f64],
    tangent_flags: &mut [bool],
    parameters_array: &[f64],
) {
    if points_array.len() < 2 {
        tangent_flags[0] = true;
        tangents_array[0] = 0.0;
    } else if !tangent_flags[0] {
        // Pour les periodiques on evalue la tangente du point de fermeture
        // par une interpolation de degre 2 entre le dernier point, le point
        // de fermeture et le deuxieme point.
        let degree = 2usize;
        let period = parameters_array[parameters_array.len() - 1] - parameters_array[0];
        let point_array = [
            points_array[points_array.len() - 1],
            points_array[0],
            points_array[1],
        ];
        let parameter_array = [
            parameters_array[parameters_array.len() - 2] - period,
            parameters_array[0],
            parameters_array[1],
        ];
        tangent_flags[0] = true;
        let mut eval_result = [0.0f64; 2];
        eval_lagrange(
            parameter_array[1],
            1,
            degree,
            1,
            &point_array,
            &parameter_array,
            &mut eval_result,
        );
        tangents_array[0] = eval_result[1];
    }
}

/// OCCT static BuildTangents (L110-160).
fn build_tangents(
    points_array: &[f64],
    tangents_array: &mut [f64],
    tangent_flags: &mut [bool],
    parameters_array: &[f64],
) {
    let mut degree = 3usize;
    if points_array.len() < 3 {
        panic!("Standard_ConstructionError");
    }
    if points_array.len() == 3 {
        degree = 2;
    }
    if !tangent_flags[0] {
        // point_array = &PointsArray(Lower) — the first degree+1 points.
        let point_array: Vec<f64> = points_array[..degree + 1].to_vec();
        let parameter_array: Vec<f64> = parameters_array[..degree + 1].to_vec();
        tangent_flags[0] = true;
        let mut eval_result = [0.0f64; 2];
        eval_lagrange(
            parameters_array[0],
            1,
            degree,
            1,
            &point_array,
            &parameter_array,
            &mut eval_result,
        );
        tangents_array[0] = eval_result[1];
    }
    let upper = tangent_flags.len() - 1;
    if !tangent_flags[upper] {
        // point_array = &PointsArray(Upper - degree) — the last degree+1
        // points.
        let point_array: Vec<f64> = points_array[points_array.len() - 1 - degree..].to_vec();
        tangent_flags[upper] = true;
        let iup = parameters_array.len() - degree;
        let parameter_array: Vec<f64> = parameters_array[iup..].to_vec();
        let mut eval_result = [0.0f64; 2];
        eval_lagrange(
            parameters_array[parameters_array.len() - 1],
            1,
            degree,
            1,
            &point_array,
            &parameter_array,
            &mut eval_result,
        );
        tangents_array[upper] = eval_result[1];
    }
}

/// OCCT Law_Interpolate.
#[derive(Debug, Clone)]
pub struct LawInterpolate {
    my_tolerance: f64,
    my_points: Vec<f64>,
    my_is_done: bool,
    my_parameters: Vec<f64>,
    my_curve: Option<Rc<RefCell<LawBSpline>>>,
    my_periodic: bool,
    my_tangents: Vec<f64>,
    my_tangent_flags: Vec<bool>,
    my_tangent_request: bool,
}

impl LawInterpolate {
    /// OCCT Law_Interpolate(PointsPtr, PeriodicFlag, Tolerance) (L162-180).
    pub fn new(points: Vec<f64>, periodic_flag: bool, tolerance: f64) -> Self {
        let my_parameters = build_parameters(periodic_flag, &points);
        let my_tangent_flags = vec![false; points.len()];
        let my_tangents = vec![0.0; points.len()];
        LawInterpolate {
            my_tolerance: tolerance,
            my_tangent_flags,
            my_is_done: false,
            my_points: points,
            my_parameters,
            my_curve: None,
            my_periodic: periodic_flag,
            my_tangents,
            my_tangent_request: false,
        }
    }

    /// OCCT Law_Interpolate(PointsPtr, ParametersPtr, PeriodicFlag,
    /// Tolerance) (L182-211).
    pub fn with_parameters(
        points: Vec<f64>,
        parameters: Vec<f64>,
        periodic_flag: bool,
        tolerance: f64,
    ) -> Self {
        if periodic_flag && points.len() + 1 != parameters.len() {
            panic!("Standard_ConstructionError");
        }
        if !check_parameters(&parameters) {
            panic!("Standard_ConstructionError");
        }
        let my_tangent_flags = vec![false; points.len()];
        let my_tangents = vec![0.0; points.len()];
        LawInterpolate {
            my_tolerance: tolerance,
            my_tangent_flags,
            my_is_done: false,
            my_points: points,
            my_parameters: parameters,
            my_curve: None,
            my_periodic: periodic_flag,
            my_tangents,
            my_tangent_request: false,
        }
    }

    /// OCCT Load(Tangents, TangentFlagsPtr) (L213-234).
    pub fn load(&mut self, tangents: &[f64], tangent_flags: &[bool]) {
        self.my_tangent_request = true;
        self.my_tangent_flags = tangent_flags.to_vec();
        if tangents.len() != self.my_points.len()
            || tangent_flags.len() != self.my_points.len()
        {
            panic!("Standard_ConstructionError");
        }
        self.my_tangents = tangents.to_vec();
    }

    /// OCCT Load(InitialTangent, FinalTangent) (L236-247).
    pub fn load_end_tangents(&mut self, initial_tangent: f64, final_tangent: f64) {
        self.my_tangent_request = true;
        self.my_tangent_flags[0] = true;
        self.my_tangents[0] = initial_tangent;
        self.my_tangent_flags[self.my_points.len() - 1] = true;
        self.my_tangents[self.my_points.len() - 1] = final_tangent;
    }

    /// OCCT Perform() (L249-256).
    pub fn perform(&mut self) {
        if self.my_periodic {
            self.perform_periodic();
        } else {
            self.perform_non_periodic();
        }
    }

    /// OCCT PerformPeriodic (L258-412).
    fn perform_periodic(&mut self) {
        let num_points = self.my_points.len();
        let period = self.my_parameters[self.my_parameters.len() - 1] - self.my_parameters[0];
        let mut num_poles = num_points + 1;
        let num_distinct_knots = num_points + 1;
        let half_order = 2usize;
        let degree = 3usize;
        num_poles += 2;
        if self.my_tangent_request {
            for flag in &self.my_tangent_flags[1..] {
                if *flag {
                    num_poles += 1;
                }
            }
        }
        let mut parameters = vec![0.0f64; num_poles];
        let mut flatknots = vec![0.0f64; num_poles + degree + 1];
        let mut mults = vec![0i32; num_distinct_knots];
        let mut contact_order_array = vec![0i32; num_poles];
        let mut poles = vec![0.0f64; num_poles];
        let n = self.my_parameters.len();
        for ii in 1..=half_order {
            flatknots[ii - 1] = self.my_parameters[n - 2] - period;
            flatknots[ii + half_order - 1] = self.my_parameters[0];
            flatknots[num_poles + ii - 1] = self.my_parameters[n - 1];
            flatknots[num_poles + half_order + ii - 1] = self.my_parameters[half_order] + period;
        }
        for ii in 1..num_distinct_knots - 1 {
            mults[ii] = 1;
        }
        mults[0] = half_order as i32;
        mults[num_distinct_knots - 1] = half_order as i32;
        {
            let (tangents, flags) = (&mut self.my_tangents, &mut self.my_tangent_flags);
            build_periodic_tangent(&self.my_points, tangents, flags, &self.my_parameters);
        }
        contact_order_array[1] = 1;
        parameters[0] = self.my_parameters[0];
        parameters[1] = self.my_parameters[0];
        poles[0] = self.my_points[0];
        poles[1] = self.my_tangents[0];
        let mut mult_index = 2usize;
        let mut index = 3usize;
        let mut index1 = degree + 2;
        if self.my_tangent_request {
            for ii in 1..self.my_tangent_flags.len() {
                // OCCT: for (ii = Lower + 1; ii <= Upper; ii++) — 1-based ii.
                let ii_1 = ii + 1; // 1-based index into my_parameters
                parameters[index - 1] = self.my_parameters[ii];
                flatknots[index1 - 1] = self.my_parameters[ii];
                poles[index - 1] = self.my_points[ii - 1];
                index += 1;
                index1 += 1;
                if self.my_tangent_flags[ii - 1] {
                    mults[mult_index - 1] += 1;
                    contact_order_array[index - 1] = 1;
                    parameters[index - 1] = self.my_parameters[ii];
                    flatknots[index1 - 1] = self.my_parameters[ii];
                    poles[index - 1] = self.my_tangents[ii - 1];
                    index += 1;
                    index1 += 1;
                }
                mult_index += 1;
                let _ = ii_1;
            }
        } else {
            index = degree + 1;
            index1 = 2;
            for ii in 0..self.my_parameters.len() {
                parameters[index1 - 1] = self.my_parameters[ii];
                flatknots[index - 1] = self.my_parameters[ii];
                index += 1;
                index1 += 1;
            }
            index = 3;
            for ii in 1..self.my_points.len() {
                // copy all the given points since the last one will be
                // initialized below by the first point in myPoints.
                poles[index - 1] = self.my_points[ii];
                index += 1;
            }
        }
        contact_order_array[num_poles - 2] = 1;
        parameters[num_poles - 2] = self.my_parameters[n - 1];
        // for the periodic curve ONLY the tangent of the first point will be
        // used since the curve should close itself at the first point.  See
        // BuildPeriodicTangent.
        poles[num_poles - 2] = self.my_tangents[0];
        parameters[num_poles - 1] = self.my_parameters[n - 1];
        poles[num_poles - 1] = self.my_points[0];
        let inversion_problem = bspl_interpolate(
            degree,
            &flatknots,
            &parameters,
            &contact_order_array,
            1,
            &mut poles,
        );
        if inversion_problem == 0 {
            let newpoles = poles[..num_poles - 2].to_vec();
            self.my_curve = Some(Rc::new(RefCell::new(LawBSpline::new(
                &newpoles,
                &self.my_parameters,
                &mults,
                degree,
                self.my_periodic,
            ))));
            self.my_is_done = true;
        }
    }

    /// OCCT PerformNonPeriodic (L414-560).
    fn perform_non_periodic(&mut self) {
        let num_points = self.my_points.len();
        let mut num_distinct_knots = num_points;
        let mut num_poles = num_points;
        let degree;
        if num_poles == 2 && !self.my_tangent_request {
            degree = 1;
        } else if num_poles == 3 && !self.my_tangent_request {
            degree = 2;
            num_distinct_knots = 2;
        } else {
            degree = 3;
            num_poles += 2;
            if self.my_tangent_request {
                for flag in &self.my_tangent_flags[1..self.my_tangent_flags.len() - 1] {
                    if *flag {
                        num_poles += 1;
                    }
                }
            }
        }
        let mut parameters = vec![0.0f64; num_poles];
        let mut flatknots = vec![0.0f64; num_poles + degree + 1];
        let mut mults = vec![0i32; num_distinct_knots];
        let mut knots = vec![0.0f64; num_distinct_knots];
        let mut contact_order_array = vec![0i32; num_poles];
        let mut poles = vec![0.0f64; num_poles];
        let n = self.my_parameters.len();
        for ii in 1..=degree + 1 {
            flatknots[ii - 1] = self.my_parameters[0];
            flatknots[ii + num_poles - 1] = self.my_parameters[n - 1];
        }
        for ii in 1..num_distinct_knots - 1 {
            mults[ii] = 1;
        }
        mults[0] = degree as i32 + 1;
        mults[num_distinct_knots - 1] = degree as i32 + 1;
        match degree {
            1 => {
                for ii in 0..num_poles {
                    poles[ii] = self.my_points[ii];
                }
                self.my_curve = Some(Rc::new(RefCell::new(LawBSpline::new(
                    &poles,
                    &self.my_parameters,
                    &mults,
                    degree,
                    false,
                ))));
                self.my_is_done = true;
            }
            2 => {
                knots[0] = self.my_parameters[0];
                knots[1] = self.my_parameters[num_poles - 1];
                for ii in 0..num_poles {
                    poles[ii] = self.my_points[ii];
                }
                let inversion_problem = bspl_interpolate(
                    degree,
                    &flatknots,
                    &self.my_parameters,
                    &contact_order_array,
                    1,
                    &mut poles,
                );
                if inversion_problem == 0 {
                    self.my_curve = Some(Rc::new(RefCell::new(LawBSpline::new(
                        &poles, &knots, &mults, degree, false,
                    ))));
                    self.my_is_done = true;
                }
            }
            _ => {
                // check if the boundary conditions are set
                if num_points >= 3 {
                    // cannot build the tangents with degree 3 with only 2
                    // points if those were not given in advance.
                    let (tangents, flags) = (&mut self.my_tangents, &mut self.my_tangent_flags);
                    build_tangents(&self.my_points, tangents, flags, &self.my_parameters);
                }
                contact_order_array[1] = 1;
                parameters[0] = self.my_parameters[0];
                parameters[1] = self.my_parameters[0];
                poles[0] = self.my_points[0];
                poles[1] = self.my_tangents[0];
                let mut mult_index = 2usize;
                let mut index = 3usize;
                let mut index1 = 2usize;
                let mut index2 = 1usize; // OCCT myPoints->Lower() + 1 (1-based)
                let mut index3 = degree + 2;
                if self.my_tangent_request {
                    let params_len = self.my_parameters.len();
                    for ii in 1..params_len - 1 {
                        // OCCT: for (ii = Lower + 1; ii < Upper; ii++)
                        parameters[index - 1] = self.my_parameters[ii];
                        poles[index - 1] = self.my_points[index2 - 1];
                        flatknots[index3 - 1] = self.my_parameters[ii];
                        index += 1;
                        index3 += 1;
                        if self.my_tangent_flags[index1 - 1] {
                            // set the multiplicities, the order of the
                            // contact, the flatknots.
                            mults[mult_index - 1] += 1;
                            contact_order_array[index - 1] = 1;
                            flatknots[index3 - 1] = self.my_parameters[ii];
                            parameters[index - 1] = self.my_parameters[ii];
                            poles[index - 1] = self.my_tangents[ii];
                            index += 1;
                            index3 += 1;
                        }
                        mult_index += 1;
                        index1 += 1;
                        index2 += 1;
                    }
                } else {
                    index1 = 2;
                    for ii in 0..self.my_parameters.len() {
                        parameters[index1 - 1] = self.my_parameters[ii];
                        index1 += 1;
                    }
                    index = degree + 1;
                    for ii in 0..self.my_parameters.len() {
                        flatknots[index - 1] = self.my_parameters[ii];
                        index += 1;
                    }
                    index = 3;
                    for ii in 1..self.my_points.len() - 1 {
                        poles[index - 1] = self.my_points[ii];
                        index += 1;
                    }
                }
                poles[num_poles - 2] = self.my_tangents[num_points - 1];
                contact_order_array[num_poles - 2] = 1;
                parameters[num_poles - 1] = self.my_parameters[n - 1];
                parameters[num_poles - 2] = self.my_parameters[n - 1];
                poles[num_poles - 1] = self.my_points[num_points - 1];
                let inversion_problem = bspl_interpolate(
                    degree,
                    &flatknots,
                    &parameters,
                    &contact_order_array,
                    1,
                    &mut poles,
                );
                if inversion_problem == 0 {
                    self.my_curve = Some(Rc::new(RefCell::new(LawBSpline::new(
                        &poles,
                        &self.my_parameters,
                        &mults,
                        degree,
                        false,
                    ))));
                    self.my_is_done = true;
                }
            }
        }
    }

    /// OCCT Curve() (L562-567) — raises StdFail_NotDone when not done.
    pub fn curve(&self) -> Rc<RefCell<LawBSpline>> {
        assert!(self.my_is_done, "StdFail_NotDone");
        self.my_curve.as_ref().unwrap().clone()
    }

    /// OCCT IsDone() (L569-572).
    pub fn is_done(&self) -> bool {
        self.my_is_done
    }
}
