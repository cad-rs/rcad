//! OCCT Geom2dAPI_Interpolate (Geom2dAPI_Interpolate.hxx/.cxx) — constrained
//! B-spline interpolation through a table of 2D points, with optional
//! tangency constraints (C1 at tangent points, C2 elsewhere) and periodic
//! closure.
//!
//! 1:1 translation of Geom2dAPI_Interpolate.cxx (whole file, L31-838) using
//! the rcad OCCT-aligned ports of its dependencies:
//!   - `crate::math::plib::eval_lagrange`  (PLib::EvalLagrange, PLib.cxx L1122-1249)
//!   - `crate::math::bspl_lib::interpolate` (BSplCLib::Interpolate, BSplCLib.cxx L3353-3395)
//!
//! The produced curve is a non-rational `BSplineCurve2` (all weights 1);
//! OCCT's periodic flag is not representable in `BSplineCurve2`, so a
//! periodic result is built with the periodic knot sequence and the flag is
//! dropped (architecture note).

use glam::DVec2;

use crate::geom::BSplineCurve2;
use crate::math::bspl_lib::interpolate as bspl_interpolate;
use crate::math::plib::eval_lagrange;

/// OCCT Geom2dAPI_Interpolate::CheckPoints (Geom2dAPI_Interpolate.cxx L33-44).
fn check_points(points: &[DVec2], tolerance: f64) -> bool {
    let tolerance_squared = tolerance * tolerance;
    let mut result = true;
    for ii in 0..points.len() - 1 {
        if !result {
            break;
        }
        let distance_squared = points[ii].distance_squared(points[ii + 1]);
        result = distance_squared >= tolerance_squared;
    }
    result
}

/// OCCT Geom2dAPI_Interpolate::CheckTangents (L48-66).
fn check_tangents(tangents: &[DVec2], tangent_flags: &[bool], tolerance: f64) -> bool {
    let tolerance_squared = tolerance * tolerance;
    let mut result = true;
    let mut index = 0;
    for ii in 0..tangents.len() {
        if !result {
            break;
        }
        if tangent_flags[index] {
            let distance_squared = tangents[ii].length_squared();
            result = distance_squared >= tolerance_squared;
        }
        index += 1;
    }
    result
}

/// OCCT Geom2dAPI_Interpolate::CheckParameters (L70-81).
fn check_parameters(parameters: &[f64]) -> bool {
    let mut result = true;
    for ii in 0..parameters.len() - 1 {
        if !result {
            break;
        }
        let distance = parameters[ii + 1] - parameters[ii];
        result = distance >= f64::MIN_POSITIVE;
    }
    result
}

/// OCCT Geom2dAPI_Interpolate::BuildParameters (L85-111) — cumulative chord
/// lengths; one extra closing parameter when periodic.
fn build_parameters(periodic: bool, points: &[DVec2]) -> Vec<f64> {
    let mut num_parameters = points.len();
    if periodic {
        num_parameters += 1;
    }
    let mut parameters = vec![0.0; num_parameters];
    parameters[0] = 0.0e0;
    let mut index = 1; // 0-based slot 2 (OCCT 1-based index 2)
    for ii in 0..points.len() - 1 {
        let distance = points[ii].distance(points[ii + 1]);
        parameters[index] = parameters[ii] + distance;
        index += 1;
    }
    if periodic {
        let distance = points[points.len() - 1].distance(points[0]);
        parameters[index] = parameters[points.len() - 1] + distance;
    }
    parameters
}

/// OCCT Geom2dAPI_Interpolate::BuildPeriodicTangent (L115-153).
fn build_periodic_tangent(
    points: &[DVec2],
    tangents: &mut [DVec2],
    tangent_flags: &mut [bool],
    parameters: &[f64],
) {
    if points.len() < 3 {
        panic!("Standard_ConstructionError: too few points");
    }

    if !tangent_flags[0] {
        let degree = if points.len() == 3 { 2 } else { 3 };
        let mut point_array = Vec::with_capacity((degree + 1) * 2);
        for p in &points[0..=degree] {
            point_array.push(p.x);
            point_array.push(p.y);
        }
        let mut eval_result = [0.0f64; 4];
        tangent_flags[0] = true;
        eval_lagrange(
            parameters[0],
            1,
            degree,
            2,
            &point_array,
            &parameters[0..=degree],
            &mut eval_result,
        );
        tangents[0] = DVec2::new(eval_result[2], eval_result[3]);
    }
}

/// OCCT Geom2dAPI_Interpolate::BuildTangents (L157-214) — estimate the
/// missing first/last tangents from a local Lagrange interpolation.
fn build_tangents(
    points: &[DVec2],
    tangents: &mut [DVec2],
    tangent_flags: &mut [bool],
    parameters: &[f64],
) {
    let mut degree = 3;
    if points.len() < 3 {
        panic!("Standard_ConstructionError: too few points");
    }
    if points.len() == 3 {
        degree = 2;
    }
    if !tangent_flags[0] {
        let mut point_array = Vec::with_capacity((degree + 1) * 2);
        for p in &points[0..=degree] {
            point_array.push(p.x);
            point_array.push(p.y);
        }
        tangent_flags[0] = true;
        let mut eval_result = [0.0f64; 4];
        eval_lagrange(
            parameters[0],
            1,
            degree,
            2,
            &point_array,
            &parameters[0..=degree],
            &mut eval_result,
        );
        tangents[0] = DVec2::new(eval_result[2], eval_result[3]);
    }
    if !tangent_flags[tangent_flags.len() - 1] {
        let start = points.len() - 1 - degree;
        let mut point_array = Vec::with_capacity((degree + 1) * 2);
        for p in &points[start..] {
            point_array.push(p.x);
            point_array.push(p.y);
        }
        tangent_flags[tangent_flags.len() - 1] = true;
        let mut eval_result = [0.0f64; 4];
        eval_lagrange(
            parameters[parameters.len() - 1],
            1,
            degree,
            2,
            &point_array,
            &parameters[start..],
            &mut eval_result,
        );
        tangents[tangents.len() - 1] = DVec2::new(eval_result[2], eval_result[3]);
    }
}

/// OCCT Geom2dAPI_Interpolate::ScaleTangents (L222-279) — scale the given
/// tangents so that they have the length of the derivative of the local
/// Lagrange interpolation.
fn scale_tangents(
    points: &[DVec2],
    tangents: &mut [DVec2],
    tangent_flags: &[bool],
    parameters: &[f64],
) {
    let num_points = points.len();
    let degree = if num_points == 2 {
        1
    } else if num_points >= 3 {
        2
    } else {
        0
    };

    let mut index = 0;
    for ii in 0..tangent_flags.len() {
        if tangent_flags[ii] {
            let mut point_array = Vec::with_capacity((degree + 1) * 2);
            for p in &points[index..=index + degree] {
                point_array.push(p.x);
                point_array.push(p.y);
            }
            let mut eval_result = [0.0f64; 4];
            eval_lagrange(
                parameters[ii],
                1,
                degree,
                2,
                &point_array,
                &parameters[index..=index + degree],
                &mut eval_result,
            );
            let mut value = [0.0f64; 2];
            for jj in 0..2 {
                value[0] += tangents[ii][jj].abs();
                value[1] += eval_result[2 + jj].abs();
            }
            let ratio = value[1] / value[0];
            tangents[ii] = DVec2::new(ratio * tangents[ii].x, ratio * tangents[ii].y);
            if ii != 0 {
                index += 1;
            }
            if index > num_points - 1 - degree {
                index = num_points - 1 - degree;
            }
        }
    }
}

/// Build a `BSplineCurve2` from poles, distinct knots and multiplicities
/// (OCCT `new Geom2d_BSplineCurve(Poles, Knots, Mults, Degree[, Periodic])`).
/// The full (expanded) knot vector is stored.
fn bspline_from_poles_knots_mults(
    poles: &[DVec2],
    knots: &[f64],
    mults: &[usize],
    degree: usize,
) -> BSplineCurve2 {
    let mut expanded = Vec::with_capacity(knots.len() + degree + 1);
    for (i, &k) in knots.iter().enumerate() {
        for _ in 0..mults[i] {
            expanded.push(k);
        }
    }
    BSplineCurve2 {
        degree,
        knots: expanded,
        control_points: poles.to_vec(),
        weights: vec![1.0; poles.len()],
    }
}

/// OCCT Geom2dAPI_Interpolate — constrained B-spline interpolation through
/// 2D points with optional tangency constraints.
pub struct Geom2dInterpolate {
    tolerance: f64,
    points: Vec<DVec2>,
    is_done: bool,
    curve: Option<BSplineCurve2>,
    parameters: Vec<f64>,
    periodic: bool,
    tangent_request: bool,
    tangents: Vec<DVec2>,
    tangent_flags: Vec<bool>,
}

impl Geom2dInterpolate {
    /// OCCT Geom2dAPI_Interpolate(Points, PeriodicFlag, Tolerance)
    /// (L283-309).  Panics (Standard_ConstructionError) when the points are
    /// too close together.
    pub fn new(points: Vec<DVec2>, periodic: bool, tolerance: f64) -> Self {
        let result = check_points(&points, tolerance);
        let mut interp = Geom2dInterpolate {
            tolerance,
            points,
            is_done: false,
            curve: None,
            parameters: Vec::new(),
            periodic,
            tangent_request: false,
            tangents: vec![DVec2::ZERO; 0],
            tangent_flags: vec![false; 0],
        };
        if !result {
            panic!("Standard_ConstructionError: points too close");
        }
        let n = interp.points.len();
        interp.tangents = vec![DVec2::ZERO; n];
        interp.tangent_flags = vec![false; n];
        interp.parameters = build_parameters(periodic, &interp.points);
        interp
    }

    /// OCCT Geom2dAPI_Interpolate(Points, Parameters, PeriodicFlag, Tolerance)
    /// (L313-354) — explicit parameter values.
    pub fn new_with_params(points: Vec<DVec2>, parameters: Vec<f64>, periodic: bool, tolerance: f64) -> Self {
        let result = check_points(&points, tolerance);
        if periodic && points.len() + 1 != parameters.len() {
            panic!("Standard_ConstructionError: parameters/points mismatch");
        }
        let n = points.len();
        let mut interp = Geom2dInterpolate {
            tolerance,
            points,
            is_done: false,
            curve: None,
            parameters,
            periodic,
            tangent_request: false,
            tangents: vec![DVec2::ZERO; n],
            tangent_flags: vec![false; n],
        };
        if !result {
            panic!("Standard_ConstructionError: points too close");
        }
        let result = check_parameters(&interp.parameters);
        if !result {
            panic!("Standard_ConstructionError: bad parameters");
        }
        interp
    }

    /// OCCT Geom2dAPI_Interpolate::Load(Tangents, TangentFlags, Scale)
    /// (L358-391).
    pub fn load(&mut self, tangents: &[DVec2], tangent_flags: &[bool], scale: bool) {
        self.tangent_request = true;
        self.tangent_flags = tangent_flags.to_vec();
        if tangents.len() != self.points.len() || tangent_flags.len() != self.points.len() {
            panic!("Standard_ConstructionError: length mismatch");
        }
        let result = check_tangents(tangents, tangent_flags, self.tolerance);
        if result {
            self.tangents = tangents.to_vec();
            if scale {
                scale_tangents(&self.points, &mut self.tangents, tangent_flags, &self.parameters);
            }
        } else {
            panic!("Standard_ConstructionError: bad tangents");
        }
    }

    /// OCCT Geom2dAPI_Interpolate::Load(InitialTangent, FinalTangent, Scale)
    /// (L395-417).
    pub fn load_initial_final(&mut self, initial_tangent: DVec2, final_tangent: DVec2, scale: bool) {
        self.tangent_request = true;
        self.tangent_flags[0] = true;
        self.tangent_flags[self.points.len() - 1] = true;
        self.tangents[0] = initial_tangent;
        self.tangents[self.points.len() - 1] = final_tangent;
        let result = check_tangents(&self.tangents, &self.tangent_flags, self.tolerance);
        if !result {
            panic!("Standard_ConstructionError: bad tangents");
        }
        if scale {
            scale_tangents(&self.points, &mut self.tangents, &self.tangent_flags, &self.parameters);
        }
    }

    /// OCCT Geom2dAPI_Interpolate::Perform (L421-431).
    pub fn perform(&mut self) {
        if self.periodic {
            self.perform_periodic();
        } else {
            self.perform_non_periodic();
        }
    }

    /// OCCT Geom2dAPI_Interpolate::IsDone (L835-838).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// OCCT Geom2dAPI_Interpolate::Curve (L817-824) — panics (StdFail_NotDone)
    /// when the interpolation failed.
    pub fn curve(&self) -> BSplineCurve2 {
        if !self.is_done {
            panic!("StdFail_NotDone: interpolation not done");
        }
        self.curve.clone().expect("curve set when done")
    }

    /// OCCT Geom2dAPI_Interpolate::PerformPeriodic (L435-613).
    fn perform_periodic(&mut self) {
        let num_points = self.points.len();
        let period = self.parameters[self.parameters.len() - 1] - self.parameters[0];
        let mut num_poles = num_points + 1;
        if num_points == 2 && !self.tangent_request {
            // Build a periodic curve of degree 1.
            let degree = 1usize;
            let mults = vec![1usize; num_poles];
            let poles: Vec<DVec2> = self.points.clone();
            self.curve = Some(bspline_from_poles_knots_mults(&poles, &self.parameters, &mults, degree));
            self.is_done = true;
        } else {
            let num_distinct_knots = num_points + 1;
            let half_order = 2usize;
            let degree = 3usize;
            num_poles += 2;
            if self.tangent_request {
                for ii in 1..self.tangent_flags.len() {
                    if self.tangent_flags[ii] {
                        num_poles += 1;
                    }
                }
            }

            let n = num_poles;
            let mut parameters = vec![0.0; n];
            let mut flatknots = vec![0.0; n + degree + 1];
            let mut mults = vec![0usize; num_distinct_knots];
            let mut contact_order_array = vec![0i32; n];
            let mut poles = vec![DVec2::ZERO; n];

            for ii in 0..half_order {
                flatknots[ii] = self.parameters[self.parameters.len() - 2] - period;
                flatknots[ii + half_order] = self.parameters[0];
                flatknots[num_poles + ii] = self.parameters[self.parameters.len() - 1];
                // OCCT myParameters->Value(half_order) — 1-based index 2.
                flatknots[num_poles + half_order + ii] = self.parameters[half_order - 1] + period;
            }
            for ii in 1..num_distinct_knots - 1 {
                mults[ii] = 1;
            }
            mults[0] = half_order;
            mults[num_distinct_knots - 1] = half_order;
            if num_points >= 3 {
                build_periodic_tangent(&self.points, &mut self.tangents, &mut self.tangent_flags, &self.parameters);
            }
            contact_order_array[1] = 1;
            parameters[0] = self.parameters[0];
            parameters[1] = self.parameters[0];
            poles[0] = self.points[0];
            poles[1] = self.tangents[0];

            let mut mult_index = 1; // OCCT 1-based 2 -> 0-based 1
            let mut index = 2; // OCCT 3 -> 0-based 2
            let mut index1 = degree + 1; // OCCT degree+2 -> 0-based degree+1
            if self.tangent_request {
                for ii in 1..self.tangent_flags.len() {
                    parameters[index] = self.parameters[ii];
                    flatknots[index1] = self.parameters[ii];
                    poles[index] = self.points[ii];
                    index += 1;
                    index1 += 1;
                    if self.tangent_flags[ii] {
                        mults[mult_index] += 1;
                        contact_order_array[index] = 1;
                        parameters[index] = self.parameters[ii];
                        flatknots[index1] = self.parameters[ii];
                        poles[index] = self.tangents[ii];
                        index += 1;
                        index1 += 1;
                    }
                    mult_index += 1;
                }
            } else {
                let mut idx = degree; // OCCT degree+1 -> 0-based degree
                let mut idx1 = 1; // OCCT 2 -> 0-based 1
                for ii in 0..self.parameters.len() {
                    parameters[idx1] = self.parameters[ii];
                    flatknots[idx] = self.parameters[ii];
                    idx += 1;
                    idx1 += 1;
                }
                let mut idx2 = 2; // OCCT 3 -> 0-based 2
                for ii in 1..self.points.len() {
                    poles[idx2] = self.points[ii];
                    idx2 += 1;
                }
            }
            contact_order_array[num_poles - 2] = 1;
            parameters[num_poles - 2] = self.parameters[self.parameters.len() - 1];
            poles[num_poles - 2] = self.tangents[0];
            parameters[num_poles - 1] = self.parameters[self.parameters.len() - 1];
            poles[num_poles - 1] = self.points[0];

            let mut poles_flat = vec![0.0; num_poles * 2];
            for (i, p) in poles.iter().enumerate() {
                poles_flat[2 * i] = p.x;
                poles_flat[2 * i + 1] = p.y;
            }
            let inversion_problem = bspl_interpolate(
                degree,
                &flatknots,
                &parameters,
                &contact_order_array,
                2,
                &mut poles_flat,
            );
            if inversion_problem == 0 {
                let new_poles: Vec<DVec2> = poles_flat
                    .chunks_exact(2)
                    .take(num_poles - 2)
                    .map(|c| DVec2::new(c[0], c[1]))
                    .collect();
                self.curve = Some(bspline_from_poles_knots_mults(
                    &new_poles,
                    &self.parameters,
                    &mults,
                    degree,
                ));
                self.is_done = true;
            }
        }
    }

    /// OCCT Geom2dAPI_Interpolate::PerformNonPeriodic (L617-810).
    fn perform_non_periodic(&mut self) {
        let num_points = self.points.len();
        let mut num_distinct_knots = num_points;
        let mut num_poles = num_points;
        let mut degree: usize;
        if num_poles == 2 && !self.tangent_request {
            degree = 1;
        } else if num_poles == 3 && !self.tangent_request {
            degree = 2;
            num_distinct_knots = 2;
        } else {
            degree = 3;
            num_poles += 2;
            if self.tangent_request {
                for ii in 1..self.tangent_flags.len() - 1 {
                    if self.tangent_flags[ii] {
                        num_poles += 1;
                    }
                }
            }
        }

        let n = num_poles;
        let mut parameters = vec![0.0; n];
        let mut flatknots = vec![0.0; n + degree + 1];
        let mut mults = vec![0usize; num_distinct_knots];
        let mut knots = vec![0.0; num_distinct_knots];
        let mut contact_order_array = vec![0i32; n];
        let mut poles = vec![DVec2::ZERO; n];

        for ii in 0..=degree {
            flatknots[ii] = self.parameters[0];
            flatknots[ii + num_poles] = self.parameters[num_points - 1];
        }
        for ii in 1..num_distinct_knots - 1 {
            mults[ii] = 1;
        }
        mults[0] = degree + 1;
        mults[num_distinct_knots - 1] = degree + 1;

        match degree {
            1 => {
                for ii in 0..num_poles {
                    poles[ii] = self.points[ii];
                }
                self.curve = Some(bspline_from_poles_knots_mults(&poles, &self.parameters, &mults, degree));
                self.is_done = true;
            }
            2 => {
                knots[0] = self.parameters[0];
                knots[1] = self.parameters[2];
                for ii in 0..num_poles {
                    poles[ii] = self.points[ii];
                }
                let mut poles_flat = vec![0.0; num_poles * 2];
                for (i, p) in poles.iter().enumerate() {
                    poles_flat[2 * i] = p.x;
                    poles_flat[2 * i + 1] = p.y;
                }
                let inversion_problem = bspl_interpolate(
                    degree,
                    &flatknots,
                    &self.parameters,
                    &contact_order_array,
                    2,
                    &mut poles_flat,
                );
                if inversion_problem == 0 {
                    let solved: Vec<DVec2> = poles_flat
                        .chunks_exact(2)
                        .map(|c| DVec2::new(c[0], c[1]))
                        .collect();
                    self.curve = Some(bspline_from_poles_knots_mults(&solved, &knots, &mults, degree));
                    self.is_done = true;
                }
            }
            _ => {
                // degree 3
                if num_points >= 3 {
                    build_tangents(&self.points, &mut self.tangents, &mut self.tangent_flags, &self.parameters);
                }
                contact_order_array[1] = 1;
                parameters[0] = self.parameters[0];
                parameters[1] = self.parameters[0];
                poles[0] = self.points[0];
                poles[1] = self.tangents[0];
                let mut mult_index = 1; // OCCT 2 -> 0-based 1
                let mut index = 2; // OCCT 3 -> 0-based 2
                let mut index1 = 1; // OCCT 2 -> 0-based 1
                let mut index2 = 1; // OCCT Lower+1 -> 0-based 1
                let mut index3 = degree + 1; // OCCT degree+2 -> 0-based degree+1
                if self.tangent_request {
                    for ii in 1..self.parameters.len() - 1 {
                        parameters[index] = self.parameters[ii];
                        poles[index] = self.points[index2];
                        flatknots[index3] = self.parameters[ii];
                        index += 1;
                        index3 += 1;
                        if self.tangent_flags[index1] {
                            mults[mult_index] += 1;
                            contact_order_array[index] = 1;
                            flatknots[index3] = self.parameters[ii];
                            parameters[index] = self.parameters[ii];
                            poles[index] = self.tangents[ii];
                            index += 1;
                            index3 += 1;
                        }
                        mult_index += 1;
                        index1 += 1;
                        index2 += 1;
                    }
                } else {
                    let mut idx1 = 1; // OCCT 2 -> 0-based 1
                    for ii in 0..self.parameters.len() {
                        parameters[idx1] = self.parameters[ii];
                        idx1 += 1;
                    }
                    let mut idx2 = 2; // OCCT 3 -> 0-based 2
                    for ii in 1..self.points.len() - 1 {
                        poles[idx2] = self.points[ii];
                        idx2 += 1;
                    }
                    let mut idx3 = degree; // OCCT degree+1 -> 0-based degree
                    for ii in 0..self.parameters.len() {
                        flatknots[idx3] = self.parameters[ii];
                        idx3 += 1;
                    }
                }
                poles[num_poles - 2] = self.tangents[num_points - 1];
                contact_order_array[num_poles - 2] = 1;
                parameters[num_poles - 1] = self.parameters[self.parameters.len() - 1];
                parameters[num_poles - 2] = self.parameters[self.parameters.len() - 1];
                poles[num_poles - 1] = self.points[num_points - 1];

                let mut poles_flat = vec![0.0; num_poles * 2];
                for (i, p) in poles.iter().enumerate() {
                    poles_flat[2 * i] = p.x;
                    poles_flat[2 * i + 1] = p.y;
                }
                let inversion_problem = bspl_interpolate(
                    degree,
                    &flatknots,
                    &parameters,
                    &contact_order_array,
                    2,
                    &mut poles_flat,
                );
                if inversion_problem == 0 {
                    let solved: Vec<DVec2> = poles_flat
                        .chunks_exact(2)
                        .map(|c| DVec2::new(c[0], c[1]))
                        .collect();
                    self.curve = Some(bspline_from_poles_knots_mults(
                        &solved,
                        &self.parameters,
                        &mults,
                        degree,
                    ));
                    self.is_done = true;
                }
            }
        }
    }
}
