use crate::tolerance::*;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3, Curve2dEval, CurveEval};
use std::collections::BTreeMap;

/// Move a 2D B-spline so that its point at `u` reaches `target`.
///
/// The edit distributes the displacement across the non-zero rational basis
/// functions at `u` using the minimum-norm update that satisfies the point
/// constraint for fixed knots and weights.
pub fn move_bspline2_point(curve: &BSplineCurve2, u: f64, target: DVec2) -> BSplineCurve2 {
    let current = curve.point_at(u);
    let delta = target - current;
    let coeffs = rational_basis_coefficients(
        curve.degree,
        &curve.knots,
        &curve.weights,
        curve.control_points.len(),
        u,
    );
    let denom: f64 = coeffs.iter().map(|(_, c)| c * c).sum();
    if denom <= TOLERANCE_VEC_SQ_MIN {
        return curve.clone();
    }

    let mut edited = curve.clone();
    for (idx, coeff) in coeffs {
        edited.control_points[idx] += delta * (coeff / denom);
    }
    edited
}

/// Move a 3D B-spline so that its point at `u` reaches `target`.
///
/// See [`move_bspline2_point`] for the edit model.
pub fn move_bspline3_point(curve: &BSplineCurve3, u: f64, target: DVec3) -> BSplineCurve3 {
    let current = curve.point_at(u);
    let delta = target - current;
    let coeffs = rational_basis_coefficients(
        curve.degree,
        &curve.knots,
        &curve.weights,
        curve.control_points.len(),
        u,
    );
    let denom: f64 = coeffs.iter().map(|(_, c)| c * c).sum();
    if denom <= TOLERANCE_VEC_SQ_MIN {
        return curve.clone();
    }

    let mut edited = curve.clone();
    for (idx, coeff) in coeffs {
        edited.control_points[idx] += delta * (coeff / denom);
    }
    edited
}

/// Move a 2D B-spline so that its point and first derivative at `u` match
/// the supplied targets as closely as possible for fixed knots and weights.
pub fn move_bspline2_tangent(
    curve: &BSplineCurve2,
    u: f64,
    target_point: DVec2,
    target_derivative: DVec2,
) -> BSplineCurve2 {
    let point_delta = target_point - curve.point_at(u);
    let derivative_delta = target_derivative - curve.derivative_at(u);
    if let Some(edited) =
        move_bspline2_tangent_occt_style(curve, u, point_delta, derivative_delta, 1, 1)
    {
        return edited;
    }

    let coeffs = rational_basis_and_derivative_coefficients(
        curve.degree,
        &curve.knots,
        &curve.weights,
        curve.control_points.len(),
        u,
    );
    let Some(update) = solve_two_constraint_update(&coeffs, point_delta, derivative_delta) else {
        return curve.clone();
    };

    let mut edited = curve.clone();
    for (idx, delta) in update {
        edited.control_points[idx] += delta;
    }
    edited
}

/// 3D analogue of [`move_bspline2_tangent`].
pub fn move_bspline3_tangent(
    curve: &BSplineCurve3,
    u: f64,
    target_point: DVec3,
    target_derivative: DVec3,
) -> BSplineCurve3 {
    let point_delta = target_point - curve.point_at(u);
    let derivative_delta = target_derivative - curve.derivative_at(u);
    if let Some(edited) =
        move_bspline3_tangent_occt_style(curve, u, point_delta, derivative_delta, 1, 1)
    {
        return edited;
    }

    let coeffs = rational_basis_and_derivative_coefficients(
        curve.degree,
        &curve.knots,
        &curve.weights,
        curve.control_points.len(),
        u,
    );
    let Some(update) = solve_two_constraint_update(&coeffs, point_delta, derivative_delta) else {
        return curve.clone();
    };

    let mut edited = curve.clone();
    for (idx, delta) in update {
        edited.control_points[idx] += delta;
    }
    edited
}

fn move_bspline2_tangent_occt_style(
    curve: &BSplineCurve2,
    u: f64,
    point_delta: DVec2,
    derivative_delta: DVec2,
    start_condition: isize,
    end_condition: isize,
) -> Option<BSplineCurve2> {
    let (first_fn, second_fn) = occt_auxiliary_functions(
        curve.degree,
        &curve.knots,
        curve.control_points.len(),
        u,
        start_condition,
        end_condition,
    )?;
    let (a, b) = occt_auxiliary_solution(
        curve.degree,
        &curve.knots,
        &curve.weights,
        u,
        &first_fn,
        &second_fn,
        point_delta,
        derivative_delta,
    )?;

    let mut edited = curve.clone();
    for i in 0..edited.control_points.len() {
        edited.control_points[i] += first_fn[i] * a + second_fn[i] * b;
    }
    Some(edited)
}

fn move_bspline3_tangent_occt_style(
    curve: &BSplineCurve3,
    u: f64,
    point_delta: DVec3,
    derivative_delta: DVec3,
    start_condition: isize,
    end_condition: isize,
) -> Option<BSplineCurve3> {
    let (first_fn, second_fn) = occt_auxiliary_functions(
        curve.degree,
        &curve.knots,
        curve.control_points.len(),
        u,
        start_condition,
        end_condition,
    )?;
    let (a, b) = occt_auxiliary_solution(
        curve.degree,
        &curve.knots,
        &curve.weights,
        u,
        &first_fn,
        &second_fn,
        point_delta,
        derivative_delta,
    )?;

    let mut edited = curve.clone();
    for i in 0..edited.control_points.len() {
        edited.control_points[i] += first_fn[i] * a + second_fn[i] * b;
    }
    Some(edited)
}

fn occt_auxiliary_functions(
    degree: usize,
    knots: &[f64],
    control_point_count: usize,
    u: f64,
    start_condition: isize,
    end_condition: isize,
) -> Option<(Vec<f64>, Vec<f64>)> {
    if start_condition < -1 || end_condition < -1 {
        return None;
    }
    let conditions = start_condition + end_condition + 4;
    if conditions as usize > control_point_count || degree == 0 {
        return None;
    }

    let greville = greville_points(degree, knots, control_point_count)?;
    let start_index = (start_condition + 1).max(0) as usize;
    let end_index = control_point_count.checked_sub((end_condition + 2).max(0) as usize)?;
    if start_index > end_index || end_index >= control_point_count {
        return None;
    }

    let index = locate_greville(&greville, u, start_index, end_index);
    let other_index = if index == start_index {
        index + 1
    } else if index == end_index {
        index.saturating_sub(1)
    } else if index + 1 < knots.len() && u - knots[index] < knots[index + 1] - u {
        index.saturating_sub(1)
    } else {
        index + 1
    };

    let start_value = if start_index == 0 {
        greville[0] - (greville[control_point_count - 1] - greville[0])
    } else {
        greville[start_index - 1]
    };
    let end_value = if end_index + 1 >= control_point_count {
        greville[control_point_count - 1] + (greville[control_point_count - 1] - greville[0])
    } else {
        greville[end_index + 1]
    };

    let first = auxiliary_function_values(
        &greville,
        start_index,
        end_index,
        index,
        start_value,
        end_value,
    )?;
    let second = auxiliary_function_values(
        &greville,
        start_index,
        end_index,
        other_index,
        start_value,
        end_value,
    )?;
    Some((first, second))
}

fn greville_points(degree: usize, knots: &[f64], control_point_count: usize) -> Option<Vec<f64>> {
    if degree == 0 || knots.len() < control_point_count + degree + 1 {
        return None;
    }
    let mut points = Vec::with_capacity(control_point_count);
    for i in 0..control_point_count {
        let sum: f64 = (1..=degree).map(|j| knots[i + j]).sum();
        points.push(sum / degree as f64);
    }
    Some(points)
}

fn locate_greville(greville: &[f64], u: f64, start: usize, end: usize) -> usize {
    let mut best = start;
    let mut best_dist = (greville[start] - u).abs();
    for (i, value) in greville.iter().enumerate().take(end + 1).skip(start) {
        let dist = (*value - u).abs();
        if dist < best_dist {
            best = i;
            best_dist = dist;
        }
    }
    best
}

fn auxiliary_function_values(
    greville: &[f64],
    start: usize,
    end: usize,
    pivot: usize,
    start_value: f64,
    end_value: f64,
) -> Option<Vec<f64>> {
    let mut values = vec![0.0; greville.len()];
    let left_den = greville[pivot] - start_value;
    let right_den = end_value - greville[pivot];
    if left_den.abs() <= TOLERANCE_FLOAT_LOOSE || right_den.abs() <= TOLERANCE_FLOAT_LOOSE {
        return None;
    }

    for i in start..=pivot {
        values[i] = ((greville[i] - start_value) / left_den).powi(3);
    }
    for i in pivot..=end {
        values[i] = ((end_value - greville[i]) / right_den).powi(3);
    }
    Some(values)
}

fn occt_auxiliary_solution<V>(
    degree: usize,
    knots: &[f64],
    weights: &[f64],
    u: f64,
    first_fn: &[f64],
    second_fn: &[f64],
    point_delta: V,
    derivative_delta: V,
) -> Option<(V, V)>
where
    V: Copy + std::ops::Add<Output = V> + std::ops::Mul<f64, Output = V>,
{
    let first = scalar_rational_value_and_derivative(degree, knots, weights, first_fn, u)?;
    let second = scalar_rational_value_and_derivative(degree, knots, weights, second_fn, u)?;
    let det = first.0 * second.1 - second.0 * first.1;
    if det.abs() <= TOLERANCE_FLOAT_LOOSE {
        return None;
    }

    let a = point_delta * (second.1 / det) + derivative_delta * (-second.0 / det);
    let b = point_delta * (-first.1 / det) + derivative_delta * (first.0 / det);
    Some((a, b))
}

fn scalar_rational_value_and_derivative(
    degree: usize,
    knots: &[f64],
    weights: &[f64],
    values: &[f64],
    u: f64,
) -> Option<(f64, f64)> {
    let coeffs =
        rational_basis_and_derivative_coefficients(degree, knots, weights, values.len(), u);
    if coeffs.is_empty() {
        return None;
    }
    let value = coeffs.iter().map(|(idx, c, _)| c * values[*idx]).sum();
    let derivative = coeffs.iter().map(|(idx, _, d)| d * values[*idx]).sum();
    Some((value, derivative))
}

fn rational_basis_coefficients(
    degree: usize,
    knots: &[f64],
    weights: &[f64],
    control_point_count: usize,
    u: f64,
) -> Vec<(usize, f64)> {
    if control_point_count == 0
        || weights.len() != control_point_count
        || knots.len() < control_point_count + degree + 1
    {
        return Vec::new();
    }

    let span = find_span(control_point_count, degree, knots, u);
    let basis = basis_functions(span, u, degree, knots);
    let first = span.saturating_sub(degree);

    let weighted: Vec<(usize, f64)> = basis
        .iter()
        .enumerate()
        .filter_map(|(local, basis_value)| {
            let idx = first + local;
            (idx < control_point_count).then_some((idx, basis_value * weights[idx]))
        })
        .collect();
    let weight_sum: f64 = weighted.iter().map(|(_, value)| *value).sum();
    if weight_sum.abs() <= TOLERANCE_VEC_SQ_MIN {
        return Vec::new();
    }

    weighted
        .into_iter()
        .map(|(idx, weighted_value)| (idx, weighted_value / weight_sum))
        .collect()
}

fn rational_basis_and_derivative_coefficients(
    degree: usize,
    knots: &[f64],
    weights: &[f64],
    control_point_count: usize,
    u: f64,
) -> Vec<(usize, f64, f64)> {
    let coeffs = rational_basis_coefficients(degree, knots, weights, control_point_count, u);
    if coeffs.is_empty() {
        return Vec::new();
    }
    let eps = ((knots[knots.len() - degree - 1] - knots[degree]).abs() * TOLERANCE_MESH_LEGACY).max(TOLERANCE_LINEAR_RELAX_8);
    let u0 = (u - eps).max(knots[degree]);
    let u1 = (u + eps).min(knots[knots.len() - degree - 1]);
    if (u1 - u0).abs() <= TOLERANCE_FLOAT_LOOSE {
        return coeffs.into_iter().map(|(idx, c)| (idx, c, 0.0)).collect();
    }

    let left = rational_basis_coefficients(degree, knots, weights, control_point_count, u0);
    let right = rational_basis_coefficients(degree, knots, weights, control_point_count, u1);
    let mut deriv = BTreeMap::new();
    for (idx, value) in left {
        *deriv.entry(idx).or_insert(0.0) -= value / (u1 - u0);
    }
    for (idx, value) in right {
        *deriv.entry(idx).or_insert(0.0) += value / (u1 - u0);
    }

    coeffs
        .into_iter()
        .map(|(idx, c)| (idx, c, deriv.get(&idx).copied().unwrap_or(0.0)))
        .collect()
}

fn solve_two_constraint_update<V>(
    coeffs: &[(usize, f64, f64)],
    point_delta: V,
    derivative_delta: V,
) -> Option<Vec<(usize, V)>>
where
    V: Copy + std::ops::Add<Output = V> + std::ops::Mul<f64, Output = V>,
{
    let aa: f64 = coeffs.iter().map(|(_, c, _)| c * c).sum();
    let ab: f64 = coeffs.iter().map(|(_, c, d)| c * d).sum();
    let bb: f64 = coeffs.iter().map(|(_, _, d)| d * d).sum();
    let det = aa * bb - ab * ab;
    if det.abs() <= TOLERANCE_VEC_SQ_MIN {
        return None;
    }

    let alpha = point_delta * (bb / det) + derivative_delta * (-ab / det);
    let beta = point_delta * (-ab / det) + derivative_delta * (aa / det);
    Some(
        coeffs
            .iter()
            .map(|(idx, c, d)| (*idx, alpha * *c + beta * *d))
            .collect(),
    )
}

fn find_span(control_point_count: usize, degree: usize, knots: &[f64], u: f64) -> usize {
    let n = control_point_count - 1;
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[degree] {
        return degree;
    }

    let mut low = degree;
    let mut high = n + 1;
    let mut mid = (low + high) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

fn basis_functions(span: usize, u: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
    let mut basis = vec![0.0; degree + 1];
    let mut left = vec![0.0; degree + 1];
    let mut right = vec![0.0; degree + 1];
    basis[0] = 1.0;

    for j in 1..=degree {
        left[j] = u - knots[span + 1 - j];
        right[j] = knots[span + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() <= TOLERANCE_VEC_SQ_MIN {
                0.0
            } else {
                basis[r] / denom
            };
            basis[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        basis[j] = saved;
    }
    basis
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn move_bspline2_point_hits_target_for_linear_curve() {
        let curve = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(2.0, 0.0)],
            weights: vec![1.0, 1.0],
        };

        let edited = move_bspline2_point(&curve, 0.5, DVec2::new(1.0, 3.0));

        assert!((edited.point_at(0.5) - DVec2::new(1.0, 3.0)).length() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn move_bspline3_point_hits_target_for_linear_curve() {
        let curve = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };

        let edited = move_bspline3_point(&curve, 0.5, DVec3::new(1.0, 3.0, -2.0));

        assert!((edited.point_at(0.5) - DVec3::new(1.0, 3.0, -2.0)).length() < TOLERANCE_LEN_MIN);
    }

    #[test]
    fn move_bspline2_tangent_hits_point_and_derivative_for_linear_curve() {
        let curve = BSplineCurve2 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec2::new(0.0, 0.0), DVec2::new(2.0, 0.0)],
            weights: vec![1.0, 1.0],
        };

        let edited = move_bspline2_tangent(&curve, 0.5, DVec2::new(1.0, 3.0), DVec2::new(4.0, 2.0));

        assert!((edited.point_at(0.5) - DVec2::new(1.0, 3.0)).length() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((edited.derivative_at(0.5) - DVec2::new(4.0, 2.0)).length() < TOLERANCE_LINEAR_RELAX_8);
    }

    #[test]
    fn move_bspline3_tangent_hits_point_and_derivative_for_linear_curve() {
        let curve = BSplineCurve3 {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };

        let edited = move_bspline3_tangent(
            &curve,
            0.5,
            DVec3::new(1.0, 3.0, -2.0),
            DVec3::new(4.0, 2.0, 1.0),
        );

        assert!((edited.point_at(0.5) - DVec3::new(1.0, 3.0, -2.0)).length() < TOLERANCE_LINEAR_ULTRA_STRICT);
        assert!((edited.derivative_at(0.5) - DVec3::new(4.0, 2.0, 1.0)).length() < TOLERANCE_LINEAR_RELAX_8);
    }

    #[test]
    fn move_bspline2_tangent_matches_occt_draw_b1_length() {
        let mut curve = BSplineCurve2 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0],
            control_points: vec![
                DVec2::new(0.0, 0.0),
                DVec2::new(1.0, 0.5),
                DVec2::new(2.0, 1.0),
                DVec2::new(2.0, 2.0),
                DVec2::new(3.0, 1.5),
                DVec2::new(4.0, 1.5),
                DVec2::new(5.0, 2.0),
            ],
            weights: vec![1.0; 7],
        };
        let u = 2.0;
        let point = curve.point_at(u);
        let tangent = curve.derivative_at(u);
        let mut dyvalue = tangent.y;
        for _ in 0..100 {
            curve = move_bspline2_tangent(&curve, u, point, DVec2::new(tangent.x, dyvalue));
            dyvalue += 0.005;
        }

        let length = sample_bspline2_length(&curve, 4096);

        assert!(
            (length - 5.9590472422107315).abs() < TOLERANCE_RETRY_LADDER_COARSE,
            "length={length}"
        );
    }

    fn sample_bspline2_length(curve: &BSplineCurve2, samples: usize) -> f64 {
        let u0 = curve.knots[curve.degree];
        let u1 = curve.knots[curve.knots.len() - curve.degree - 1];
        let mut length = 0.0;
        let mut prev = curve.point_at(u0);
        for i in 1..=samples {
            let u = u0 + (u1 - u0) * (i as f64 / samples as f64);
            let current = curve.point_at(u);
            length += (current - prev).length();
            prev = current;
        }
        length
    }
}
