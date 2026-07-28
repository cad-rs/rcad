//! OCCT MathInteg: numerical integration.
//!
//! Corresponds to OCCT `math_Integration`.
//!
//! Functions: simpson_integrate, gaussian_quadrature

// =============================================================================
// math_Integration — numerical integration
// =============================================================================

/// Simpson's rule for numerical integration.
///
/// Integrates f from a to b using n subintervals (n must be even).
pub fn simpson_integrate(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let n = if !n.is_multiple_of(2) { n + 1 } else { n };
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        let coef = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += coef * f(x);
    }
    sum * h / 3.0
}

/// Nodes and weights for Gaussian quadrature of various orders.
fn gaussian_nodes_weights(n: usize) -> Vec<(f64, f64)> {
    match n {
        1 => vec![(0.0, 2.0)],
        2 => vec![
            (-0.707_106_781_186_547_5, 1.0),
            (0.707_106_781_186_547_5, 1.0),
        ],
        3 => vec![
            (0.0, 0.888_888_888_888_888_8),
            (-0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
            (0.774_596_669_241_483_4, 0.555_555_555_555_555_6),
        ],
        4 => vec![
            (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
            (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
            (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
            (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
        ],
        5 => vec![
            (0.0, 0.568_888_888_888_888_9),
            (-0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
            (0.538_469_310_105_683_1, 0.478_628_670_499_366_5),
            (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
            (0.906_179_845_938_664, 0.236_926_885_056_189_1),
        ],
        6 => vec![
            (-0.932_469_514_203_152, 0.171_324_492_379_170_4),
            (-0.661_209_386_466_264_5, 0.360_761_573_048_138_6),
            (-0.238_619_186_083_196_9, 0.467_913_934_572_691),
            (0.238_619_186_083_196_9, 0.467_913_934_572_691),
            (0.661_209_386_466_264_5, 0.360_761_573_048_138_6),
            (0.932_469_514_203_152, 0.171_324_492_379_170_4),
        ],
        _ => gaussian_nodes_weights(6),
    }
}

/// Gaussian quadrature for numerical integration.
pub fn gaussian_quadrature(f: fn(f64) -> f64, a: f64, b: f64, n_points: usize) -> f64 {
    let nodes_weights = gaussian_nodes_weights(n_points);
    let scale = (b - a) / 2.0;
    let shift = (a + b) / 2.0;
    let mut sum = 0.0;
    for (node, weight) in nodes_weights {
        sum += weight * f(shift + scale * node);
    }
    sum * scale
}
