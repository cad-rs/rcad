//! OCCT AdvApprox_EvaluatorFunction (TKG3d).

/// OCCT AdvApprox_EvaluatorFunction — interface for a class implementing a
/// function to be approximated by [`ApproxAFunction`](super::ApproxAFunction).
///
/// OCCT uses raw in/out pointers; the Rust form passes the dimension by the
/// `result` slice length and returns the error code (`*ErrorCode`).
pub trait EvaluatorFunction {
    /// OCCT `Evaluate(Dimension, StartEnd, Parameter, DerivativeRequest,
    /// Result, ErrorCode)`.  `result` has `dimension` slots; `result[i]`
    /// corresponds to `Result[i]`.
    fn evaluate(&mut self, start_end: &[f64; 2], parameter: f64, derivative_request: i32, result: &mut [f64]) -> i32;
}
