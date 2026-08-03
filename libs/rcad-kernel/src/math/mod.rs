// OCCT Math* packages
pub mod bnd;
pub mod bspl;
pub mod bvh;
pub mod direct_polynomial_roots;
pub mod el;
pub mod math_poly;
pub mod newton_function_root;
pub mod plib;
pub mod root;
pub mod poly;
pub mod opt;
pub mod lin;
pub mod integ;
pub mod sys;

// keep flat — OCCT LProp package (partial)
pub mod top_loc;
pub mod curvature;

// OCCT CSLib package (surface normal computation)
pub mod cs_lib;

// legacy flat modules (keep during migration)
pub mod arc_length;
pub mod fit;
pub mod math_utils;
pub mod projection;
pub mod properties;

#[cfg(test)]
pub mod math_gtests;
#[cfg(test)]
pub mod tkmath_gtests;
