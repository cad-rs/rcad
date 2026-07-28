// OCCT Math* packages
pub mod bnd;
pub mod bspl;
pub mod bvh;
pub mod el;
pub mod root;
pub mod poly;
pub mod opt;
pub mod lin;
pub mod integ;
pub mod sys;

// OCCT GCPnts package
pub mod gcpnts;

// OCCT GProp package
pub mod gprop;

// OCCT BRepExtrema package
pub mod brep_extrema;

// keep flat — OCCT LProp package (partial)
pub mod curvature;

// legacy flat modules (keep during migration)
pub mod arc_length;
pub mod distance;
pub mod extrema;
pub mod fit;
pub mod math_utils;
pub mod projection;
pub mod properties;

#[cfg(test)]
pub mod math_gtests;
#[cfg(test)]
pub mod tkmath_gtests;
