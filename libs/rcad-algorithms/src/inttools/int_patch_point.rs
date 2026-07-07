//! OCCT IntPatch_Point.hxx / IntSurf_PntOn2S.hxx

/// OCCT IntPatch_Point.hxx — single intersection point on two surfaces.
#[derive(Debug, Clone)]
pub struct IntPatchPoint {
    pub p1: glam::DVec3,
    pub p2: glam::DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
    pub tolerance: f64,
}

/// OCCT IntSurf_PntOn2S.hxx — point on two surfaces (used as parameter).
#[derive(Debug, Clone)]
pub struct PntOn2S {
    pub p1: glam::DVec3,
    pub p2: glam::DVec3,
    pub u1: f64,
    pub v1: f64,
    pub u2: f64,
    pub v2: f64,
}
