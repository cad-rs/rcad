//! OCCT Plate_D1 / Plate_D2 / Plate_D3 (TKGeomAlgo/Plate) — 1:1 port.
//!
//! Carriers for the first, second and third derivatives of a surface at a
//! point.  OCCT grants friend access to Plate_GtoCConstraint and
//! Plate_FreeGtoCConstraint; in Rust the same members are reachable through
//! the pub `du()`-style accessors below.
//! gp_XYZ -> DVec3 (architecture mapping).

use glam::DVec3;

/// OCCT Plate_D1 (D1.hxx L28-47, D1.cxx L192-200).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateD1 {
    du: DVec3,
    dv: DVec3,
}

impl PlateD1 {
    /// OCCT Plate_D1(du, dv) (D1.cxx L192-196).
    pub fn new(du: DVec3, dv: DVec3) -> Self {
        PlateD1 { du, dv }
    }

    /// OCCT DU() (D1.lxx L68-71) — also the OCCT friend member `Du`.
    pub fn du(&self) -> DVec3 {
        self.du
    }

    /// OCCT DV() (D1.lxx L73-76) — also the OCCT friend member `Dv`.
    pub fn dv(&self) -> DVec3 {
        self.dv
    }
}

/// OCCT Plate_D2 (D2.hxx L104-120, D2.cxx L219-228).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateD2 {
    duu: DVec3,
    duv: DVec3,
    dvv: DVec3,
}

impl PlateD2 {
    /// OCCT Plate_D2(duu, duv, dvv) (D2.cxx L219-224).
    pub fn new(duu: DVec3, duv: DVec3, dvv: DVec3) -> Self {
        PlateD2 { duu, duv, dvv }
    }

    /// OCCT DUU() — also the OCCT friend member `Duu`.
    pub fn duu(&self) -> DVec3 {
        self.duu
    }

    /// OCCT DUV() — also the OCCT friend member `Duv`.
    pub fn duv(&self) -> DVec3 {
        self.duv
    }

    /// OCCT DVV() — also the OCCT friend member `Dvv`.
    pub fn dvv(&self) -> DVec3 {
        self.dvv
    }
}

/// OCCT Plate_D3 (D3.hxx L150-170, D3.cxx L247-257).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateD3 {
    duuu: DVec3,
    duuv: DVec3,
    duvv: DVec3,
    dvvv: DVec3,
}

impl PlateD3 {
    /// OCCT Plate_D3(duuu, duuv, duvv, dvvv) (D3.cxx L247-253).
    pub fn new(duuu: DVec3, duuv: DVec3, duvv: DVec3, dvvv: DVec3) -> Self {
        PlateD3 {
            duuu,
            duuv,
            duvv,
            dvvv,
        }
    }

    /// OCCT DUUU() — also the OCCT friend member `Duuu`.
    pub fn duuu(&self) -> DVec3 {
        self.duuu
    }

    /// OCCT DUUV() — also the OCCT friend member `Duuv`.
    pub fn duuv(&self) -> DVec3 {
        self.duuv
    }

    /// OCCT DUVV() — also the OCCT friend member `Duvv`.
    pub fn duvv(&self) -> DVec3 {
        self.duvv
    }

    /// OCCT DVVV() — also the OCCT friend member `Dvvv`.
    pub fn dvvv(&self) -> DVec3 {
        self.dvvv
    }
}
