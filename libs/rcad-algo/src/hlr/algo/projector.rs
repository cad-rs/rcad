//! OCCT HLRAlgo_Projector (TKHLR HLRAlgo package).
//!
//! 1:1 translation of `HLRAlgo_Projector.hxx` (L44-124) +
//! `HLRAlgo_Projector.cxx` (L39-396) + `HLRAlgo_Projector.lxx` (L23-77).

use glam::{DVec2, DVec3};

use rcad_kernel::math::gp::{Ax2, Ax3, Lin, Trsf};

/// OCCT Precision::Angular() (Precision.hxx — 1.e-12).
const PRECISION_ANGULAR: f64 = 1.0e-12;

/// OCCT HLRAlgo_Projector — transforms and projects points and planes for
/// hidden line removal.
#[derive(Debug, Clone)]
pub struct Projector {
    my_type: i32,
    my_persp: bool,
    my_focus: f64,
    my_scaled_trsf: Trsf,
    my_trsf: Trsf,
    my_inv_trsf: Trsf,
    my_d1: DVec2,
    my_d2: DVec2,
    my_d3: DVec2,
}

impl Default for Projector {
    fn default() -> Self {
        Self::new()
    }
}

impl Projector {
    /// OCCT HLRAlgo_Projector() — cxx L39-44.
    pub fn new() -> Self {
        let mut p = Projector {
            my_type: 0,
            my_persp: false,
            my_focus: 0.0,
            my_scaled_trsf: Trsf::identity(),
            my_trsf: Trsf::identity(),
            my_inv_trsf: Trsf::identity(),
            my_d1: DVec2::ZERO,
            my_d2: DVec2::ZERO,
            my_d3: DVec2::ZERO,
        };
        p.scaled(false);
        p
    }

    /// OCCT HLRAlgo_Projector(CS) — cxx L48-55 (axonometric projector).
    pub fn from_ax2(cs: &Ax2) -> Self {
        let mut p = Projector {
            my_type: 0,
            my_persp: false,
            my_focus: 0.0,
            my_scaled_trsf: Trsf::identity(),
            my_trsf: Trsf::identity(),
            my_inv_trsf: Trsf::identity(),
            my_d1: DVec2::ZERO,
            my_d2: DVec2::ZERO,
            my_d3: DVec2::ZERO,
        };
        // myScaledTrsf.SetTransformation(CS) — the gp_Ax2 converts to gp_Ax3.
        p.my_scaled_trsf = Trsf::set_transformation_ax3(&Ax3::from_ax2(cs));
        p.scaled(false);
        p.set_direction();
        p
    }

    /// OCCT HLRAlgo_Projector(CS, Focus) — cxx L59-66 (perspective
    /// projector).
    pub fn from_ax2_focus(cs: &Ax2, focus: f64) -> Self {
        let mut p = Projector {
            my_type: 0,
            my_persp: true,
            my_focus: focus,
            my_scaled_trsf: Trsf::identity(),
            my_trsf: Trsf::identity(),
            my_inv_trsf: Trsf::identity(),
            my_d1: DVec2::ZERO,
            my_d2: DVec2::ZERO,
            my_d3: DVec2::ZERO,
        };
        p.my_scaled_trsf = Trsf::set_transformation_ax3(&Ax3::from_ax2(cs));
        p.scaled(false);
        p.set_direction();
        p
    }

    /// OCCT HLRAlgo_Projector(T, Persp, Focus) — cxx L70-77 (automatic
    /// minmax directions).
    pub fn from_trsf(t: &Trsf, persp: bool, focus: f64) -> Self {
        let mut p = Projector {
            my_type: 0,
            my_persp: persp,
            my_focus: focus,
            my_scaled_trsf: *t,
            my_trsf: Trsf::identity(),
            my_inv_trsf: Trsf::identity(),
            my_d1: DVec2::ZERO,
            my_d2: DVec2::ZERO,
            my_d3: DVec2::ZERO,
        };
        p.scaled(false);
        p.set_direction();
        p
    }

    /// OCCT HLRAlgo_Projector(T, Persp, Focus, v1, v2, v3) — cxx L81-95
    /// (given minmax directions).
    pub fn from_trsf_dirs(t: &Trsf, persp: bool, focus: f64, v1: DVec2, v2: DVec2, v3: DVec2) -> Self {
        let mut p = Projector {
            my_type: 0,
            my_persp: persp,
            my_focus: focus,
            my_scaled_trsf: *t,
            my_trsf: Trsf::identity(),
            my_inv_trsf: Trsf::identity(),
            my_d1: v1,
            my_d2: v2,
            my_d3: v3,
        };
        p.scaled(false);
        p
    }

    /// OCCT Set(T, Persp, Focus) — cxx L99-106.
    pub fn set(&mut self, t: &Trsf, persp: bool, focus: f64) {
        self.my_persp = persp;
        self.my_focus = focus;
        self.my_scaled_trsf = *t;
        self.scaled(false);
        self.set_direction();
    }

    /// OCCT Scaled(On) — cxx L152-167: to compute with the given scale and
    /// translation.
    pub fn scaled(&mut self, on: bool) {
        self.my_type = -1;
        self.my_trsf = self.my_scaled_trsf;
        if !on {
            self.my_trsf.set_scale_factor(1.0);
            if !self.my_persp {
                self.my_trsf.set_translation_part(DVec3::ZERO);
                self.my_type = trsf_type(&self.my_trsf);
            }
        }
        self.my_inv_trsf = self.my_trsf;
        self.my_inv_trsf.invert();
    }

    /// OCCT Directions(D1, D2, D3) — lxx L23-28.
    pub fn directions(&self) -> (DVec2, DVec2, DVec2) {
        (self.my_d1, self.my_d2, self.my_d3)
    }

    /// OCCT Perspective() — lxx L32-35.
    pub fn perspective(&self) -> bool {
        self.my_persp
    }

    /// OCCT Transformation() — cxx L393-396.
    pub fn transformation(&self) -> &Trsf {
        &self.my_trsf
    }

    /// OCCT InvertedTransformation() — lxx L45-48.
    pub fn inverted_transformation(&self) -> &Trsf {
        &self.my_inv_trsf
    }

    /// OCCT FullTransformation() — lxx L52-55.
    pub fn full_transformation(&self) -> &Trsf {
        &self.my_scaled_trsf
    }

    /// OCCT Focus() — lxx L59-63 (Standard_NoSuchObject if not a
    /// perspective).
    pub fn focus(&self) -> f64 {
        if !self.my_persp {
            panic!("HLRAlgo_Projector::Not a Perspective");
        }
        self.my_focus
    }

    /// OCCT Transform(gp_Vec& D) — lxx L67-70.
    pub fn transform_vec(&self, d: DVec3) -> DVec3 {
        self.my_trsf.transform_vec(d)
    }

    /// OCCT Transform(gp_Pnt& Pnt) — lxx L74-77.
    pub fn transform_pnt(&self, pnt: DVec3) -> DVec3 {
        self.my_trsf.apply(pnt)
    }

    /// OCCT Project(P, Pout) — cxx L171-238.  NOTE: OCCT case 0 computes X/Y
    /// but never writes Pout (kept verbatim).
    pub fn project_pnt(&self, p: DVec3, pout: &mut DVec2) {
        if self.my_type != -1 {
            let x: f64;
            let y: f64;
            match self.my_type {
                0 => {
                    //-- axono standard
                    let x07 = p.x * 0.7071067811865475;
                    let y05 = p.y * 0.5;
                    let z05 = p.z * 0.5;
                    x = x07 - y05 + z05;
                    y = x07 + y05 - z05;
                    //-- Z=0.7071067811865475*(P.Y()+P.Z());
                    // OCCT case 0 never writes Pout.
                    let _ = (x, y);
                    return;
                }
                1 => {
                    //-- top
                    x = p.x;
                    y = p.y; //-- Z=P.Z();
                    *pout = DVec2::new(x, y);
                    return;
                }
                2 => {
                    x = p.x;
                    y = p.z; //-- Z=-P.Y();
                    *pout = DVec2::new(x, y);
                    return;
                }
                3 => {
                    let xmy05 = (p.x - p.y) * 0.5;
                    let z07 = p.z * 0.7071067811865476;
                    x = 0.7071067811865476 * (p.x + p.y);
                    y = -xmy05 + z07;
                    *pout = DVec2::new(x, y);
                    //-- Z= xmy05+z07;
                    return;
                }
                _ => {
                    let p2 = self.transform_pnt(p);
                    if self.my_persp {
                        let r = 1.0 - p2.z / self.my_focus;
                        *pout = DVec2::new(p2.x / r, p2.y / r);
                    } else {
                        *pout = DVec2::new(p2.x, p2.y);
                    }
                    return;
                }
            }
        }
        let p2 = self.transform_pnt(p);
        if self.my_persp {
            let r = 1.0 - p2.z / self.my_focus;
            *pout = DVec2::new(p2.x / r, p2.y / r);
        } else {
            *pout = DVec2::new(p2.x, p2.y);
        }
    }

    /// OCCT Project(P, X, Y, Z) — cxx L262-317.
    pub fn project_xyz(&self, p: DVec3, x: &mut f64, y: &mut f64, z: &mut f64) {
        if self.my_type != -1 {
            match self.my_type {
                0 => {
                    //-- axono standard
                    let x07 = p.x * 0.7071067811865475;
                    let y05 = p.y * 0.5;
                    let z05 = p.z * 0.5;
                    *x = x07 - y05 + z05;
                    *y = x07 + y05 - z05;
                    *z = 0.7071067811865475 * (p.y + p.z);
                }
                1 => {
                    //-- top
                    *x = p.x;
                    *y = p.y;
                    *z = p.z;
                }
                2 => {
                    *x = p.x;
                    *y = p.z;
                    *z = -p.y;
                }
                3 => {
                    let xmy05 = (p.x - p.y) * 0.5;
                    let z07 = p.z * 0.7071067811865476;
                    *x = 0.7071067811865476 * (p.x + p.y);
                    *y = -xmy05 + z07;
                    *z = xmy05 + z07;
                }
                _ => {
                    let p2 = self.transform_pnt(p);
                    *x = p2.x;
                    *y = p2.y;
                    *z = p2.z;
                }
            }
        } else {
            let p2 = self.transform_pnt(p);
            *x = p2.x;
            *y = p2.y;
            *z = p2.z;
            if self.my_persp {
                let r = 1.0 - *z / self.my_focus;
                *x /= r;
                *y /= r;
            }
        }
    }

    /// OCCT Project(P, D1, Pout, D1out) — cxx L321-342.
    pub fn project_pnt_dir(&self, p: DVec3, d1: DVec3, pout: &mut DVec2, d1out: &mut DVec2) {
        let pp = self.my_trsf.apply(p);
        let dd1 = self.my_trsf.transform_vec(d1);
        if self.my_persp {
            let r = 1.0 - pp.z / self.my_focus;
            *pout = DVec2::new(pp.x / r, pp.y / r);
            *d1out = DVec2::new(
                dd1.x / r + pp.x * dd1.z / (self.my_focus * r * r),
                dd1.y / r + pp.y * dd1.z / (self.my_focus * r * r),
            );
        } else {
            *pout = DVec2::new(pp.x, pp.y);
            *d1out = DVec2::new(dd1.x, dd1.y);
        }
    }

    /// OCCT Shoot(X, Y) — cxx L346-359: a line going through the eye towards
    /// the 2d point (X, Y).
    pub fn shoot(&self, x: f64, y: f64) -> Lin {
        let mut l = if self.my_persp {
            Lin::from_pnt_dir(DVec3::new(0.0, 0.0, self.my_focus), DVec3::new(x, y, -self.my_focus))
        } else {
            Lin::from_pnt_dir(DVec3::new(x, y, 0.0), DVec3::new(0.0, 0.0, -1.0))
        };
        l.transform(&self.my_inv_trsf);
        l
    }

    /// OCCT SetDirection() — cxx L363-389.
    fn set_direction(&mut self) {
        let mut v1 = self.transform_vec(DVec3::new(1.0, 0.0, 0.0));
        if (v1.x.abs() + v1.y.abs()) < PRECISION_ANGULAR {
            v1 = DVec3::new(1.0, 1.0, 0.0);
        }
        let d1 = DVec2::new(v1.x, v1.y);
        self.my_d1 = DVec2::new(-d1.y, d1.x);
        let mut v2 = self.transform_vec(DVec3::new(0.0, 1.0, 0.0));
        if (v2.x.abs() + v2.y.abs()) < PRECISION_ANGULAR {
            v2 = DVec3::new(1.0, 1.0, 0.0);
        }
        let d2 = DVec2::new(v2.x, v2.y);
        self.my_d2 = DVec2::new(-d2.y, d2.x);
        let mut v3 = self.transform_vec(DVec3::new(0.0, 0.0, 1.0));
        if (v3.x.abs() + v3.y.abs()) < PRECISION_ANGULAR {
            v3 = DVec3::new(1.0, 1.0, 0.0);
        }
        let d3 = DVec2::new(v3.x, v3.y);
        self.my_d3 = DVec2::new(-d3.y, d3.x);
    }
}

/// OCCT static TrsfType(Trsf) — HLRAlgo_Projector.cxx L112-150: recognizes
/// the standard views by their (scale-1) vectorial part.
fn trsf_type(trsf: &Trsf) -> i32 {
    let mat = trsf.vectorial_part();
    let v = |r: usize, c: usize| mat[r - 1][c - 1];
    if (v(1, 1) - 1.0).abs() < 1e-15
        && (v(2, 2) - 1.0).abs() < 1e-15
        && (v(3, 3) - 1.0).abs() < 1e-15
    {
        1 //-- top
    } else if (v(1, 1) - 0.7071067811865476).abs() < 1e-15
        && (v(1, 2) + 0.5).abs() < 1e-15
        && (v(1, 3) - 0.5).abs() < 1e-15
        && (v(2, 1) - 0.7071067811865476).abs() < 1e-15
        && (v(2, 2) - 0.5).abs() < 1e-15
        && (v(2, 3) + 0.5).abs() < 1e-15
        && v(3, 1).abs() < 1e-15
        && (v(3, 2) - 0.7071067811865476).abs() < 1e-15
        && (v(3, 3) - 0.7071067811865476).abs() < 1e-15
    {
        0 //--
    } else if (v(1, 1) - 1.0).abs() < 1e-15
        && (v(2, 3) - 1.0).abs() < 1e-15
        && (v(3, 2) + 1.0).abs() < 1e-15
    {
        2 //-- front
    } else if (v(1, 1) - 0.7071067811865476).abs() < 1e-15
        && (v(1, 2) - 0.7071067811865476).abs() < 1e-15
        && v(1, 3).abs() < 1e-15
        && (v(2, 1) + 0.5).abs() < 1e-15
        && (v(2, 2) - 0.5).abs() < 1e-15
        && (v(2, 3) - 0.7071067811865476).abs() < 1e-15
        && (v(3, 1) - 0.5).abs() < 1e-15
        && (v(3, 2) + 0.5).abs() < 1e-15
        && (v(3, 3) - 0.7071067811865476).abs() < 1e-15
    {
        3 //-- axo
    } else {
        -1
    }
}
