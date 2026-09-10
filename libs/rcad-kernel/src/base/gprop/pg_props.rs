//! OCCT GProp_PGProps (TKGeomBase GProp package).
//!
//! 1:1 translation of `GProp_PGProps.cxx` (L24-221, the 1D-array and
//! static-barycentre members used by GProp_PEquation) over the stored
//! mass/inertia state of GProp_GProps.  The 2D-array overloads are not
//! ported (no consumer in rcad yet).
//!
//! OCCT keeps `g`, `dim`, `inertia` in the GProp_GProps base; the same
//! members live here.

use glam::DVec3;

/// OCCT gp_Mat as a row-major 3x3.
type Mat3 = [[f64; 3]; 3];

/// OCCT GProp_PGProps — computation of the global properties of a point
/// cloud (centre of mass and inertia matrix).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GPropPGProps {
    g: DVec3,
    dim: f64,
    inertia: Mat3,
}

impl Default for GPropPGProps {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT gp_Mat(a11..a13, a21..a23, a31..a33) rows.
fn mat3(r1: [f64; 3], r2: [f64; 3], r3: [f64; 3]) -> Mat3 {
    [r1, r2, r3]
}

impl GPropPGProps {
    /// OCCT GProp_PGProps() — cxx L24-30.
    pub fn new() -> Self {
        GPropPGProps {
            g: DVec3::ZERO,
            dim: 0.0,
            inertia: mat3([0.0; 3], [0.0; 3], [0.0; 3]),
        }
    }

    /// OCCT GProp_PGProps(thePnts) — cxx L34-41.
    pub fn from_points(the_pnts: &[DVec3]) -> Self {
        let mut props = GPropPGProps::new();
        for p in the_pnts {
            props.add_point(*p);
        }
        props
    }

    /// OCCT GProp_PGProps(thePnts, theDensity) — cxx L59-73.
    pub fn from_points_density(the_pnts: &[DVec3], the_density: &[f64]) -> Self {
        if the_pnts.len() != the_density.len() {
            panic!("GProp_PGProps: points and density arrays have different lengths");
        }
        let mut props = GPropPGProps::new();
        for (p, d) in the_pnts.iter().zip(the_density.iter()) {
            props.add_point_density(*p, *d);
        }
        props
    }

    /// OCCT AddPoint(thePnt) — cxx L102-128.
    pub fn add_point(&mut self, the_pnt: DVec3) {
        let (a_x, a_y, a_z) = (the_pnt.x, the_pnt.y, the_pnt.z);

        let a_mp = mat3(
            [a_y * a_y + a_z * a_z, -a_x * a_y, -a_x * a_z],
            [-a_x * a_y, a_x * a_x + a_z * a_z, -a_y * a_z],
            [-a_x * a_z, -a_y * a_z, a_x * a_x + a_y * a_y],
        );

        if self.dim == 0.0 {
            self.dim = 1.0;
            self.g = the_pnt;
            self.inertia = a_mp;
        } else {
            let a_new_mass = self.dim + 1.0;
            let a_new_x = (self.g.x * self.dim + a_x) / a_new_mass;
            let a_new_y = (self.g.y * self.dim + a_y) / a_new_mass;
            let a_new_z = (self.g.z * self.dim + a_z) / a_new_mass;
            self.g = DVec3::new(a_new_x, a_new_y, a_new_z);
            self.dim = a_new_mass;
            for (row, mrow) in self.inertia.iter_mut().zip(a_mp.iter()) {
                for (v, mv) in row.iter_mut().zip(mrow.iter()) {
                    *v += mv;
                }
            }
        }
    }

    /// OCCT AddPoint(thePnt, theDensity) — cxx L132-163.
    pub fn add_point_density(&mut self, the_pnt: DVec3, the_density: f64) {
        if the_density <= 0.0 {
            panic!("GProp_PGProps::AddPoint: density must be positive");
        }

        let (a_x, a_y, a_z) = (the_pnt.x, the_pnt.y, the_pnt.z);

        let a_mp = mat3(
            [a_y * a_y + a_z * a_z, -a_x * a_y, -a_x * a_z],
            [-a_x * a_y, a_x * a_x + a_z * a_z, -a_y * a_z],
            [-a_x * a_z, -a_y * a_z, a_x * a_x + a_y * a_y],
        );

        if self.dim == 0.0 {
            self.dim = the_density;
            self.g = the_pnt;
            self.inertia = a_mp;
            for row in self.inertia.iter_mut() {
                for v in row.iter_mut() {
                    *v *= the_density;
                }
            }
        } else {
            let a_new_mass = self.dim + the_density;
            let a_new_x = (self.g.x * self.dim + a_x * the_density) / a_new_mass;
            let a_new_y = (self.g.y * self.dim + a_y * the_density) / a_new_mass;
            let a_new_z = (self.g.z * self.dim + a_z * the_density) / a_new_mass;
            self.g = DVec3::new(a_new_x, a_new_y, a_new_z);
            self.dim = a_new_mass;
            // inertia = inertia + aMp * theDensity.
            for (row, mrow) in self.inertia.iter_mut().zip(a_mp.iter()) {
                for (v, mv) in row.iter_mut().zip(mrow.iter()) {
                    *v += mv * the_density;
                }
            }
        }
    }

    /// OCCT CentreOfMass() (GProp_GProps accessor) — the stored barycentre.
    pub fn centre_of_mass(&self) -> DVec3 {
        self.g
    }

    /// OCCT MatrixOfInertia() — GProp_GProps.cxx L110-115:
    /// `HOperator(g, Origin, dim, HMat); return inertia - HMat;` — the
    /// Huygens transfer of the raw (origin-referenced) inertia to the
    /// centre of mass.  GProp::HOperator (GProp.cxx L20-39) builds
    /// `Mass * (|QG|^2*I - QG QG^T)` with QG = G - Q.
    pub fn matrix_of_inertia(&self) -> Mat3 {
        // GProp::HOperator(G = g, Q = Origin).
        let qg = self.g;
        let ixx = qg.y * qg.y + qg.z * qg.z;
        let iyy = qg.x * qg.x + qg.z * qg.z;
        let izz = qg.y * qg.y + qg.x * qg.x;
        let ixy = -qg.x * qg.y;
        let iyz = -qg.y * qg.z;
        let ixz = -qg.x * qg.z;
        // SetCols(Ixx Ixy Ixz, Ixy Iyy Iyz, Ixz Iyz Izz) * Mass.
        // Column-major SetCols maps to the row-major matrix below.
        let h_mat = [
            [ixx * self.dim, ixy * self.dim, ixz * self.dim],
            [ixy * self.dim, iyy * self.dim, iyz * self.dim],
            [ixz * self.dim, iyz * self.dim, izz * self.dim],
        ];
        let mut result = self.inertia;
        for (row, hrow) in result.iter_mut().zip(h_mat.iter()) {
            for (v, hv) in row.iter_mut().zip(hrow.iter()) {
                *v -= hv;
            }
        }
        result
    }

    /// OCCT Mass() (GProp_GProps accessor).
    pub fn mass(&self) -> f64 {
        self.dim
    }

    /// OCCT static Barycentre(thePnts) — cxx L167-176.
    pub fn barycentre(the_pnts: &[DVec3]) -> DVec3 {
        let mut a_sum = the_pnts[0];
        for p in &the_pnts[1..] {
            a_sum += *p;
        }
        a_sum /= the_pnts.len() as f64;
        a_sum
    }

    /// OCCT static Barycentre(thePnts, theDensity, theMass, theG) — cxx
    /// L196-221.
    pub fn barycentre_density(the_pnts: &[DVec3], the_density: &[f64]) -> (f64, DVec3) {
        if the_pnts.len() != the_density.len() {
            panic!("GProp_PGProps::Barycentre: points and density arrays have different lengths");
        }

        let mut the_mass = the_density[0];
        let mut a_sum = the_pnts[0] * the_mass;
        for i in 1..the_pnts.len() {
            the_mass += the_density[i];
            a_sum += the_pnts[i] * the_density[i];
        }
        a_sum /= the_mass;
        (the_mass, a_sum)
    }
}
