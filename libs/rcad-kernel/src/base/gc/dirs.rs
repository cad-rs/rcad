//! OCCT gce_MakeDir (TKGeomBase/gce/gce_MakeDir.cxx L27-80) — the unit
//! direction constructors.  A gp_Dir is rcad's normalized `Vec3`
//! (= `DVec3`).

use glam::DVec3;

use super::GceError;

/// OCCT gp::Resolution() (gp.hxx) — DBL_MIN.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT gce_MakeDir::gce_MakeDir(const gp_Pnt& P1, const gp_Pnt& P2)
/// (gce_MakeDir.cxx L29-41): the unit direction P1 -> P2; ConfusedPoints
/// when the distance is at most gp::Resolution().
pub fn make_dir_2p(p1: DVec3, p2: DVec3) -> Result<DVec3, GceError> {
    if p1.distance(p2) <= GP_RESOLUTION {
        Err(GceError::ConfusedPoints)
    } else {
        // OCCT L37: TheDir = gp_Dir(P2.XYZ() - P1.XYZ()) — the normalized
        // difference vector.
        Ok((p2 - p1).normalize_or_zero())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The OCCT (P1, P2) ctor: a unit direction P1 -> P2; ConfusedPoints at
    // a sub-Resolution distance.
    #[test]
    fn make_dir_2p_unit_and_confused() {
        let d = make_dir_2p(DVec3::ZERO, DVec3::new(3.0, 0.0, 0.0)).unwrap();
        assert!((d - DVec3::X).length() < 1e-15);
        // The ConfusedPoints arm: identical points.
        assert!(matches!(
            make_dir_2p(DVec3::ZERO, DVec3::ZERO),
            Err(GceError::ConfusedPoints)
        ));
    }
}
