//! OCCT GeomFill_Fixed (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_Fixed.hxx (members) + GeomFill_Fixed.cxx (whole file L26-119).

use glam::DVec3;

use rcad_kernel::math::GeomAbsShape;

use super::trihedron_law::{TrihedronLaw, TrihedronLawBase};

/// OCCT GeomFill_Fixed (GeomFill_Fixed.hxx L100-106): the fixed members
/// T, N, B plus the inherited GeomFill_TrihedronLaw base.
#[derive(Debug, Clone)]
pub struct Fixed {
    /// OCCT GeomFill_TrihedronLaw protected base (myCurve / myTrimmed).
    pub(crate) base: TrihedronLawBase,
    /// OCCT gp_Vec T.
    t: DVec3,
    /// OCCT gp_Vec N.
    n: DVec3,
    /// OCCT gp_Vec B.
    b: DVec3,
}

impl Fixed {
    /// OCCT GeomFill_Fixed::GeomFill_Fixed (L26-39).
    pub fn new(tangent: DVec3, normal: DVec3) -> Self {
        if super::constant_bi_normal::gp_vec_is_parallel(tangent, normal, 0.01) {
            panic!("Standard_ConstructionError: GeomFill_Fixed : Two parallel vectors !");
        }
        let mut t = tangent;
        t = t.normalize_or_zero();
        let mut n = normal;
        n = n.normalize_or_zero();
        let mut b = t.cross(n);
        b = b.normalize_or_zero();
        Fixed {
            base: TrihedronLawBase::default(),
            t,
            n,
            b,
        }
    }
}

impl TrihedronLaw for Fixed {
    fn my_curve(&self) -> &Option<rcad_kernel::geom::Curve3> {
        &self.base.my_curve
    }

    fn my_trimmed(&self) -> &Option<rcad_kernel::geom::Curve3> {
        &self.base.my_trimmed
    }

    fn set_my_curve(&mut self, c: rcad_kernel::geom::Curve3) {
        self.base.my_curve = Some(c);
    }

    fn set_my_trimmed(&mut self, c: Option<rcad_kernel::geom::Curve3>) {
        self.base.my_trimmed = c;
    }

    /// OCCT Copy (L40-46).
    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let mut copy = Fixed::new(self.t, self.n);
        if let Some(curve) = &self.base.my_curve {
            TrihedronLaw::set_curve(&mut copy, curve.clone());
        }
        Box::new(copy)
    }

    /// OCCT D0 (L47-55).
    fn d0(
        &self,
        _param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool {
        *tangent = self.t;
        *normal = self.n;
        *binormal = self.b;
        true
    }

    /// OCCT D1 (L56-73).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        _param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
    ) -> bool {
        *tangent = self.t;
        *normal = self.n;
        *binormal = self.b;

        // OCCT gp_Vec V0(0, 0, 0); DTangent = DNormal = DBiNormal = V0.
        *dtangent = DVec3::ZERO;
        *dnormal = DVec3::ZERO;
        *dbinormal = DVec3::ZERO;
        true
    }

    /// OCCT D2 (L74-96).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        _param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        d2tangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        d2normal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
        d2binormal: &mut DVec3,
    ) -> bool {
        *tangent = self.t;
        *normal = self.n;
        *binormal = self.b;

        // OCCT gp_Vec V0(0, 0, 0); DTangent = D2Tangent = V0;
        // DNormal = D2Normal = V0; DBiNormal = D2BiNormal = V0.
        *dtangent = DVec3::ZERO;
        *d2tangent = DVec3::ZERO;
        *dnormal = DVec3::ZERO;
        *d2normal = DVec3::ZERO;
        *dbinormal = DVec3::ZERO;
        *d2binormal = DVec3::ZERO;
        true
    }

    /// OCCT NbIntervals (L97-101).
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// OCCT Intervals (L102-107).
    fn intervals(&self, t: &mut Vec<f64>, _s: GeomAbsShape) {
        let n = t.len();
        t[0] = f64::NEG_INFINITY;
        t[n - 1] = f64::INFINITY;
    }

    /// OCCT GetAverageLaw (L108-114).
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3) {
        *atangent = self.t;
        *anormal = self.n;
        *abinormal = self.b;
    }

    /// OCCT IsConstant (L115-119).
    fn is_constant(&self) -> bool {
        true
    }
}
