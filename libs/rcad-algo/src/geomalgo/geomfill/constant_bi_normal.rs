//! OCCT GeomFill_ConstantBiNormal (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_ConstantBiNormal.hxx (members) + GeomFill_ConstantBiNormal.cxx
//! (whole file L34-295, including the file-local FDeriv / DDeriv helpers).

use glam::DVec3;

use rcad_kernel::geom::Curve3;
use rcad_kernel::math::GeomAbsShape;

use super::frenet::Frenet;
use super::trihedron_law::{TrihedronLaw, TrihedronLawBase};

/// OCCT Precision::Confusion().
const CONFUSION: f64 = 1.0e-12;
/// OCCT Precision::Angular().
const ANGULAR: f64 = 1.0e-12;

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx) — |V1 x V2| <= AngularTolerance *
/// |V1| * |V2| (pure math helper, no dedicated OCCT function file).
pub(crate) fn gp_vec_is_parallel(v1: DVec3, v2: DVec3, angular_tolerance: f64) -> bool {
    let cross_mag = v1.cross(v2).length();
    cross_mag <= angular_tolerance * v1.length() * v2.length()
}

/// OCCT gp_Vec::IsNormal (gp_Vec.hxx) — |V1 . V2| <= AngularTolerance *
/// |V1| * |V2| (pure math helper, no dedicated OCCT function file).
pub(crate) fn gp_vec_is_normal(v1: DVec3, v2: DVec3, angular_tolerance: f64) -> bool {
    let dot = v1.dot(v2).abs();
    dot <= angular_tolerance * v1.length() * v2.length()
}

/// OCCT GeomFill_ConstantBiNormal (GeomFill_ConstantBiNormal.hxx L99-102):
/// BN + the inherited GeomFill_TrihedronLaw base + the Frenet handle.
#[derive(Debug, Clone)]
pub struct ConstantBiNormal {
    /// OCCT GeomFill_TrihedronLaw protected base (myCurve / myTrimmed).
    pub(crate) base: TrihedronLawBase,
    /// OCCT gp_Vec BN.
    bn: DVec3,
    /// OCCT handle(GeomFill_Frenet) frenet.
    frenet: Frenet,
}

/// OCCT static FDeriv (L32-43) — computes (F/|F|)'.
fn f_deriv(f: DVec3, df: DVec3) -> DVec3 {
    let norma = f.length();
    (df - f * (f.dot(df)) / (norma * norma)) / norma
}

/// OCCT static DDeriv (L45-54) — computes (F/|F|)''.
fn d_deriv(f: DVec3, df: DVec3, d2f: DVec3) -> DVec3 {
    let norma = f.length();
    (d2f - 2.0 * df * (f.dot(df)) / (norma * norma)) / norma
        - f * ((df.dot(df) + f.dot(d2f) - 3.0 * (f.dot(df)) * (f.dot(df)) / (norma * norma))
            / (norma * norma * norma))
}

impl ConstantBiNormal {
    /// OCCT GeomFill_ConstantBiNormal::GeomFill_ConstantBiNormal (L56-60).
    pub fn new(bi_normal: DVec3) -> Self {
        ConstantBiNormal {
            base: TrihedronLawBase::default(),
            bn: bi_normal,
            frenet: Frenet::new(),
        }
    }
}

impl TrihedronLaw for ConstantBiNormal {
    fn my_curve(&self) -> &Option<Curve3> {
        &self.base.my_curve
    }

    fn my_trimmed(&self) -> &Option<Curve3> {
        &self.base.my_trimmed
    }

    fn set_my_curve(&mut self, c: Curve3) {
        self.base.my_curve = Some(c);
    }

    fn set_my_trimmed(&mut self, c: Option<Curve3>) {
        self.base.my_trimmed = c;
    }

    /// OCCT Copy (L62-70).
    fn copy_law(&self) -> Box<dyn TrihedronLaw> {
        let mut copy = ConstantBiNormal::new(self.bn);
        if let Some(curve) = &self.base.my_curve {
            TrihedronLaw::set_curve(&mut copy, curve.clone());
        }
        Box::new(copy)
    }

    /// OCCT SetCurve (L72-81).
    fn set_curve(&mut self, c: Curve3) -> bool {
        super::trihedron_law::trihedron_law_base_set_curve(self, c.clone());
        // OCCT: if (!C.IsNull()) { isOK = frenet->SetCurve(C); } — the rcad
        // Curve3 mapping has no null state, so the branch always holds.
        self.frenet.set_curve(c)
    }

    /// OCCT D0 (L83-109).
    fn d0(
        &self,
        param: f64,
        tangent: &mut DVec3,
        normal: &mut DVec3,
        binormal: &mut DVec3,
    ) -> bool {
        // if BN^T != 0 then N = (BN^T).Normalized ; T = N^BN
        // else T = (N^BN).Normalized ; N = BN^T

        self.frenet.d0(param, tangent, normal, binormal);
        *binormal = self.bn;
        if binormal.cross(*tangent).length() > CONFUSION {
            *normal = binormal.cross(*tangent).normalize_or_zero();
            *tangent = normal.cross(*binormal);
        } else {
            *tangent = normal.cross(*binormal).normalize_or_zero();
            *normal = binormal.cross(*tangent);
        }
        true
    }

    /// OCCT D1 (L111-160).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        tangent: &mut DVec3,
        dtangent: &mut DVec3,
        normal: &mut DVec3,
        dnormal: &mut DVec3,
        binormal: &mut DVec3,
        dbinormal: &mut DVec3,
    ) -> bool {
        self.frenet.d1(
            param, tangent, dtangent, normal, dnormal, binormal, dbinormal,
        );
        *binormal = self.bn;
        *dbinormal = DVec3::ZERO;
        if binormal.cross(*tangent).length() > CONFUSION {
            let f = binormal.cross(*tangent);
            let df = binormal.cross(*dtangent);
            *normal = f.normalize_or_zero();
            *dnormal = f_deriv(f, df);

            *tangent = normal.cross(*binormal);
            *dtangent = dnormal.cross(*binormal);
        } else {
            let f = normal.cross(*binormal);
            let df = dnormal.cross(*binormal);
            *tangent = f.normalize_or_zero();
            *dtangent = f_deriv(f, df);

            *normal = binormal.cross(*tangent);
            *dnormal = binormal.cross(*dtangent);
        }
        true
    }

    /// OCCT D2 (L162-225).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
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
        self.frenet.d2(
            param,
            tangent,
            dtangent,
            d2tangent,
            normal,
            dnormal,
            d2normal,
            binormal,
            dbinormal,
            d2binormal,
        );
        *binormal = self.bn;
        *dbinormal = DVec3::ZERO;
        *d2binormal = DVec3::ZERO;
        if binormal.cross(*tangent).length() > CONFUSION {
            let f = binormal.cross(*tangent);
            let df = binormal.cross(*dtangent);
            let d2f = binormal.cross(*d2tangent);
            *normal = f.normalize_or_zero();
            *dnormal = f_deriv(f, df);
            *d2normal = d_deriv(f, df, d2f);

            *tangent = normal.cross(*binormal);
            *dtangent = dnormal.cross(*binormal);
            *d2tangent = d2normal.cross(*binormal);
        } else {
            let f = normal.cross(*binormal);
            let df = dnormal.cross(*binormal);
            let d2f = d2normal.cross(*binormal);
            *tangent = f.normalize_or_zero();
            *dtangent = f_deriv(f, df);
            *d2tangent = d_deriv(f, df, d2f);

            *normal = binormal.cross(*tangent);
            *dnormal = binormal.cross(*dtangent);
            *d2normal = binormal.cross(*d2tangent);
        }
        true
    }

    /// OCCT NbIntervals (L227-230).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        TrihedronLaw::nb_intervals(&self.frenet, s)
    }

    /// OCCT Intervals (L232-236).
    fn intervals(&self, t: &mut Vec<f64>, s: GeomAbsShape) {
        TrihedronLaw::intervals(&self.frenet, t, s);
    }

    /// OCCT GetAverageLaw (L238-252).
    fn get_average_law(&self, atangent: &mut DVec3, anormal: &mut DVec3, abinormal: &mut DVec3) {
        self.frenet
            .get_average_law(atangent, anormal, abinormal);
        *abinormal = self.bn;
        if abinormal.cross(*atangent).length() > CONFUSION {
            *anormal = abinormal.cross(*atangent).normalize_or_zero();
            *atangent = anormal.cross(*abinormal);
        } else {
            *atangent = anormal.cross(*abinormal).normalize_or_zero();
            *anormal = abinormal.cross(*atangent);
        }
    }

    /// OCCT IsConstant (L254-257).
    fn is_constant(&self) -> bool {
        TrihedronLaw::is_constant(&self.frenet)
    }

    /// OCCT IsOnlyBy3dCurve (L259-295).
    fn is_only_by3d_curve(&self) -> bool {
        let my_curve = self.base.my_curve.as_ref().expect("null myCurve");

        // OCCT GeomAbs_CurveType TheType = myCurve->GetType(); the switch
        // assigns TheAxe per analytic type; the Line case returns early and
        // the default returns false (Rust match carries the same branches).
        let the_axe = match my_curve {
            Curve3::Circle(circle) => {
                // TheAxe = myCurve->Circle().Axis().
                circle.normal
            }
            Curve3::Ellipse(ellipse) => {
                // TheAxe = myCurve->Ellipse().Axis().
                ellipse.normal
            }
            Curve3::Hyperbola(hyperbola) => {
                // TheAxe = myCurve->Hyperbola().Axis().
                hyperbola.normal
            }
            Curve3::Parabola(parabola) => {
                // TheAxe = myCurve->Parabola().Axis().
                parabola.axis_dir
            }
            // La normale du plan de la courbe est il perpendiculaire a la
            // BiNormale ?
            Curve3::Line(line) => {
                let v = line.direction;
                return gp_vec_is_normal(v, self.bn, ANGULAR);
            }
            _ => return false, // pas de risques
        };

        // La normale du plan de la courbe est il // a la BiNormale ?
        gp_vec_is_parallel(the_axe, self.bn, ANGULAR)
    }
}
