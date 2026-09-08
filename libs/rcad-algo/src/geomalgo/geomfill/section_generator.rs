//! OCCT GeomFill_SectionGenerator (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_SectionGenerator.hxx (L36-73) + GeomFill_SectionGenerator.cxx
//! (whole file L22-109).
//!
//! Architecture mapping: OCCT inherits GeomFill_Profiler; Rust composes the
//! base (the profiler translation lives in `super::profiler`).

use glam::{DVec2, DVec3};

use super::profiler::Profiler;

/// OCCT GeomFill_SectionGenerator — gives the functions needed for
/// instantiation from AppSurf in AppBlend. Allows to evaluate a surface
/// passing by all the curves if the Profiler.
pub struct SectionGenerator {
    /// OCCT base class GeomFill_Profiler.
    pub base: Profiler,
    /// OCCT myParams (handle(NCollection_HArray1<double>), 1-based values).
    my_params: Option<Vec<f64>>,
}

impl SectionGenerator {
    /// OCCT GeomFill_SectionGenerator::GeomFill_SectionGenerator (L22-34).
    pub fn new() -> Self {
        let mut section_gen = SectionGenerator {
            base: Profiler::new(),
            my_params: None,
        };
        // OCCT: if (mySequence.Length() > 1) { HPar(i) = i - 1; SetParam(HPar); }
        if section_gen.base.my_sequence.len() > 1 {
            let h_par: Vec<f64> =
                (0..section_gen.base.my_sequence.len()).map(|i| i as f64).collect();
            section_gen.set_param(h_par);
        }
        section_gen
    }

    /// OCCT SetParam (L38-46).
    pub fn set_param(&mut self, params: Vec<f64>) {
        let l = params.len();
        // OCCT: myParams = Params; then re-numbers the 1-based values from
        // Params->Lower()..Upper() into 1..L.
        let mut my_params = vec![0.0f64; l];
        for ii in 1..=l {
            my_params[ii - 1] = params[ii - 1];
        }
        self.my_params = Some(my_params);
    }

    /// OCCT GetShape (L50-60).
    pub fn get_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32, nb_poles2d: &mut i32) {
        let c = &self.base.my_sequence[0];
        *nb_poles = c.control_points.len() as i32;
        *nb_knots = c.knots_mults().0.len() as i32;
        *degree = c.degree as i32;
        *nb_poles2d = 0;
    }

    /// OCCT Knots (L64-67).
    pub fn knots(&self, t_knots: &mut [f64]) {
        let (knots, _) = self.base.my_sequence[0].knots_mults();
        t_knots.copy_from_slice(&knots);
    }

    /// OCCT Mults (L71-74).
    pub fn mults(&self, t_mults: &mut [i32]) {
        let (_, mults) = self.base.my_sequence[0].knots_mults();
        t_mults.copy_from_slice(&mults);
    }

    /// OCCT Section (L78-89) — used for the first and last section. The
    /// method returns true if the derivatives are computed, otherwise it
    /// returns false.
    #[allow(clippy::too_many_arguments)]
    pub fn section_d1(
        &self,
        p: i32,
        poles: &mut [DVec3],
        _dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        _dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        _dweigths: &mut [f64],
    ) -> bool {
        self.section(p, poles, poles2d, weigths);
        false
    }

    /// OCCT Section (L93-102).
    pub fn section(&self, p: i32, poles: &mut [DVec3], _poles2d: &mut [DVec2], weigths: &mut [f64]) {
        let c = &self.base.my_sequence[(p - 1) as usize];
        poles.copy_from_slice(&c.control_points);
        weigths.copy_from_slice(&c.weights);
    }

    /// OCCT Parameter (L106-109) — returns the parameter of Section<P>, to
    /// impose it for the approximation.
    pub fn parameter(&self, p: i32) -> f64 {
        self.my_params.as_ref().expect("GeomFill_SectionGenerator::Parameter")[(p - 1) as usize]
    }
}

impl Default for SectionGenerator {
    fn default() -> Self {
        Self::new()
    }
}
