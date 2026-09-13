//! OCCT GeomLib_MakeCurvefromApprox (TKGeomBase/GeomLib) — the BSpline curve
//! construction from an approximation.
//!
//! 1:1 translation of `GeomLib_MakeCurvefromApprox.hxx` (L25-89),
//! `GeomLib_MakeCurvefromApprox.lxx` (L24-27, the inline IsDone) and
//! `GeomLib_MakeCurvefromApprox.cxx` (L29-203).
//!
//! The OCCT class holds the `AdvApprox_ApproxAFunction` by value; the rcad
//! translation holds the kernel `ApproxAFunction` by reference for the same
//! reason the OCCT member is a value copy of the argument (the accessors read
//! through to the same approximation state).
//!
//! Architecture difference: the rcad `BSplineCurve3` stores the full
//! (multiplicity-expanded) knot vector, while the OCCT `Geom_BSplineCurve`
//! constructor takes the distinct knots plus multiplicities; the expansion is
//! the representation change of that mapping.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{BSplineCurve2, BSplineCurve3};
use rcad_kernel::math::adv_approx::ApproxAFunction;

/// OCCT GeomLib_MakeCurvefromApprox (hxx L25-89) — this class is used to
/// construct the BSpline curve from an Approximation (ApproxAFunction from
/// AdvApprox).
pub struct GeomLibMakeCurvefromApprox<'a> {
    /// OCCT hxx L87: AdvApprox_ApproxAFunction myApprox.
    my_approx: &'a ApproxAFunction,
}

impl<'a> GeomLibMakeCurvefromApprox<'a> {
    /// OCCT GeomLib_MakeCurvefromApprox(const AdvApprox_ApproxAFunction&)
    /// (cxx L33-36).
    pub fn new(approx: &'a ApproxAFunction) -> Self {
        GeomLibMakeCurvefromApprox { my_approx: approx }
    }

    /// OCCT GeomLib_MakeCurvefromApprox::IsDone (lxx L24-27).
    pub fn is_done(&self) -> bool {
        self.my_approx.is_done()
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Nb1DSpaces (cxx L42-45).
    pub fn nb1d_spaces(&self) -> i32 {
        self.my_approx.num_sub_spaces_of(1) as i32
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Nb2DSpaces (cxx L51-54).
    pub fn nb2d_spaces(&self) -> i32 {
        self.my_approx.num_sub_spaces_of(2) as i32
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Nb3DSpaces (cxx L60-63).
    pub fn nb3d_spaces(&self) -> i32 {
        self.my_approx.num_sub_spaces_of(3) as i32
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Curve2d(Index2d) (cxx L68-85) —
    /// returns a polynomial curve whose poles correspond to the Index2d 2D
    /// space.
    pub fn curve2d(&self, index2d: i32) -> Option<BSplineCurve2> {
        if index2d < 0 || index2d > self.nb2d_spaces() {
            panic!("Standard_OutOfRange: GeomLib_MakeCurvefromApprox : Curve2d");
        }
        if !self.is_done() {
            panic!("StdFail_NotDone: GeomLib_MakeCurvefromApprox : Curve2d");
        }

        let nb_poles = self.my_approx.nb_poles();
        let poles_flat = self.my_approx.poles2d_flat(index2d as usize);
        let poles: Vec<DVec2> = (0..nb_poles)
            .map(|i| DVec2::new(poles_flat[i * 2], poles_flat[i * 2 + 1]))
            .collect();

        Some(self.make_curve2d(poles))
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Curve2d(Index1d, Index2d) (cxx
    /// L91-126) — returns a rational curve whose poles correspond to the
    /// index2d of the 2D space and whose weights correspond to one
    /// dimensional space of index 1d.
    pub fn curve2d_1d_2d(&self, index1d: i32, index2d: i32) -> Option<BSplineCurve2> {
        if index1d < 0 || index1d > self.nb1d_spaces() || index2d < 0 || index2d > self.nb2d_spaces()
        {
            panic!("Standard_OutOfRange: GeomLib_MakeCurvefromApprox : Curve2d");
        }
        if !self.is_done() && !self.my_approx.has_result() {
            panic!("StdFail_NotDone: GeomLib_MakeCurvefromApprox : Curve2d");
        }

        let nb_poles = self.my_approx.nb_poles();
        let mut poles_flat = self.my_approx.poles2d_flat(index2d as usize);
        let weights = self.my_approx.poles1d_flat(index1d as usize);

        // OCCT: Poles(i).SetCoord(X / W, Y / W).
        for i in 0..nb_poles {
            let w = weights[i];
            poles_flat[i * 2] /= w;
            poles_flat[i * 2 + 1] /= w;
        }
        let poles: Vec<DVec2> = (0..nb_poles)
            .map(|i| DVec2::new(poles_flat[i * 2], poles_flat[i * 2 + 1]))
            .collect();

        Some(self.make_curve2d(poles))
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Curve2dFromTwo1d(Index1d, Index2d)
    /// (cxx L132-171) — returns a 2D curve building it from the 1D curve in
    /// x at Index1d and y at Index2d amongst the 1D curves.
    pub fn curve2d_from_two1d(&self, index1d: i32, index2d: i32) -> Option<BSplineCurve2> {
        if index1d < 0 || index1d > self.nb1d_spaces() || index2d < 0 || index2d > self.nb1d_spaces()
        {
            panic!("Standard_OutOfRange: GeomLib_MakeCurvefromApprox : Curve2d");
        }
        if !self.is_done() && !self.my_approx.has_result() {
            panic!("StdFail_NotDone: GeomLib_MakeCurvefromApprox : Curve2d");
        }

        let nb_poles = self.my_approx.nb_poles();
        let poles1d1 = self.my_approx.poles1d_flat(index1d as usize);
        let poles1d2 = self.my_approx.poles1d_flat(index2d as usize);

        // OCCT: Poles(i).SetCoord(Poles1d1.Value(i), Poles1d2.Value(i)).
        let poles: Vec<DVec2> = (0..nb_poles)
            .map(|i| DVec2::new(poles1d1[i], poles1d2[i]))
            .collect();

        Some(self.make_curve2d(poles))
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Curve(Index3d) (cxx L177-195) —
    /// returns a polynomial curve whose poles correspond to the Index3D 3D
    /// space.
    pub fn curve(&self, index3d: i32) -> Option<BSplineCurve3> {
        if index3d < 0 || index3d > self.nb3d_spaces() {
            panic!("Standard_OutOfRange: GeomLib_MakeCurvefromApprox : Curve");
        }
        if !self.is_done() && !self.my_approx.has_result() {
            panic!("StdFail_NotDone: GeomLib_MakeCurvefromApprox : Curve");
        }

        let nb_poles = self.my_approx.nb_poles();
        let poles_flat = self.my_approx.poles_flat(index3d as usize);
        let poles: Vec<DVec3> = (0..nb_poles)
            .map(|i| {
                DVec3::new(
                    poles_flat[i * 3],
                    poles_flat[i * 3 + 1],
                    poles_flat[i * 3 + 2],
                )
            })
            .collect();

        Some(self.make_curve3d(poles))
    }

    /// OCCT GeomLib_MakeCurvefromApprox::Curve(Index1d, Index3d) (cxx
    /// L201-241) — returns a rational curve whose poles correspond to the
    /// index3D of the 3D space and whose weights correspond to the index1d 1D
    /// space.
    pub fn curve_1d_3d(&self, index1d: i32, index3d: i32) -> Option<BSplineCurve3> {
        if index1d < 0 || index1d > self.nb1d_spaces() || index3d < 0 || index3d > self.nb3d_spaces()
        {
            panic!("Standard_OutOfRange: GeomLib_MakeCurvefromApprox : Curve3d");
        }
        if !self.is_done() {
            panic!("StdFail_NotDone: GeomLib_MakeCurvefromApprox : Curve3d");
        }

        let nb_poles = self.my_approx.nb_poles();
        let mut poles_flat = self.my_approx.poles_flat(index3d as usize);
        let weights = self.my_approx.poles1d_flat(index1d as usize);

        // OCCT: Poles(i).SetCoord(X / W, Y / W, Z / W).
        for i in 0..nb_poles {
            let w = weights[i];
            poles_flat[i * 3] /= w;
            poles_flat[i * 3 + 1] /= w;
            poles_flat[i * 3 + 2] /= w;
        }
        let poles: Vec<DVec3> = (0..nb_poles)
            .map(|i| {
                DVec3::new(
                    poles_flat[i * 3],
                    poles_flat[i * 3 + 1],
                    poles_flat[i * 3 + 2],
                )
            })
            .collect();

        Some(self.make_curve3d(poles))
    }

    /// OCCT: new Geom2d_BSplineCurve(Poles, Knots, Mults, Degree)
    /// (cxx L83-84 / L123-124 / L168-169).
    fn make_curve2d(&self, poles: Vec<DVec2>) -> BSplineCurve2 {
        BSplineCurve2 {
            degree: self.my_approx.degree() as usize,
            knots: full_knots(
                self.my_approx.knots_vec(),
                self.my_approx.multiplicities_vec(),
            ),
            control_points: poles,
            weights: vec![1.0; self.my_approx.nb_poles()],
        }
    }

    /// OCCT: new Geom_BSplineCurve(Poles, Knots, Mults, Degree)
    /// (cxx L193 / L239).
    fn make_curve3d(&self, poles: Vec<DVec3>) -> BSplineCurve3 {
        BSplineCurve3 {
            degree: self.my_approx.degree() as usize,
            knots: full_knots(
                self.my_approx.knots_vec(),
                self.my_approx.multiplicities_vec(),
            ),
            control_points: poles,
            weights: vec![1.0; self.my_approx.nb_poles()],
            is_periodic: false,
        }
    }
}

/// Expand the distinct knots + multiplicities into the full knot vector the
/// rcad BSplineCurve representation stores (the OCCT constructor consumes the
/// pair directly).
fn full_knots(knots: &[f64], mults: &[i32]) -> Vec<f64> {
    let mut out = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            out.push(*k);
        }
    }
    out
}
