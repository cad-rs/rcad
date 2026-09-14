//! OCCT Geom_OsculatingSurface (TKG3d/Geom) — 1:1 translation of
//! `Geom_OsculatingSurface.pxx` (L38-130) + `Geom_OsculatingSurface.cxx`
//! (L32-839).
//!
//! Internal helper class for Geom_OffsetSurface: it detects whether a
//! (BSpline / Bezier) surface has a punctual U or V isoparametric curve on
//! one of its bounds within a tolerance, and builds the corresponding
//! osculating surfaces.
//!
//! Architecture differences (rcad value model):
//!   - `occ::handle<Geom_Surface> myBasisSurf` maps to the kernel surface
//!     value [`Surface3`] (`theBS->Copy()` of `Init` L111 is the rcad clone);
//!     the OCCT default-constructed null handle is `None`;
//!   - `NCollection_Sequence<handle(Geom_BSplineSurface)> myOsculSurf1/2`
//!     map to `Vec<BSplineSurface>` and `NCollection_Sequence<int> myKdeg`
//!     to `Vec<i32>`;
//!   - `std::array<bool, 4> myAlong` maps to `[bool; 4]`;
//!   - the OCCT copy/move ctors and assignment operators (pxx L49-101) are
//!     the derived Rust `Clone` (the class holds plain values).
//!
//! The class drives three OCCT leaf families that live in other packages:
//! `BSplSLib::PrepareEval` / `BSplSLib::BuildCache` (TKM/BSplSLib),
//! `BSplCLib::Bohm` / `BSplCLib::BuildKnots` / `BSplCLib::LocateParameter`
//! (TKM/BSplCLib) and `PLib::Trimming` / `UTrimming` / `VTrimming`
//! (TKM/PLib).  Their OCCT homes are rcad `math/bspl_lib.rs` and
//! `math/plib.rs`, which are outside this file's ownership; the leaf bodies
//! used by this class are therefore re-hosted in this file (its only
//! consumer) with their OCCT line anchors.  The already-translated
//! `BSplCLib::LocateParameter` (`math/bspl_lib.rs::locate_parameter_flat`),
//! the knot-sequence construction (`math/bspl_lib.rs::knot_sequence`) and
//! `BSplCLib::Hunt` (`math/bspl_lib.rs::hunt`) are CALLED, not duplicated.

use glam::DVec3;

use crate::base::proj_lib::proj_lib_projected_curve::IsoType;
use crate::geom::{BSplineSurface, Surface3, SurfaceEval};
use crate::math::bspl_lib::{hunt, knot_sequence, knot_sequence_length, locate_parameter_flat};
use crate::math::convert_grid_polynomial_to_poles::ConvertGridPolynomialToPoles;
use crate::math::direct_polynomial_roots::epsilon;

// =========================================================================
// The OCCT leaf helpers this class drives
// =========================================================================

/// OCCT BSplCLib::BuildKnots (BSplCLib.cxx L1555-1770) — the local
/// `2*Degree` knot window, `Mults == nullptr` arm: the unrolled OCCT switch
/// copies `Knots(Index - Degree + i)`, i = 0..2*Degree (the same body as the
/// private `math/bspl_lib.rs::build_knots_local`).
///
/// The `Mults != nullptr` arm is NOT reachable from Geom_OsculatingSurface:
/// `BSplSLib::BuildCache` is called with `BSplCLib::NoMults()`
/// (Geom_OsculatingSurface.cxx L646-661).
fn bspl_clib_build_knots(degree: i32, index: i32, knots: &[f64], out: &mut [f64]) {
    let mut j = index - degree;
    for i in 0..(2 * degree as usize) {
        j += 1;
        out[i] = knots[(j - 1) as usize];
    }
}

/// OCCT BSplCLib::Bohm (BSplCLib.cxx L1197-1541) — the in-place Taylor
/// (divided-difference) scheme over a flat pole array with element stride
/// `dimension`, on the 0-based local knot window `knot` of `2*Degree`
/// entries.
///
/// OCCT hoists the `Dimension` 1/2/3/4 bodies into cases 1-4; the `default`
/// body is the same algorithm for an arbitrary stride (the OCCT cases carry
/// the identical pointer arithmetic, the only rounding difference being the
/// explicit `0.0` of case 1 versus `x *= 0.0` in the generic body, which
/// cannot change any downstream sum).  This transcription keeps the generic
/// body for every stride.
pub(crate) fn bspl_clib_bohm(u: f64, degree: i32, n: i32, knot: &[f64], dimension: usize, poles: &mut [f64]) {
    let min = if n < degree { n } else { degree };
    let degm1 = degree - 1;
    let mut ddmi = (degree << 1) + 1;
    let dim = dimension as i32;
    let dim2 = dim << 1;

    // First phase, independent of U: the poles of the derivatives.
    let ps_dd = degree * dim;
    for i in 0..degree {
        ddmi -= 1;
        let mut pole = ps_dd;
        let mut tbis = ps_dd - dim;
        let mut jdmi = ddmi;
        let mut j = degm1;
        while j >= i {
            jdmi -= 1;
            let coef = if knot[jdmi as usize] == knot[j as usize] {
                0.0
            } else {
                1.0 / (knot[jdmi as usize] - knot[j as usize])
            };
            for _ in 0..dim {
                poles[pole as usize] -= poles[tbis as usize];
                poles[pole as usize] *= coef;
                pole += 1;
                tbis += 1;
            }
            pole -= dim2;
            tbis -= dim2;
            j -= 1;
        }
    }

    // Second phase, dependent on U.
    let mut i_dim = -dim;
    for i in 0..degree {
        i_dim += dim;
        let mut pole = i_dim;
        let mut tbis = pole + dim;
        let coef = u - knot[i as usize];
        let mut j = i;
        while j >= 0 {
            for _ in 0..dim {
                poles[pole as usize] += coef * poles[tbis as usize];
                pole += 1;
                tbis += 1;
            }
            pole -= dim2;
            tbis -= dim2;
            j -= 1;
        }
    }

    // Multiply by the degrees.
    let mut coef = degree as f64;
    let mut dmi = degree;
    let mut pole = dim;
    for _ in 1..=min {
        for _ in 0..dim {
            poles[pole as usize] *= coef;
            pole += 1;
        }
        dmi -= 1;
        coef *= dmi as f64;
    }
}

/// OCCT PLib::Trimming(U1, U2, dim, Coefs, WCoefs) (PLib.cxx L1632-1720) —
/// the polynomial re-parameterisation `u -> v = (u - U1) / (U2 - U1)` in the
/// Horner scheme over a flat coefficient array with element stride `dim`.
///
/// Architecture difference: OCCT reads/writes the 1-based
/// `NCollection_Array1<double> Coefs`; the rcad slice is 0-based, so every
/// OCCT index `k` appears here as `k - 1`.  The `WCoefs` arm is absent
/// because Geom_OsculatingSurface passes `PLib::NoWeights2()`
/// (Geom_OsculatingSurface.cxx L689/L728).
fn plib_trimming(u1: f64, u2: f64, dim: usize, coefs: &mut [f64]) {
    let lsp = u2 - u1;
    let dim_i = dim as i32;
    // OCCT: upc = Coefs.Upper() - dim + 1, len = Coefs.Length() / dim - 1.
    let upc = coefs.len() as i32 - dim_i + 1;
    let len = coefs.len() as i32 / dim_i - 1;

    for i in 1..=len {
        let mut indc = upc - dim_i * (i - 1);
        // Lowest-degree coefficient of iteration i.
        for j in 0..dim_i {
            let w = indc - dim_i + j - 1;
            let r = indc + j - 1;
            coefs[w as usize] += u1 * coefs[r as usize];
        }
        // Intermediate coefficients.
        while indc < upc {
            indc += dim_i;
            for k in 0..dim_i {
                let w = indc - dim_i + k - 1;
                let r = indc + k - 1;
                coefs[w as usize] = u1 * coefs[r as usize] + lsp * coefs[w as usize];
            }
        }
        // Highest-degree coefficient.
        for j in 0..dim_i {
            let w = upc + j - 1;
            coefs[w as usize] *= lsp;
        }
    }
}

/// OCCT PLib::UTrimming(U1, U2, Coeffs, WCoeffs) (PLib.cxx L1839-1883) — the
/// per-COLUMN 1D trim.  `PLib::SetPoles` / `GetPoles` (PLib.cxx L115-131,
/// L145-161) marshall the `gp_Pnt` array into the flat stride-3 double array;
/// the rcad grid is already a `DVec3` grid, so the same marshalling is the
/// x/y/z gather below.
fn plib_u_trimming(u1: f64, u2: f64, coeffs: &mut [Vec<DVec3>]) {
    let rows = coeffs.len();
    if rows == 0 {
        return;
    }
    let cols = coeffs[0].len();
    for icol in 0..cols {
        let mut temp = vec![0.0f64; 3 * rows];
        for irow in 0..rows {
            let p = coeffs[irow][icol];
            temp[3 * irow] = p.x;
            temp[3 * irow + 1] = p.y;
            temp[3 * irow + 2] = p.z;
        }
        plib_trimming(u1, u2, 3, &mut temp);
        for irow in 0..rows {
            coeffs[irow][icol] =
                DVec3::new(temp[3 * irow], temp[3 * irow + 1], temp[3 * irow + 2]);
        }
    }
}

/// OCCT PLib::VTrimming(V1, V2, Coeffs, WCoeffs) (PLib.cxx L1885-1934) — the
/// per-ROW 1D trim (same marshalling as [`plib_u_trimming`]).
fn plib_v_trimming(v1: f64, v2: f64, coeffs: &mut [Vec<DVec3>]) {
    for row in coeffs.iter_mut() {
        let cols = row.len();
        let mut temp = vec![0.0f64; 3 * cols];
        for icol in 0..cols {
            let p = row[icol];
            temp[3 * icol] = p.x;
            temp[3 * icol + 1] = p.y;
            temp[3 * icol + 2] = p.z;
        }
        plib_trimming(v1, v2, 3, &mut temp);
        for icol in 0..cols {
            row[icol] = DVec3::new(temp[3 * icol], temp[3 * icol + 1], temp[3 * icol + 2]);
        }
    }
}

/// The distinct knots of a flat (multiplicity-expanded) knot vector — the
/// rcad counterpart of the OCCT `Geom_BSplineSurface::UKnots()` /
/// `VKnots()` arrays.
fn distinct_knots(flat: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i == 0 || *k != out[out.len() - 1] {
            out.push(*k);
        }
    }
    out
}

/// OCCT `Geom_BSplineSurface::UKnot(Index)` / `VKnot(Index)`.
fn distinct_knot(flat: &[f64], index: i32) -> f64 {
    distinct_knots(flat)[(index - 1) as usize]
}

/// OCCT `Geom_BSplineSurface::NbUKnots()` / `NbVKnots()`.
fn distinct_knot_count(flat: &[f64]) -> i32 {
    distinct_knots(flat).len() as i32
}

/// OCCT `Geom_Surface::EvalD1(U, V)` over the rcad surface value — the union
/// of the translated per-type bodies: the ElSLib::DN arms for the quadrics,
/// the `Geom_BSplineSurface::EvalDN` / `Geom_BezierSurface::EvalDN`
/// (`BSplSLib::DN`) arms for the polynomial kinds
/// ([`crate::geom::offset_surface_utils::eval_dn`], i.e.
/// [`crate::geom::eval_b`]), and the swept / GeomEval / offset arms translated
/// with them ([`crate::geom::extrusion_utils`],
/// [`crate::geom::revolution_utils`], [`crate::geom::eval_c`],
/// [`crate::geom::offset_surface_utils`]).
pub(crate) fn surface_d1(the_s: &Surface3, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
    match the_s {
        // The exact per-type D1 of the rcad GeomAdaptor_Surface DN engine
        // (Plane/Cylinder/Cone/Sphere/Torus = ElSLib::DN).
        Surface3::Plane(_)
        | Surface3::Cylinder(_)
        | Surface3::Cone(_)
        | Surface3::Sphere(_)
        | Surface3::Torus(_) => (
            SurfaceEval::point_at(the_s, u, v),
            the_s.dn(u, v, 1, 0),
            the_s.dn(u, v, 0, 1),
        ),
        // Geom_BSplineSurface::EvalDN / Geom_BezierSurface::EvalDN
        // (BSplSLib::DN).
        Surface3::BSpline(_) | Surface3::Bezier(_) => (
            SurfaceEval::point_at(the_s, u, v),
            crate::geom::offset_surface_utils::eval_dn(the_s, u, v, 1, 0),
            crate::geom::offset_surface_utils::eval_dn(the_s, u, v, 0, 1),
        ),
        // Geom_SurfaceOfLinearExtrusion::EvalD1 / Geom_SurfaceOfRevolution::
        // EvalD1 / GeomEval_EllipsoidSurface::EvalD1 /
        // GeomEval_CircularHelicoidSurface::EvalD1 / Geom_OffsetSurface::EvalD1
        // (the leaves are called directly: `eval_d1` itself routes here).
        Surface3::LinearExtrusion(le) => {
            let d1 = crate::geom::extrusion_utils::linear_extrusion_eval_d1(le, u, v);
            (d1.point, d1.d1u, d1.d1v)
        }
        Surface3::Revolution(rev) => {
            let d1 = crate::geom::revolution_utils::revolution_eval_d1(rev, u, v);
            (d1.point, d1.d1u, d1.d1v)
        }
        Surface3::Ellipsoid(el) => {
            let d1 = crate::geom::eval_c::ellipsoid_eval_d1(el, u, v);
            (d1.point, d1.d1u, d1.d1v)
        }
        Surface3::Helicoid(h) => {
            let d1 = crate::geom::eval_c::helicoid_eval_d1(h, u, v);
            (d1.point, d1.d1u, d1.d1v)
        }
        Surface3::Offset(of) => {
            let d1 = crate::geom::offset_surface_utils::offset_payload_eval_d1(of, u, v);
            (d1.point, d1.d1u, d1.d1v)
        }
        _ => panic!(
            "GAP: Geom_Surface::EvalD1 (TKG3d/Geom) is not translated for this surface \
             type (the rcad GeomAdaptor_Surface DN engine covers the ElSLib surfaces, \
             Geom_BSplineSurface, Geom_BezierSurface, Geom_SurfaceOfLinearExtrusion, \
             Geom_SurfaceOfRevolution, GeomEval_EllipsoidSurface, \
             GeomEval_CircularHelicoidSurface and Geom_OffsetSurface) — \
             Geom_OsculatingSurface::isQPunctual"
        ),
    }
}

// =========================================================================
// Geom_OsculatingSurface
// =========================================================================

/// OCCT Geom_OsculatingSurface (Geom_OsculatingSurface.pxx L38-130).
#[derive(Clone, Debug)]
pub struct OsculatingSurface {
    /// OCCT: occ::handle<Geom_Surface> myBasisSurf (null before Init).
    my_basis_surf: Option<Surface3>,
    /// OCCT: double myTol.
    my_tol: f64,
    /// OCCT: NCollection_Sequence<handle(Geom_BSplineSurface)> myOsculSurf1.
    my_oscul_surf1: Vec<BSplineSurface>,
    /// OCCT: NCollection_Sequence<handle(Geom_BSplineSurface)> myOsculSurf2.
    my_oscul_surf2: Vec<BSplineSurface>,
    /// OCCT: NCollection_Sequence<int> myKdeg.
    my_kdeg: Vec<i32>,
    /// OCCT: std::array<bool, 4> myAlong.
    my_along: [bool; 4],
}

impl Default for OsculatingSurface {
    fn default() -> Self {
        OsculatingSurface::new()
    }
}

impl OsculatingSurface {
    /// OCCT Geom_OsculatingSurface() (cxx L32-36).
    pub fn new() -> Self {
        OsculatingSurface {
            my_basis_surf: None,
            my_tol: 0.0,
            my_oscul_surf1: Vec::new(),
            my_oscul_surf2: Vec::new(),
            my_kdeg: Vec::new(),
            my_along: [false, false, false, false],
        }
    }

    /// OCCT Geom_OsculatingSurface(theBS, theTol) (cxx L40-45).
    pub fn with_surface(the_bs: &Surface3, the_tol: f64) -> Self {
        let mut a = OsculatingSurface {
            my_basis_surf: None,
            my_tol: 0.0,
            my_oscul_surf1: Vec::new(),
            my_oscul_surf2: Vec::new(),
            my_kdeg: Vec::new(),
            my_along: [false, false, false, false],
        };
        a.init(the_bs, the_tol);
        a
    }

    /// OCCT Geom_OsculatingSurface::Init (cxx L105-391).
    ///
    /// The `unused_assignments` allowance covers the OCCT locals (`k`, `s`,
    /// `u_knot`, `v_knot`) which are re-initialised at the top of every
    /// `while (IsQPunc)` entry exactly as in the C++ body.
    #[allow(unused_assignments)]
    pub fn init(&mut self, the_bs: &Surface3, the_tol: f64) {
        self.clear_oscul_flags();
        self.my_tol = the_tol;
        // OCCT L109: consider all singularities below Tol, not just above
        // 1.e-12.
        let tol_min = 0.0;
        let mut oscul_surf = true;
        self.my_basis_surf = Some(the_bs.clone());
        self.my_oscul_surf1.clear();
        self.my_oscul_surf2.clear();
        self.my_kdeg.clear();

        if !matches!(the_bs, Surface3::BSpline(_) | Surface3::Bezier(_)) {
            self.clear_oscul_flags();
            return;
        }

        let bounds = SurfaceEval::default_domain(the_bs);
        let (u1, u2, v1, v2) = (bounds[0], bounds[1], bounds[2], bounds[3]);

        self.my_along[0] = self.is_q_punctual(the_bs, v1, IsoType::IsoV, tol_min, the_tol);
        self.my_along[1] = self.is_q_punctual(the_bs, v2, IsoType::IsoV, tol_min, the_tol);
        self.my_along[2] = self.is_q_punctual(the_bs, u1, IsoType::IsoU, tol_min, the_tol);
        self.my_along[3] = self.is_q_punctual(the_bs, u2, IsoType::IsoU, tol_min, the_tol);

        if !(self.my_along[0] || self.my_along[1] || self.my_along[2] || self.my_along[3]) {
            return;
        }

        // OCCT: occ::handle<Geom_BSplineSurface> InitSurf, L, S.
        let init_surf: BSplineSurface = match the_bs {
            Surface3::Bezier(bzs) => {
                let u_degree = bzs.control_points.len() - 1;
                let v_degree = bzs.control_points.first().map(|r| r.len()).unwrap_or(1) - 1;
                let mut knots_u = vec![0.0f64; u_degree + 1];
                knots_u.extend(std::iter::repeat(1.0).take(u_degree + 1));
                let mut knots_v = vec![0.0f64; v_degree + 1];
                knots_v.extend(std::iter::repeat(1.0).take(v_degree + 1));
                BSplineSurface {
                    degree_u: u_degree,
                    degree_v: v_degree,
                    knots_u,
                    knots_v,
                    control_points: bzs.control_points.clone(),
                    weights: bzs.weights.clone(),
                    // OCCT Geom_BezierSurface::IsUPeriodic / IsVPeriodic
                    // return false (Geom_BezierSurface.cxx L1977-1986).
                    is_periodic_u: false,
                    is_periodic_v: false,
                }
            }
            Surface3::BSpline(b) => b.clone(),
            _ => unreachable!(),
        };

        let mut l: Option<BSplineSurface> = None;
        let mut s: BSplineSurface = init_surf.clone();
        let mut k = 0i32;
        let mut u_knot = 1i32;
        let mut v_knot = 1i32;
        let nb_u_knots = distinct_knot_count(&init_surf.knots_u);
        let nb_v_knots = distinct_knot_count(&init_surf.knots_v);

        if self.is_along_u() && self.is_along_v() {
            self.clear_oscul_flags();
        }

        if (self.is_along_u() && init_surf.degree_v > 1)
            || (self.is_along_v() && init_surf.degree_u > 1)
        {
            if self.my_along[0] || self.my_along[1] {
                for i in 1..nb_u_knots {
                    if self.my_along[0] {
                        s = init_surf.clone();
                        k = 0;
                        let mut is_q_punc = true;
                        u_knot = i;
                        v_knot = 1;
                        while is_q_punc {
                            oscul_surf = match self.build_osculating_surface(v1, u_knot, v_knot, &s)
                            {
                                Some(built) => {
                                    l = Some(built);
                                    true
                                }
                                None => false,
                            };
                            if !oscul_surf {
                                break;
                            }
                            k += 1;
                            is_q_punc = self.is_q_punctual(
                                &Surface3::BSpline(l.clone().unwrap()),
                                v1,
                                IsoType::IsoV,
                                0.0,
                                the_tol,
                            );
                            u_knot = 1;
                            v_knot = 1;
                            s = l.clone().unwrap();
                        }
                        if oscul_surf {
                            self.my_oscul_surf1.push(l.clone().unwrap());
                        } else {
                            self.clear_oscul_flags();
                        }
                        if self.my_along[1] && oscul_surf {
                            s = init_surf.clone();
                            k = 0;
                            let mut is_q_punc = true;
                            u_knot = i;
                            v_knot = nb_v_knots - 1;
                            while is_q_punc {
                                oscul_surf =
                                    match self.build_osculating_surface(v2, u_knot, v_knot, &s) {
                                        Some(built) => {
                                            l = Some(built);
                                            true
                                        }
                                        None => false,
                                    };
                                if !oscul_surf {
                                    break;
                                }
                                k += 1;
                                is_q_punc = self.is_q_punctual(
                                    &Surface3::BSpline(l.clone().unwrap()),
                                    v2,
                                    IsoType::IsoV,
                                    0.0,
                                    the_tol,
                                );
                                u_knot = 1;
                                v_knot = 1;
                                s = l.clone().unwrap();
                            }
                            if oscul_surf {
                                self.my_oscul_surf2.push(l.clone().unwrap());
                                self.my_kdeg.push(k);
                            }
                        }
                    } else {
                        s = init_surf.clone();
                        k = 0;
                        let mut is_q_punc = true;
                        u_knot = i;
                        v_knot = nb_v_knots - 1;
                        while is_q_punc {
                            oscul_surf =
                                match self.build_osculating_surface(v2, u_knot, v_knot, &s) {
                                    Some(built) => {
                                        l = Some(built);
                                        true
                                    }
                                    None => false,
                                };
                            if !oscul_surf {
                                break;
                            }
                            k += 1;
                            is_q_punc = self.is_q_punctual(
                                &Surface3::BSpline(l.clone().unwrap()),
                                v2,
                                IsoType::IsoV,
                                0.0,
                                the_tol,
                            );
                            u_knot = 1;
                            v_knot = 1;
                            s = l.clone().unwrap();
                        }
                        if oscul_surf {
                            self.my_oscul_surf2.push(l.clone().unwrap());
                            self.my_kdeg.push(k);
                        } else {
                            self.clear_oscul_flags();
                        }
                    }
                }
            }

            if self.my_along[2] || self.my_along[3] {
                for i in 1..nb_v_knots {
                    if self.my_along[2] {
                        s = init_surf.clone();
                        k = 0;
                        let mut is_q_punc = true;
                        u_knot = 1;
                        v_knot = i;
                        while is_q_punc {
                            oscul_surf = match self.build_osculating_surface(u1, u_knot, v_knot, &s)
                            {
                                Some(built) => {
                                    l = Some(built);
                                    true
                                }
                                None => false,
                            };
                            if !oscul_surf {
                                break;
                            }
                            k += 1;
                            is_q_punc = self.is_q_punctual(
                                &Surface3::BSpline(l.clone().unwrap()),
                                u1,
                                IsoType::IsoU,
                                0.0,
                                the_tol,
                            );
                            u_knot = 1;
                            v_knot = 1;
                            s = l.clone().unwrap();
                        }
                        if oscul_surf {
                            self.my_oscul_surf1.push(l.clone().unwrap());
                        } else {
                            self.clear_oscul_flags();
                        }
                        if self.my_along[3] && oscul_surf {
                            s = init_surf.clone();
                            k = 0;
                            let mut is_q_punc = true;
                            u_knot = nb_u_knots - 1;
                            v_knot = i;
                            while is_q_punc {
                                oscul_surf =
                                    match self.build_osculating_surface(u2, u_knot, v_knot, &s) {
                                        Some(built) => {
                                            l = Some(built);
                                            true
                                        }
                                        None => false,
                                    };
                                if !oscul_surf {
                                    break;
                                }
                                k += 1;
                                is_q_punc = self.is_q_punctual(
                                    &Surface3::BSpline(l.clone().unwrap()),
                                    u2,
                                    IsoType::IsoU,
                                    0.0,
                                    the_tol,
                                );
                                u_knot = 1;
                                v_knot = 1;
                                s = l.clone().unwrap();
                            }
                            if oscul_surf {
                                self.my_oscul_surf2.push(l.clone().unwrap());
                                self.my_kdeg.push(k);
                            }
                        }
                    } else {
                        s = init_surf.clone();
                        k = 0;
                        let mut is_q_punc = true;
                        u_knot = nb_u_knots - 1;
                        v_knot = i;
                        while is_q_punc {
                            oscul_surf = match self.build_osculating_surface(u2, u_knot, v_knot, &s)
                            {
                                Some(built) => {
                                    l = Some(built);
                                    true
                                }
                                None => false,
                            };
                            if !oscul_surf {
                                break;
                            }
                            k += 1;
                            is_q_punc = self.is_q_punctual(
                                &Surface3::BSpline(l.clone().unwrap()),
                                u2,
                                IsoType::IsoU,
                                0.0,
                                the_tol,
                            );
                            u_knot = 1;
                            v_knot = 1;
                            s = l.clone().unwrap();
                        }
                        if oscul_surf {
                            self.my_oscul_surf2.push(l.clone().unwrap());
                            self.my_kdeg.push(k);
                        } else {
                            self.clear_oscul_flags();
                        }
                    }
                }
            }
        } else {
            self.clear_oscul_flags();
        }
    }

    /// OCCT Geom_OsculatingSurface::BasisSurface (pxx L69).
    pub fn basis_surface(&self) -> Option<&Surface3> {
        self.my_basis_surf.as_ref()
    }

    /// OCCT Geom_OsculatingSurface::Tolerance (pxx L72).
    pub fn tolerance(&self) -> f64 {
        self.my_tol
    }

    /// OCCT Geom_OsculatingSurface::HasOscSurf (pxx L75).
    pub fn has_osc_surf(&self) -> bool {
        self.my_along.iter().any(|a| *a)
    }

    /// OCCT Geom_OsculatingSurface::IsAlongU (pxx L78).
    pub fn is_along_u(&self) -> bool {
        self.my_along[0] || self.my_along[1]
    }

    /// OCCT Geom_OsculatingSurface::IsAlongV (pxx L81).
    pub fn is_along_v(&self) -> bool {
        self.my_along[2] || self.my_along[3]
    }

    /// OCCT Geom_OsculatingSurface::UOsculatingSurface (cxx L395-461) —
    /// returns `(along, theT, theL)`: `along == false` is the OCCT `false`
    /// return, and `theL` is only meaningful when `along` is true.
    #[allow(unused_assignments)]
    pub fn u_osculating_surface(
        &self,
        the_u: f64,
        the_v: f64,
    ) -> (bool, bool, Option<BSplineSurface>) {
        let mut along = false;
        let mut the_t = false;
        let mut the_l: Option<BSplineSurface> = None;
        if self.my_along[0] || self.my_along[1] {
            let mut nu = 1i32;
            let mut nv = 1i32;
            // OCCT L406: myBasisSurf->Bounds(u1, u2, v1, v2) — the values are
            // not consumed by this body.
            let _bounds = SurfaceEval::default_domain(
                self.my_basis_surf
                    .as_ref()
                    .expect("Geom_OsculatingSurface: null myBasisSurf"),
            );
            let mut nb_uk = 1i32;
            let mut nb_vk = 2i32;
            let mut is_to_skip_second = false;
            if let Surface3::BSpline(bsur) = self
                .my_basis_surf
                .as_ref()
                .expect("Geom_OsculatingSurface: null myBasisSurf")
            {
                let u_knots = distinct_knots(&bsur.knots_u);
                let v_knots = distinct_knots(&bsur.knots_v);
                nb_uk = u_knots.len() as i32;
                nb_vk = v_knots.len() as i32;
                hunt(&u_knots, the_u, &mut nu);
                hunt(&v_knots, the_v, &mut nv);
                if nu < 1 {
                    nu = 1;
                }
                if nu >= nb_uk {
                    nu = nb_uk - 1;
                }
                if nb_vk == 2 && nv == 1 {
                    // Need to find the closest end.
                    if v_knots[(nb_vk - 1) as usize] - the_v > the_v - v_knots[0] {
                        is_to_skip_second = true;
                    }
                }
            } else {
                nu = 1;
                nv = 1;
                nb_vk = 2;
            }

            if self.my_along[0] && nv == 1 {
                the_l = Some(self.my_oscul_surf1[(nu - 1) as usize].clone());
                along = true;
            }
            if self.my_along[1] && (nv == nb_vk - 1) && !is_to_skip_second {
                // theT means that derivative vector of osculating surface is
                // opposite to the original. This happens when (v-t)^k is
                // negative, i.e. difference between degrees (k) is odd and t
                // is the last parameter.
                if self.my_kdeg[(nu - 1) as usize] % 2 != 0 {
                    the_t = true;
                }
                the_l = Some(self.my_oscul_surf2[(nu - 1) as usize].clone());
                along = true;
            }
        }
        (along, the_t, the_l)
    }

    /// OCCT Geom_OsculatingSurface::VOsculatingSurface (cxx L465-528).
    #[allow(unused_assignments)]
    pub fn v_osculating_surface(
        &self,
        the_u: f64,
        the_v: f64,
    ) -> (bool, bool, Option<BSplineSurface>) {
        let mut along = false;
        let mut the_t = false;
        let mut the_l: Option<BSplineSurface> = None;
        if self.my_along[2] || self.my_along[3] {
            let mut nu = 1i32;
            let mut nv = 1i32;
            // OCCT L476: myBasisSurf->Bounds(u1, u2, v1, v2) — not consumed.
            let _bounds = SurfaceEval::default_domain(
                self.my_basis_surf
                    .as_ref()
                    .expect("Geom_OsculatingSurface: null myBasisSurf"),
            );
            let mut nb_uk = 2i32;
            let mut nb_vk = 1i32;
            let mut is_to_skip_second = false;
            if let Surface3::BSpline(bsur) = self
                .my_basis_surf
                .as_ref()
                .expect("Geom_OsculatingSurface: null myBasisSurf")
            {
                let u_knots = distinct_knots(&bsur.knots_u);
                let v_knots = distinct_knots(&bsur.knots_v);
                nb_uk = u_knots.len() as i32;
                nb_vk = v_knots.len() as i32;
                hunt(&u_knots, the_u, &mut nu);
                hunt(&v_knots, the_v, &mut nv);
                if nv < 1 {
                    nv = 1;
                }
                if nv >= nb_vk {
                    nv = nb_vk - 1;
                }
                if nb_uk == 2 && nu == 1 {
                    // Need to find the closest end.
                    if u_knots[(nb_uk - 1) as usize] - the_u > the_u - u_knots[0] {
                        is_to_skip_second = true;
                    }
                }
            } else {
                nu = 1;
                nv = 1;
                nb_uk = 2;
            }

            if self.my_along[2] && nu == 1 {
                the_l = Some(self.my_oscul_surf1[(nv - 1) as usize].clone());
                along = true;
            }
            if self.my_along[3] && (nu == nb_uk - 1) && !is_to_skip_second {
                if self.my_kdeg[(nv - 1) as usize] % 2 != 0 {
                    the_t = true;
                }
                the_l = Some(self.my_oscul_surf2[(nv - 1) as usize].clone());
                along = true;
            }
        }
        (along, the_t, the_l)
    }

    /// OCCT Geom_OsculatingSurface::buildOsculatingSurface (cxx L532-780) —
    /// `None` is the OCCT `false` return ("the osculating surface can't be
    /// built"); the OCCT output handle is only written on success.
    fn build_osculating_surface(
        &self,
        the_param: f64,
        the_su_knot: i32,
        the_sv_knot: i32,
        the_bs: &BSplineSurface,
    ) -> Option<BSplineSurface> {
        let udeg = the_bs.degree_u as i32;
        let vdeg = the_bs.degree_v as i32;
        // OCCT L549-555: "surface osculatrice nulle".
        if (self.is_along_u() && vdeg <= 1) || (self.is_along_v() && udeg <= 1) {
            return None;
        }

        let min_degree = udeg.min(vdeg);
        let max_degree = udeg.max(vdeg);

        // NCollection_Array2<gp_Pnt> cachepoles(1, MaxDegree + 1, 1,
        // MinDegree + 1).
        let mut cache_poles =
            vec![vec![DVec3::ZERO; (min_degree + 1) as usize]; (max_degree + 1) as usize];

        // For polynomial grid.
        let mut max_u_degree = udeg;
        let mut max_v_degree = vdeg;
        // NCollection_HArray2<int> NumCoeffPerSurface(1, 1, 1, 2).
        let mut num_coeff_per_surface = [0i32; 2];
        // PolynomialUIntervals / PolynomialVIntervals / TrueUIntervals /
        // TrueVIntervals: NCollection_HArray1<double>(1, 2).
        let mut polynomial_u_intervals = [0.0f64; 2];
        let mut polynomial_v_intervals = [0.0f64; 2];
        let mut true_u_intervals = [0.0f64; 2];
        let mut true_v_intervals = [0.0f64; 2];

        for i in 1..=2 {
            polynomial_u_intervals[(i - 1) as usize] = (i - 1) as f64;
            polynomial_v_intervals[(i - 1) as usize] = (i - 1) as f64;
            true_u_intervals[(i - 1) as usize] = distinct_knot(&the_bs.knots_u, the_su_knot + i - 1);
            true_v_intervals[(i - 1) as usize] = distinct_knot(&the_bs.knots_v, the_sv_knot + i - 1);
        }

        let mut osc_u_num_coeff = 0i32;
        let mut osc_v_num_coeff = 0i32;
        if self.is_along_u() {
            osc_u_num_coeff = udeg + 1;
            osc_v_num_coeff = vdeg;
        }
        if self.is_along_v() {
            osc_u_num_coeff = udeg;
            osc_v_num_coeff = vdeg + 1;
        }
        num_coeff_per_surface[0] = osc_u_num_coeff;
        num_coeff_per_surface[1] = osc_v_num_coeff;
        let nbc = num_coeff_per_surface[0] * num_coeff_per_surface[1] * 3;
        if nbc == 0 {
            return None;
        }
        let mut coefficients = vec![0.0f64; nbc as usize];

        // Building the cache.
        let vcacheparameter0 = distinct_knot(&the_bs.knots_v, the_sv_knot);
        let ucacheparameter0 = distinct_knot(&the_bs.knots_u, the_su_knot);
        let vspanlength = distinct_knot(&the_bs.knots_v, the_sv_knot + 1) - vcacheparameter0;
        let uspanlength = distinct_knot(&the_bs.knots_u, the_su_knot + 1) - ucacheparameter0;

        // Always reduce to a parametrization such that locally it is the iso
        // u=0 or v=0 that is degenerate.
        let is_v_negative = the_param > vcacheparameter0 + vspanlength / 2.0;
        let is_u_negative = the_param > ucacheparameter0 + uspanlength / 2.0;

        let mut vcacheparameter = vcacheparameter0;
        let mut ucacheparameter = ucacheparameter0;
        if self.is_along_u() && (the_param > vcacheparameter + vspanlength / 2.0) {
            vcacheparameter += vspanlength;
        }
        if self.is_along_v() && (the_param > ucacheparameter + uspanlength / 2.0) {
            ucacheparameter += uspanlength;
        }

        // OCCT L646-661: BSplSLib::BuildCache(ucacheparameter,
        // vcacheparameter, uspanlength, vspanlength, IsUPeriodic(),
        // IsVPeriodic(), UDegree(), VDegree(), ULocalIndex, VLocalIndex,
        // aUFlatKnts, aVFlatKnts, aPoles, NoWeights(), cachepoles,
        // NoWeights()) — ULocalIndex/VLocalIndex are 0 on entry.
        bspl_slib_build_cache(
            ucacheparameter,
            vcacheparameter,
            uspanlength,
            vspanlength,
            the_bs.is_periodic_u,
            the_bs.is_periodic_v,
            udeg,
            vdeg,
            0,
            0,
            &the_bs.knots_u,
            &the_bs.knots_v,
            &the_bs.control_points,
            None,
            None,
            &mut cache_poles,
        );

        let mut osc_coeff = vec![
            vec![DVec3::ZERO; osc_v_num_coeff.max(0) as usize];
            osc_u_num_coeff.max(0) as usize
        ];
        let mut index;

        if self.is_along_u() {
            if udeg > vdeg {
                for n in 1..=udeg + 1 {
                    for m in 1..=vdeg {
                        osc_coeff[(n - 1) as usize][(m - 1) as usize] =
                            cache_poles[(n - 1) as usize][m as usize];
                    }
                }
            } else {
                for n in 1..=udeg + 1 {
                    for m in 1..=vdeg {
                        osc_coeff[(n - 1) as usize][(m - 1) as usize] =
                            cache_poles[m as usize][(n - 1) as usize];
                    }
                }
            }
            if is_v_negative {
                plib_v_trimming(-1.0, 0.0, &mut osc_coeff);
            }

            index = 1;
            for n in 1..=udeg + 1 {
                for m in 1..=vdeg {
                    let p = osc_coeff[(n - 1) as usize][(m - 1) as usize];
                    coefficients[(index - 1) as usize] = p.x;
                    index += 1;
                    coefficients[(index - 1) as usize] = p.y;
                    index += 1;
                    coefficients[(index - 1) as usize] = p.z;
                    index += 1;
                }
            }
        }

        if self.is_along_v() {
            if udeg > vdeg {
                for n in 1..=udeg {
                    for m in 1..=vdeg + 1 {
                        osc_coeff[(n - 1) as usize][(m - 1) as usize] =
                            cache_poles[n as usize][(m - 1) as usize];
                    }
                }
            } else {
                for n in 1..=udeg {
                    for m in 1..=vdeg + 1 {
                        osc_coeff[(n - 1) as usize][(m - 1) as usize] =
                            cache_poles[(m - 1) as usize][n as usize];
                    }
                }
            }
            if is_u_negative {
                plib_u_trimming(-1.0, 0.0, &mut osc_coeff);
            }
            index = 1;
            for n in 1..=udeg {
                for m in 1..=vdeg + 1 {
                    let p = osc_coeff[(n - 1) as usize][(m - 1) as usize];
                    coefficients[(index - 1) as usize] = p.x;
                    index += 1;
                    coefficients[(index - 1) as usize] = p.y;
                    index += 1;
                    coefficients[(index - 1) as usize] = p.z;
                    index += 1;
                }
            }
        }

        if self.is_along_u() {
            max_v_degree -= 1;
        }
        if self.is_along_v() {
            max_u_degree -= 1;
        }
        let u_continuity = -1;
        let v_continuity = -1;

        let data = ConvertGridPolynomialToPoles::from_grid(
            1,
            1,
            u_continuity,
            v_continuity,
            max_u_degree,
            max_v_degree,
            &num_coeff_per_surface,
            &coefficients,
            &polynomial_u_intervals,
            &polynomial_v_intervals,
            &true_u_intervals,
            &true_v_intervals,
        );

        // theBSpl = new Geom_BSplineSurface(Data.Poles(), Data.UKnots(),
        // Data.VKnots(), Data.UMultiplicities(), Data.VMultiplicities(),
        // Data.UDegree(), Data.VDegree(), false, false).
        let u_degree = data.u_degree().max(0) as usize;
        let v_degree = data.v_degree().max(0) as usize;
        let nu = knot_sequence_length(data.u_multiplicities(), u_degree, false);
        let nv = knot_sequence_length(data.v_multiplicities(), v_degree, false);
        let mut knots_u = vec![0.0f64; nu];
        knot_sequence(
            data.u_knots(),
            data.u_multiplicities(),
            u_degree,
            false,
            &mut knots_u,
        );
        let mut knots_v = vec![0.0f64; nv];
        knot_sequence(
            data.v_knots(),
            data.v_multiplicities(),
            v_degree,
            false,
            &mut knots_v,
        );
        let nb_u = data.nb_u_poles();
        let nb_v = data.nb_v_poles();
        Some(BSplineSurface {
            degree_u: u_degree,
            degree_v: v_degree,
            knots_u,
            knots_v,
            control_points: data.poles().clone(),
            weights: vec![vec![1.0f64; nb_v]; nb_u],
            is_periodic_u: false,
            is_periodic_v: false,
        })
    }

    /// OCCT Geom_OsculatingSurface::isQPunctual (cxx L784-832) — true when the
    /// isoparametric `the_it` at `the_param` keeps a `D1` norm inside
    /// [the_tol_min, the_tol_max] over the whole opposite domain.
    fn is_q_punctual(
        &self,
        the_s: &Surface3,
        the_param: f64,
        the_it: IsoType,
        the_tol_min: f64,
        the_tol_max: f64,
    ) -> bool {
        let mut along = true;
        let bounds = SurfaceEval::default_domain(the_s);
        let (u1, u2, v1, v2) = (bounds[0], bounds[1], bounds[2], bounds[3]);
        let d1_norm_max;
        if the_it == IsoType::IsoV {
            let step = (u2 - u1) / 10.0;
            let mut m = 0.0f64;
            let mut t = u1;
            while t <= u2 {
                let (_p, d1u, _d1v) = surface_d1(the_s, t, the_param);
                m = m.max(d1u.length());
                t += step;
            }
            d1_norm_max = m;
        } else {
            let step = (v2 - v1) / 10.0;
            let mut m = 0.0f64;
            let mut t = v1;
            while t <= v2 {
                let (_p, _d1u, d1v) = surface_d1(the_s, the_param, t);
                m = m.max(d1v.length());
                t += step;
            }
            d1_norm_max = m;
        }
        if d1_norm_max > the_tol_max || d1_norm_max < the_tol_min {
            along = false;
        }
        along
    }

    /// OCCT Geom_OsculatingSurface::clearOsculFlags (cxx L836-839).
    fn clear_oscul_flags(&mut self) {
        self.my_along = [false, false, false, false];
    }
}

// =========================================================================
// BSplSLib::PrepareEval / BSplSLib::BuildCache
// =========================================================================

/// OCCT `BSplSLib_DataContainer` (BSplSLib.hxx) — the scratch buffers of
/// `PrepareEval`, re-hosted as owned vectors (`knots1` = `dc.knots1`,
/// `knots2` = `dc.knots2`, `poles` = `dc.poles`).
pub(crate) struct PrepareEvalData {
    /// The local `2*d1` knot window.
    pub(crate) knots1: Vec<f64>,
    /// The local `2*d2` knot window.
    pub(crate) knots2: Vec<f64>,
    /// The `(d1+1) * (d2+1)` local poles, element stride
    /// `rational ? 4 : 3`.
    pub(crate) poles: Vec<f64>,
}

/// The "verify if locally non rational" block of OCCT PrepareEval
/// (BSplSLib.cxx L477-520 for the UDegree <= VDegree arm, L613-655 for the
/// other): it re-reads the (UDegree+1) x (VDegree+1) local weight block and
/// clears `rational` when every weight equals the first one within
/// `Epsilon(w)`.
#[allow(clippy::too_many_arguments)]
fn check_locally_rational(
    rational: &mut bool,
    u_degree: i32,
    v_degree: i32,
    p_lower_row: i32,
    p_upper_row: i32,
    p_lower_col: i32,
    p_upper_col: i32,
    weights: &[Vec<f64>],
    ip0: i32,
    jp0: i32,
) {
    if !*rational {
        return;
    }
    *rational = false;
    let mut ip = ip0;
    let mut jp = jp0;
    if ip < p_lower_row {
        ip = p_upper_row;
    }
    if jp < p_lower_col {
        jp = p_upper_col;
    }
    let w = weights[(ip - p_lower_row) as usize][(jp - p_lower_col) as usize];
    let eps = epsilon(w);
    for _ in 0..=u_degree {
        if *rational {
            break;
        }
        jp = jp0;
        if jp < p_lower_col {
            jp = p_upper_col;
        }
        for _ in 0..=v_degree {
            let mut dw = weights[(ip - p_lower_row) as usize][(jp - p_lower_col) as usize] - w;
            if dw < 0.0 {
                dw = -dw;
            }
            *rational = dw > eps;
            if *rational {
                break;
            }
            jp += 1;
            if jp > p_upper_col {
                jp = p_lower_col;
            }
        }
        ip += 1;
        if ip > p_upper_row {
            ip = p_lower_row;
        }
    }
}

/// The out-parameters of `PrepareEval` (OCCT passes them by reference) plus
/// the scratch container.  `flag_u_or_v` is the OCCT return value.
pub(crate) struct PrepareEvalResult {
    pub(crate) flag_u_or_v: bool,
    pub(crate) u1: f64,
    pub(crate) u2: f64,
    pub(crate) d1: i32,
    pub(crate) d2: i32,
    pub(crate) rational: bool,
    pub(crate) dc: PrepareEvalData,
}

/// OCCT BSplSLib::PrepareEval (BSplSLib.cxx L313-728) — builds the local knot
/// windows and the local pole block of the current span.
///
/// `u_knots` / `v_knots` are the flat knot sequences (`UKnotSequence()` /
/// `VKnotSequence()`); the OCCT call convention `Mults ==
/// BSplCLib::NoMults()` is assumed for the `Mults` arguments, so both
/// branches take the OCCT `nullptr` arm (`uindex -= UKLower + UDegree`) and
/// `BSplCLib::LocateParameter` (the `Mults == nullptr` form) is delegated to
/// `math/bspl_lib.rs::locate_parameter_flat`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn prepare_eval(
    u: f64,
    v: f64,
    u_index: i32,
    v_index: i32,
    u_degree: i32,
    v_degree: i32,
    u_rat: bool,
    v_rat: bool,
    u_per: bool,
    v_per: bool,
    poles: &[Vec<DVec3>],
    weights: Option<&[Vec<f64>]>,
    u_knots: &[f64],
    v_knots: &[f64],
) -> PrepareEvalResult {
    let mut rational = u_rat || v_rat;
    let u_k_lower = 1i32;
    let u_k_upper = u_knots.len() as i32;
    let v_k_lower = 1i32;
    let v_k_upper = v_knots.len() as i32;

    let p_lower_row = 1i32;
    let p_upper_row = poles.len() as i32;
    let p_lower_col = 1i32;
    let p_upper_col = poles.first().map(|r| r.len()).unwrap_or(0) as i32;

    // OCCT `Poles.Value(i, j)` / `Weights->Value(i, j)` on the 1-based
    // (row = U index, col = V index) pole grid.
    let pole_at = |i: i32, j: i32| poles[(i - p_lower_row) as usize][(j - p_lower_col) as usize];
    let weight_at = |i: i32, j: i32| -> f64 {
        match weights {
            Some(w) => w[(i - p_lower_row) as usize][(j - p_lower_col) as usize],
            None => 1.0,
        }
    };
    let _ = &weight_at;
    if u_degree <= v_degree {
        // Compute the indices.
        let mut uindex = u_index;
        let mut vindex = v_index;
        let u1;
        let u2;
        if uindex < u_k_lower || uindex > u_k_upper {
            let mut new_u = u;
            locate_parameter_flat(
                u_degree as usize,
                u_knots,
                u,
                u_per,
                u_k_lower + u_degree,
                u_k_upper - u_degree,
                &mut uindex,
                &mut new_u,
            );
            u1 = new_u;
        } else {
            u1 = u;
        }
        if vindex < v_k_lower || vindex > v_k_upper {
            let mut new_u = v;
            locate_parameter_flat(
                v_degree as usize,
                v_knots,
                v,
                v_per,
                v_k_lower + v_degree,
                v_k_upper - v_degree,
                &mut vindex,
                &mut new_u,
            );
            u2 = new_u;
        } else {
            u2 = v;
        }

        // Get the knots.
        let d1 = u_degree;
        let d2 = v_degree;
        let mut knots1 = vec![0.0f64; 2 * d1 as usize];
        let mut knots2 = vec![0.0f64; 2 * d2 as usize];
        bspl_clib_build_knots(d1, uindex, u_knots, &mut knots1);
        bspl_clib_build_knots(d2, vindex, v_knots, &mut knots2);

        // OCCT: UMults == nullptr -> uindex -= UKLower + UDegree.
        uindex -= u_k_lower + u_degree;
        vindex -= v_k_lower + v_degree;

        // Verify if locally non rational.
        if let Some(w) = weights {
            check_locally_rational(
                &mut rational,
                u_degree,
                v_degree,
                p_lower_row,
                p_upper_row,
                p_lower_col,
                p_upper_col,
                w,
                p_lower_row + uindex,
                p_lower_col + vindex,
            );
        }

        // Copy the poles (U is the major index).
        let dim: usize = if rational { 4 } else { 3 };
        let mut flat = vec![0.0f64; ((d1 + 1) * (d2 + 1)) as usize * dim];
        let mut k = 0usize;
        let mut ip = p_lower_row + uindex;
        if ip < p_lower_row {
            ip = p_upper_row;
        }
        for _ in 0..=d1 {
            let mut jp = p_lower_col + vindex;
            if jp < p_lower_col {
                jp = p_upper_col;
            }
            for _ in 0..=d2 {
                let p = pole_at(ip, jp);
                if rational {
                    let w = weight_at(ip, jp);
                    flat[k] = p.x * w;
                    flat[k + 1] = p.y * w;
                    flat[k + 2] = p.z * w;
                    flat[k + 3] = w;
                } else {
                    flat[k] = p.x;
                    flat[k + 1] = p.y;
                    flat[k + 2] = p.z;
                }
                k += dim;
                jp += 1;
                if jp > p_upper_col {
                    jp = p_lower_col;
                }
            }
            ip += 1;
            if ip > p_upper_row {
                ip = p_lower_row;
            }
        }

        return PrepareEvalResult {
            flag_u_or_v: true,
            u1,
            u2,
            d1,
            d2,
            rational,
            dc: PrepareEvalData {
                knots1,
                knots2,
                poles: flat,
            },
        };
    }

    // UDegree > VDegree: the first direction computed is V.
    let mut uindex = u_index;
    let mut vindex = v_index;
    let u1;
    let u2;
    if uindex < u_k_lower || uindex > u_k_upper {
        let mut new_u = u;
        locate_parameter_flat(
            u_degree as usize,
            u_knots,
            u,
            u_per,
            u_k_lower + u_degree,
            u_k_upper - u_degree,
            &mut uindex,
            &mut new_u,
        );
        u2 = new_u;
    } else {
        u2 = u;
    }
    if vindex < v_k_lower || vindex > v_k_upper {
        let mut new_u = v;
        locate_parameter_flat(
            v_degree as usize,
            v_knots,
            v,
            v_per,
            v_k_lower + v_degree,
            v_k_upper - v_degree,
            &mut vindex,
            &mut new_u,
        );
        u1 = new_u;
    } else {
        u1 = v;
    }

    // Get the knots.
    let d2 = u_degree;
    let d1 = v_degree;
    let mut knots1 = vec![0.0f64; 2 * d1 as usize];
    let mut knots2 = vec![0.0f64; 2 * d2 as usize];
    bspl_clib_build_knots(d2, uindex, u_knots, &mut knots2);
    bspl_clib_build_knots(d1, vindex, v_knots, &mut knots1);

    // OCCT: UMults == nullptr -> uindex -= UKLower + UDegree.
    uindex -= u_k_lower + u_degree;
    vindex -= v_k_lower + v_degree;

    // Verify if locally non rational.
    if let Some(w) = weights {
        check_locally_rational(
            &mut rational,
            u_degree,
            v_degree,
            p_lower_row,
            p_upper_row,
            p_lower_col,
            p_upper_col,
            w,
            p_lower_row + uindex,
            p_lower_col + vindex,
        );
    }

    // Copy the poles (V is the major index).
    let dim: usize = if rational { 4 } else { 3 };
    let mut flat = vec![0.0f64; ((d1 + 1) * (d2 + 1)) as usize * dim];
    let mut k = 0usize;
    let mut jp = p_lower_col + vindex;
    if jp < p_lower_col {
        jp = p_upper_col;
    }
    for _ in 0..=d1 {
        let mut ip = p_lower_row + uindex;
        if ip < p_lower_row {
            ip = p_upper_row;
        }
        for _ in 0..=d2 {
            let p = pole_at(ip, jp);
            if rational {
                let w = weight_at(ip, jp);
                flat[k] = p.x * w;
                flat[k + 1] = p.y * w;
                flat[k + 2] = p.z * w;
                flat[k + 3] = w;
            } else {
                flat[k] = p.x;
                flat[k + 1] = p.y;
                flat[k + 2] = p.z;
            }
            k += dim;
            ip += 1;
            if ip > p_upper_row {
                ip = p_lower_row;
            }
        }
        jp += 1;
        if jp > p_upper_col {
            jp = p_lower_col;
        }
    }

    PrepareEvalResult {
        flag_u_or_v: false,
        u1,
        u2,
        d1,
        d2,
        rational,
        dc: PrepareEvalData {
            knots1,
            knots2,
            poles: flat,
        },
    }
}

/// OCCT BSplSLib::BuildCache (BSplSLib.cxx L2389-2547) — the Taylor expansion
/// of the surface on the span of the given indices, normalised to [0, 1] and
/// written into `cache_poles` (the OCCT caller allocates
/// `NCollection_Array2<gp_Pnt>(1, d2+1, 1, d1+1)`, i.e.
/// `cache_poles[(iii-1)][(jjj-1)]` here).
///
/// `u_index` / `v_index` are the OCCT by-value span hints (both 0 at the
/// Geom_OsculatingSurface call site, hence the LocateParameter arm inside
/// [`prepare_eval`]).  `weights == None` is `BSplSLib::NoWeights()`; the
/// rational arm additionally needs `cache_weights` (`CacheWeights != nullptr`).
#[allow(clippy::too_many_arguments)]
fn bspl_slib_build_cache(
    u: f64,
    v: f64,
    u_span_domain: f64,
    v_span_domain: f64,
    u_periodic: bool,
    v_periodic: bool,
    u_degree: i32,
    v_degree: i32,
    u_index: i32,
    v_index: i32,
    u_flat_knots: &[f64],
    v_flat_knots: &[f64],
    poles: &[Vec<DVec3>],
    weights: Option<&[Vec<f64>]>,
    cache_weights: Option<&mut [Vec<f64>]>,
    cache_poles: &mut [Vec<DVec3>],
) {
    let rational_u = weights.is_some();
    let rational_v = weights.is_some();
    let pe = prepare_eval(
        u,
        v,
        u_index,
        v_index,
        u_degree,
        v_degree,
        rational_u,
        rational_v,
        u_periodic,
        v_periodic,
        poles,
        weights,
        u_flat_knots,
        v_flat_knots,
    );
    let flag_u_or_v = pe.flag_u_or_v;
    let u1 = pe.u1;
    let u2 = pe.u2;
    let d1 = pe.d1;
    let d2 = pe.d2;
    let rational = pe.rational;
    let d1p1 = (d1 + 1) as usize;
    let d2p1 = (d2 + 1) as usize;

    let (min_degree_domain, max_degree_domain) = if flag_u_or_v {
        (u_span_domain, v_span_domain)
    } else {
        (v_span_domain, u_span_domain)
    };

    let dim: usize = if rational { 4 } else { 3 };
    let mut fpoles = pe.dc.poles;
    let mut cache_weights = cache_weights;

    // OCCT L2440-2445: the two Bohm passes (U window then, per column, the V
    // window).
    bspl_clib_bohm(u1, d1, d1, &pe.dc.knots1, dim * d2p1, &mut fpoles);
    for kk in 0..=d1 {
        let off = kk as usize * dim * d2p1;
        let end = off + dim * d2p1;
        bspl_clib_bohm(u2, d2, d2, &pe.dc.knots2, dim, &mut fpoles[off..end]);
    }

    let mut factor0 = 1.0f64;
    for ii in 0..=d2 {
        let iii = ii + 1;
        let mut factor1 = 1.0f64;
        for jj in 0..=d1 {
            let jjj = jj + 1;
            let mut index = jj * d2p1 as i32 + ii;
            if rational {
                index <<= 2;
            } else {
                index = (index << 1) + index;
            }
            let f = factor0 * factor1;
            let ix = index as usize;
            cache_poles[(iii - 1) as usize][(jjj - 1) as usize] =
                DVec3::new(f * fpoles[ix], f * fpoles[ix + 1], f * fpoles[ix + 2]);
            if rational {
                let cw = cache_weights
                    .as_deref_mut()
                    .expect("BSplSLib::BuildCache: CacheWeights is null");
                cw[(iii - 1) as usize][(jjj - 1) as usize] = f * fpoles[ix + 3];
            }
            factor1 *= min_degree_domain / jjj as f64;
        }
        factor0 *= max_degree_domain / iii as f64;
    }

    // OCCT L2531-2547: the `Weights != nullptr` locally-polynomial arm — the
    // weight polynomial is set to the constant 1.
    if weights.is_some() && !rational {
        let cw = cache_weights
            .as_deref_mut()
            .expect("BSplSLib::BuildCache: CacheWeights is null");
        for ii in 1..=d2p1 {
            for jj in 1..=d1p1 {
                cw[ii - 1][jj - 1] = 0.0;
            }
        }
        cw[0][0] = 1.0;
    }
}
