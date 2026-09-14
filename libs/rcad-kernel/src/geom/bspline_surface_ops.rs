//! OCCT Geom_BSplineSurface surface-level methods (TKG3d/Geom), translated
//! over the rcad `BSplineSurface` data carrier.
//!
//! Sources:
//! - $OCCT_SRC/src/ModelingData/TKG3d/Geom/Geom_BSplineSurface_1.cxx
//!   (UIso/VIso L598-859, UKnot/VKnot L686-698, FirstVKnotIndex/LastVKnotIndex
//!   L1422-1460, LocateU L1464-1514, SetUOrigin/SetVOrigin L1026-1234,
//!   SetUNotPeriodic/SetVNotPeriodic L1238-1346, InsertUKnots/InsertVKnots
//!   L1865-2023, Bounds L1771-1777)
//! - $OCCT_SRC/src/ModelingData/TKG3d/Geom/Geom_BSplineSurface.cxx
//!   (segment L548-853, Segment L857-876, updateUKnots/updateVKnots
//!   L1148-1240, PeriodicNormalization L1244-1296)
//! - $OCCT_SRC/src/FoundationClasses/TKMath/BSplSLib/BSplSLib.cxx
//!   (SetPoles L1916-2017, GetPoles L2021-2124, InsertKnots L2126-2191,
//!   Unperiodize L2331-2380)
//!
//! Architecture difference (data model): OCCT Geom_BSplineSurface stores the
//! compressed (myUKnots, myUMults) arrays as primary data and derives the
//! expanded myUFlatKnots mirror through updateUKnots()/BSplCLib::KnotSequence.
//! The rcad `BSplineSurface` stores only the expanded flat knot vectors
//! (knots_u/knots_v).  The methods that OCCT reads on the compressed arrays
//! derive them here (see `knots_mults_of_flat`); the methods that OCCT reads
//! on the flat mirror (Bounds, the Iso knot sequences) read the stored field
//! directly.  The myUKnotSet / myUSmooth distribution-continuity caches of
//! updateUKnots have no rcad storage (they are recomputable views, not part
//! of the B-spline data), so only the flat mirror rebuild is kept.

use glam::DVec3;

use super::{BSplineCurve3, BSplineSurface};
use crate::base::extrema_ext_elc::epsilon_of;
use crate::core::precision::PCONFUSION;
use crate::math::bspl_lib::{
    at, ati, bspl_slib_iso, first_uknot_index_mults, hunt, insert_knots, knot_analysis,
    knot_sequence, knot_sequence_length, last_uknot_index_mults, locate_parameter_knots_mults,
    pole_index, prepare_insert_knots, prepare_unperiodize, unperiodize, GeomAbsKnotDistribution,
};

// ---------------------------------------------------------------------------
// Architecture shims (the compressed/flat knot duality).
// ---------------------------------------------------------------------------

/// The compressed (myUKnots, myUMults) arrays of one parameter direction,
/// derived from the stored flat knot vector.  OCCT keeps them as member data
/// (Geom_BSplineSurface.hxx myUKnots/myUMults); the rcad carrier stores the
/// flat expansion only, so the compression is the inverse of the
/// updateUKnots()/BSplCLib::KnotSequence expansion (a run-length encode of
/// the flat sequence — exact for the clamped and the periodic forms as
/// stored by the rcad constructors).
fn knots_mults_of_flat(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for (i, k) in flat.iter().enumerate() {
        if i > 0 && *k == flat[i - 1] {
            let last = mults.len() - 1;
            mults[last] += 1;
        } else {
            knots.push(*k);
            mults.push(1);
        }
    }
    (knots, mults)
}

/// The per-direction working state of the surface mutation methods — the
/// OCCT (myUKnots, myUMults, myUPeriodic, myUDeg) member group held locally
/// over the rcad flat-knot carrier (see the module architecture note).
struct SurfDir {
    knots: Vec<f64>, // OCCT myUKnots / myVKnots (compressed)
    mults: Vec<i32>, // OCCT myUMults / myVMults
    periodic: bool,  // OCCT myUPeriodic / myVPeriodic
    degree: usize,   // OCCT myUDeg / myVDeg
}

impl SurfDir {
    /// The compressed state of the stored surface direction.
    fn of(flat: &[f64], degree: usize, periodic: bool) -> Self {
        let (knots, mults) = knots_mults_of_flat(flat);
        SurfDir {
            knots,
            mults,
            periodic,
            degree,
        }
    }
}

/// OCCT Geom_BSplineSurface::FirstUKnotIndex() / FirstVKnotIndex() over the
/// direction state (Geom_BSplineSurface_1.cxx L1408-1432: the periodic arm
/// returns 1, the non-periodic arm is BSplCLib::FirstUKnotIndex).
fn first_knot_index(sd: &SurfDir) -> i32 {
    if sd.periodic {
        1
    } else {
        first_uknot_index_mults(sd.degree, &sd.mults)
    }
}

/// OCCT Geom_BSplineSurface::LastUKnotIndex() / LastVKnotIndex()
/// (Geom_BSplineSurface_1.cxx L1436-1460: the periodic arm returns
/// myUKnots.Length(), the non-periodic arm is BSplCLib::LastUKnotIndex).
fn last_knot_index(sd: &SurfDir) -> i32 {
    if sd.periodic {
        sd.knots.len() as i32
    } else {
        last_uknot_index_mults(sd.degree, &sd.mults)
    }
}

// ---------------------------------------------------------------------------
// BSplSLib pole plumbing (BSplSLib.cxx L1916-2380).
// ---------------------------------------------------------------------------

/// OCCT BSplSLib::SetPoles(Poles, FP, UDirection) (BSplSLib.cxx L1916-1961).
fn set_poles_3d(poles: &[Vec<DVec3>], fp: &mut [f64], u_direction: bool) {
    let mut l = 0usize; // OCCT: l = FP.Lower().
    if u_direction {
        // OCCT L1928-1941: rows outer, cols inner.
        for row in poles {
            for p in row {
                fp[l] = p.x;
                l += 1;
                fp[l] = p.y;
                l += 1;
                fp[l] = p.z;
                l += 1;
            }
        }
    } else {
        // OCCT L1946-1959: cols outer, rows inner.
        for j in 0..poles.first().map(|r| r.len()).unwrap_or(0) {
            for row in poles {
                let p = &row[j];
                fp[l] = p.x;
                l += 1;
                fp[l] = p.y;
                l += 1;
                fp[l] = p.z;
                l += 1;
            }
        }
    }
}

/// OCCT BSplSLib::SetPoles(Poles, Weights, FP, UDirection) (BSplSLib.cxx
/// L1965-2017) — the homogeneous form (dim = 4).
fn set_poles_rational(poles: &[Vec<DVec3>], weights: &[Vec<f64>], fp: &mut [f64], u_direction: bool) {
    let mut l = 0usize;
    if u_direction {
        for (i, row) in poles.iter().enumerate() {
            for (j, p) in row.iter().enumerate() {
                let w = weights[i][j];
                fp[l] = p.x * w;
                l += 1;
                fp[l] = p.y * w;
                l += 1;
                fp[l] = p.z * w;
                l += 1;
                fp[l] = w;
                l += 1;
            }
        }
    } else {
        for j in 0..poles.first().map(|r| r.len()).unwrap_or(0) {
            for (i, row) in poles.iter().enumerate() {
                let p = &row[j];
                let w = weights[i][j];
                fp[l] = p.x * w;
                l += 1;
                fp[l] = p.y * w;
                l += 1;
                fp[l] = p.z * w;
                l += 1;
                fp[l] = w;
                l += 1;
            }
        }
    }
}

/// OCCT BSplSLib::GetPoles(FP, Poles, UDirection) (BSplSLib.cxx L2021-2066).
fn get_poles_3d(fp: &[f64], poles: &mut [Vec<DVec3>], u_direction: bool) {
    let mut l = 0usize;
    if u_direction {
        for row in poles.iter_mut() {
            for p in row.iter_mut() {
                p.x = fp[l];
                l += 1;
                p.y = fp[l];
                l += 1;
                p.z = fp[l];
                l += 1;
            }
        }
    } else {
        for j in 0..poles.first().map(|r| r.len()).unwrap_or(0) {
            for row in poles.iter_mut() {
                let p = &mut row[j];
                p.x = fp[l];
                l += 1;
                p.y = fp[l];
                l += 1;
                p.z = fp[l];
                l += 1;
            }
        }
    }
}

/// OCCT BSplSLib::GetPoles(FP, Poles, Weights, UDirection) (BSplSLib.cxx
/// L2070-2124) — the homogeneous form.
fn get_poles_rational(
    fp: &[f64],
    poles: &mut [Vec<DVec3>],
    weights: &mut [Vec<f64>],
    u_direction: bool,
) {
    let mut l = 0usize;
    if u_direction {
        for (i, row) in poles.iter_mut().enumerate() {
            for (j, p) in row.iter_mut().enumerate() {
                let w = fp[l + 3];
                weights[i][j] = w;
                p.x = fp[l] / w;
                l += 1;
                p.y = fp[l] / w;
                l += 1;
                p.z = fp[l] / w;
                l += 1;
                l += 1;
            }
        }
    } else {
        for j in 0..poles.first().map(|r| r.len()).unwrap_or(0) {
            for (i, row) in poles.iter_mut().enumerate() {
                let p = &mut row[j];
                let w = fp[l + 3];
                weights[i][j] = w;
                p.x = fp[l] / w;
                l += 1;
                p.y = fp[l] / w;
                l += 1;
                p.z = fp[l] / w;
                l += 1;
                l += 1;
            }
        }
    }
}

/// OCCT BSplSLib::InsertKnots(UDirection, Degree, Periodic, Poles, Weights,
/// Knots, Mults, AddKnots, AddMults, NewPoles, NewWeights, NewKnots,
/// NewMults, Epsilon, Add) (BSplSLib.cxx L2126-2191) — flattens the pole
/// grid (SetPoles), runs the curve-level BSplCLib::InsertKnots with the
/// inflated element stride and unflattens (GetPoles).
#[allow(clippy::too_many_arguments)]
fn bspl_slib_insert_knots(
    u_direction: bool,
    degree: usize,
    periodic: bool,
    poles: &[Vec<DVec3>],
    weights: Option<&Vec<Vec<f64>>>,
    knots: &[f64],
    mults: &[i32],
    add_knots: &[f64],
    add_mults: Option<&[i32]>,
    new_poles: &mut [Vec<DVec3>],
    mut new_weights: Option<&mut Vec<Vec<f64>>>,
    new_knots: &mut [f64],
    new_mults: &mut [i32],
    tolerance: f64,
    add: bool,
) {
    let rational = weights.is_some();
    let mut dim: usize = if rational { 4 } else { 3 };

    // OCCT L1949-2150: the flat pole arrays sized dim * RowLength * ColLength
    // (the NewPoles grid already carries the new sizes).
    let row_length = poles.first().map(|r| r.len()).unwrap_or(0);
    let col_length = poles.len();
    let mut fp = vec![0.0f64; dim * row_length * col_length];
    let new_row_length = new_poles.first().map(|r| r.len()).unwrap_or(0);
    let new_col_length = new_poles.len();
    let mut nfp = vec![0.0f64; dim * new_row_length * new_col_length];

    // OCCT L2152-2159.
    if rational {
        set_poles_rational(poles, weights.unwrap(), &mut fp, u_direction);
    } else {
        set_poles_3d(poles, &mut fp, u_direction);
    }

    // OCCT L2161-2168: one curve-level element is a whole pole row/column.
    if u_direction {
        dim *= row_length;
    } else {
        dim *= col_length;
    }

    // OCCT L2169-2181: BSplCLib::InsertKnots(Degree, Periodic, dim, ...).
    insert_knots(
        degree,
        periodic,
        dim,
        &fp,
        knots,
        mults,
        add_knots,
        add_mults,
        &mut nfp,
        new_knots,
        new_mults,
        tolerance,
        add,
    );

    // OCCT L2183-2190.
    if rational {
        get_poles_rational(&nfp, new_poles, new_weights.take().unwrap(), u_direction);
    } else {
        get_poles_3d(&nfp, new_poles, u_direction);
    }
}

/// OCCT BSplSLib::Unperiodize(UDirection, Degree, Mults, Knots, Poles,
/// Weights, NewMults, NewKnots, NewPoles, NewWeights) (BSplSLib.cxx
/// L2331-2380) — the same flatten / BSplCLib::Unperiodize / unflatten form.
/// The OCCT stride computation is kept even though BSplCLib::Unperiodize
/// ignores the Dimension parameter (plain pole copy).
#[allow(clippy::too_many_arguments)]
#[allow(unused_assignments)]
fn bspl_slib_unperiodize(
    u_direction: bool,
    degree: usize,
    mults: &[i32],
    knots: &[f64],
    poles: &[Vec<DVec3>],
    weights: Option<&Vec<Vec<f64>>>,
    new_mults: &mut [i32],
    new_knots: &mut [f64],
    new_poles: &mut [Vec<DVec3>],
    mut new_weights: Option<&mut Vec<Vec<f64>>>,
) {
    let rational = weights.is_some();
    let mut dim: usize = if rational { 4 } else { 3 };

    let row_length = poles.first().map(|r| r.len()).unwrap_or(0);
    let col_length = poles.len();
    let mut fp = vec![0.0f64; dim * row_length * col_length];
    let new_row_length = new_poles.first().map(|r| r.len()).unwrap_or(0);
    let new_col_length = new_poles.len();
    let mut nfp = vec![0.0f64; dim * new_row_length * new_col_length];

    if rational {
        set_poles_rational(poles, weights.unwrap(), &mut fp, u_direction);
    } else {
        set_poles_3d(poles, &mut fp, u_direction);
    }

    if u_direction {
        dim *= row_length;
    } else {
        dim *= col_length;
    }

    unperiodize(degree, mults, knots, &fp, new_mults, new_knots, &mut nfp);

    if rational {
        get_poles_rational(&nfp, new_poles, new_weights.take().unwrap(), u_direction);
    } else {
        get_poles_3d(&nfp, new_poles, u_direction);
    }
}

/// OCCT Geom_BSplineSurface::updateUKnots / updateVKnots
/// (Geom_BSplineSurface.cxx L1148-1240) — rebuilds the flat knot sequence
/// mirror over the compressed arrays.  The KnotAnalysis distribution result
/// and the MaxKnotMult smoothness switch (myUKnotSet / myUSmooth) have no
/// rcad storage (architecture note); only the mirror rebuild is kept.
fn update_knots_flat(sd: &SurfDir) -> Vec<f64> {
    let mut knot_set = GeomAbsKnotDistribution::NonUniform;
    let mut max_knot_mult = 0i32;
    knot_analysis(
        sd.degree,
        sd.periodic,
        &sd.knots,
        &sd.mults,
        &mut knot_set,
        &mut max_knot_mult,
    );

    if knot_set == GeomAbsKnotDistribution::Uniform && !sd.periodic {
        // OCCT L1155-1159: the uniform form stores the compressed array as is.
        sd.knots.clone()
    } else {
        // OCCT L1162-1164: BSplCLib::KnotSequence.
        let mut flat = vec![0.0f64; knot_sequence_length(&sd.mults, sd.degree, sd.periodic)];
        knot_sequence(&sd.knots, &sd.mults, sd.degree, sd.periodic, &mut flat);
        flat
    }
}

// ---------------------------------------------------------------------------
// The surface mutation internals (per direction).
// ---------------------------------------------------------------------------

/// OCCT Geom_BSplineSurface::InsertUKnots (Geom_BSplineSurface_1.cxx
/// L1865-1942) — `u_direction` selects the U form; the V form is the
/// InsertVKnots body (L1946-2023) with the column dimension.  The Add flag
/// is false at the Segment call sites.
#[allow(clippy::too_many_arguments)]
fn insert_dir_knots(
    u_direction: bool,
    sd: &mut SurfDir,
    poles: &mut Vec<Vec<DVec3>>,
    weights: &mut Vec<Vec<f64>>,
    rational: bool,
    the_knots: &[f64],
    the_mults: &[i32],
    the_parametric_tolerance: f64,
) {
    // Check and compute new sizes (OCCT L1870-1885 / L1951-1966).
    let mut nbpoles = 0i32;
    let mut nbknots = 0i32;
    if !prepare_insert_knots(
        sd.degree,
        sd.periodic,
        &sd.knots,
        &sd.mults,
        the_knots,
        Some(the_mults),
        &mut nbpoles,
        &mut nbknots,
        the_parametric_tolerance,
        false,
    ) {
        panic!("Standard_ConstructionError: Geom_BSplineSurface::InsertKnots");
    }

    // myPoles.ColLength() is the U pole count; RowLength the V pole count.
    let nb_u = poles.len() as i32;
    let nb_v = poles.first().map(|r| r.len()).unwrap_or(0) as i32;
    let nb_new = if u_direction { nbpoles == nb_u } else { nbpoles == nb_v };
    if nb_new {
        // OCCT L1887-1890 / L1968-1971.
        return;
    }

    // OCCT L1894-1896 / L1975-1977: the result arrays.
    let (nb_rows, nb_cols) = if u_direction {
        (nbpoles as usize, nb_v as usize)
    } else {
        (nb_u as usize, nbpoles as usize)
    };
    let mut npoles = vec![vec![DVec3::ZERO; nb_cols]; nb_rows];
    let mut nweights = vec![vec![0.0f64; nb_cols]; nb_rows];
    let mut nknots = vec![0.0f64; nbknots as usize];
    let mut nmults = vec![0i32; nbknots as usize];

    // OCCT L1898-1936 / L1979-2017: BSplSLib::InsertKnots on the direction.
    bspl_slib_insert_knots(
        u_direction,
        sd.degree,
        sd.periodic,
        poles,
        if rational { Some(weights) } else { None },
        &sd.knots,
        &sd.mults,
        the_knots,
        Some(the_mults),
        &mut npoles,
        if rational { Some(&mut nweights) } else { None },
        &mut nknots,
        &mut nmults,
        the_parametric_tolerance,
        false,
    );
    // OCCT L1935 / L2016: myWeights = BSplSLib::UnitWeights(...) for the
    // non-rational arm — the rcad grid keeps all-1.0 weights.
    if !rational {
        for row in nweights.iter_mut() {
            for w in row.iter_mut() {
                *w = 1.0;
            }
        }
    }

    *poles = npoles;
    *weights = nweights;
    sd.knots = nknots;
    sd.mults = nmults;
    // OCCT L1941 / L2022: updateUKnots()/updateVKnots() — the flat mirror is
    // rebuilt once at the end of Segment (module architecture note).
}

/// OCCT Geom_BSplineSurface::SetUOrigin (Geom_BSplineSurface_1.cxx
/// L1026-1128) — `u_direction` selects the U form; the V form is the
/// SetVOrigin body (L1132-1234) over the pole columns.
fn set_dir_origin(
    u_direction: bool,
    sd: &mut SurfDir,
    poles: &mut Vec<Vec<DVec3>>,
    weights: &mut Vec<Vec<f64>>,
    rational: bool,
    the_index: i32,
) {
    if !sd.periodic {
        // OCCT L1028-1031 / L1134-1137.
        panic!("Standard_NoSuchObject: SetOrigin: surface is not periodic");
    }

    // OCCT L1036-1037 / L1142-1143.
    let first = first_knot_index(sd);
    let last = last_knot_index(sd);

    if the_index < first || the_index > last {
        // OCCT L1039-1042 / L1145-1148.
        panic!("Standard_DomainError: SetOrigin: Index out of range");
    }

    let nbknots = sd.knots.len(); // myUKnots.Length() / myVKnots.Length()
    // myPoles.ColLength() is the U pole count; RowLength the V pole count.
    let nb_u = poles.len();
    let nb_v = poles.first().map(|r| r.len()).unwrap_or(0);

    // set the knots and mults (OCCT L1051-1065 / L1157-1171).
    let period = at(&sd.knots, last) - at(&sd.knots, first);
    let mut newknots = vec![0.0f64; nbknots];
    let mut newmults = vec![0i32; nbknots];
    let mut k = 1i32;
    for i in the_index..=last {
        set_at_knot(&mut newknots, k, at(&sd.knots, i));
        set_at_mult(&mut newmults, k, ati(&sd.mults, i));
        k += 1;
    }
    for i in (first + 1)..=the_index {
        set_at_knot(&mut newknots, k, at(&sd.knots, i) + period);
        set_at_mult(&mut newmults, k, ati(&sd.mults, i));
        k += 1;
    }

    let mut index = 1i32;
    for i in (first + 1)..=the_index {
        index += ati(&sd.mults, i);
    }

    // set the poles and weights (OCCT L1073-1122 / L1179-1228).  The U form
    // rotates the pole rows, the V form the pole columns.  The OCCT pole
    // `index` is 1-based; the rcad grid is 0-based.
    let index0 = index as usize - 1;
    if u_direction {
        let mut newpoles = vec![vec![DVec3::ZERO; nb_v]; nb_u];
        let mut newweights = vec![vec![0.0f64; nb_v]; nb_u];
        // OCCT L1082-1090: for i = index..=last (rows) / j = 1..=nbvp.
        let mut k = 0usize;
        for i in index0..nb_u {
            for j in 0..nb_v {
                newpoles[k][j] = poles[i][j];
                newweights[k][j] = weights[i][j];
            }
            k += 1;
        }
        // OCCT L1091-1099: for i = first..<index (first = myPoles.LowerRow()).
        for i in 0..index0 {
            for j in 0..nb_v {
                newpoles[k][j] = poles[i][j];
                newweights[k][j] = weights[i][j];
            }
            k += 1;
        }
        *poles = newpoles;
        *weights = newweights;
    } else {
        let mut newpoles = vec![vec![DVec3::ZERO; nb_v]; nb_u];
        let mut newweights = vec![vec![0.0f64; nb_v]; nb_u];
        // OCCT L1188-1205: for j = index..=last (cols) / i = 1..=nbup.
        let mut k = 0usize;
        for j in index0..nb_v {
            for i in 0..nb_u {
                newpoles[i][k] = poles[i][j];
                newweights[i][k] = weights[i][j];
            }
            k += 1;
        }
        // OCCT L1197-1205: for j = first..<index (first = myPoles.LowerCol()).
        for j in 0..index0 {
            for i in 0..nb_u {
                newpoles[i][k] = poles[i][j];
                newweights[i][k] = weights[i][j];
            }
            k += 1;
        }
        *poles = newpoles;
        *weights = newweights;
    }
    if !rational {
        // OCCT L1121 / L1227: myWeights = BSplSLib::UnitWeights(...).
        for row in weights.iter_mut() {
            for w in row.iter_mut() {
                *w = 1.0;
            }
        }
    }

    sd.knots = newknots;
    sd.mults = newmults;
    // OCCT L1127 / L1233: updateUKnots()/updateVKnots() — rebuilt at the end
    // of Segment (module architecture note).
}

/// OCCT Geom_BSplineSurface::SetUNotPeriodic (Geom_BSplineSurface_1.cxx
/// L1238-1290) — `u_direction` selects the U form; the V form is the
/// SetVNotPeriodic body (L1294-1346).
fn set_dir_not_periodic(
    u_direction: bool,
    sd: &mut SurfDir,
    poles: &mut Vec<Vec<DVec3>>,
    weights: &mut Vec<Vec<f64>>,
    rational: bool,
) {
    // OCCT L1240 / L1296: the if (myUPeriodic) / if (myVPeriodic) guard.
    if !sd.periodic {
        return;
    }

    // OCCT L1243-1250 / L1299-1306.
    let mut nb_knots = 0i32;
    let mut nb_poles = 0i32;
    prepare_unperiodize(sd.degree, &sd.mults, &mut nb_knots, &mut nb_poles);

    let nb_u = poles.len();
    let nb_v = poles.first().map(|r| r.len()).unwrap_or(0);
    let (nb_rows, nb_cols) = if u_direction {
        (nb_poles as usize, nb_v)
    } else {
        (nb_u, nb_poles as usize)
    };
    let mut npoles = vec![vec![DVec3::ZERO; nb_cols]; nb_rows];
    let mut nweights = vec![vec![0.0f64; nb_cols]; nb_rows];
    let mut nknots = vec![0.0f64; nb_knots as usize];
    let mut nmults = vec![0i32; nb_knots as usize];

    // OCCT L1252-1281 / L1308-1337.
    bspl_slib_unperiodize(
        u_direction,
        sd.degree,
        &sd.mults,
        &sd.knots,
        poles,
        if rational { Some(weights) } else { None },
        &mut nmults,
        &mut nknots,
        &mut npoles,
        if rational { Some(&mut nweights) } else { None },
    );
    if !rational {
        // OCCT L1280 / L1336: BSplSLib::UnitWeights(...).
        for row in nweights.iter_mut() {
            for w in row.iter_mut() {
                *w = 1.0;
            }
        }
    }

    *poles = npoles;
    *weights = nweights;
    sd.mults = nmults;
    sd.knots = nknots;
    sd.periodic = false;
    // OCCT L1288 / L1344: updateUKnots()/updateVKnots() — rebuilt at the end
    // of Segment (module architecture note).
}

// 1-based knot/mult writes (the OCCT SetValue forms of SetUOrigin).
fn set_at_knot(arr: &mut [f64], i: i32, v: f64) {
    arr[(i - 1) as usize] = v;
}
fn set_at_mult(arr: &mut [i32], i: i32, v: i32) {
    arr[(i - 1) as usize] = v;
}

impl BSplineSurface {
    /// OCCT Geom_BSplineSurface::Bounds(U1, U2, V1, V2)
    /// (Geom_BSplineSurface_1.cxx L1771-1777) — reads on the flat knot
    /// mirrors (myUFlatKnots/myVFlatKnots).
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        let u1 = self.knots_u[self.degree_u]; // Value(myUDeg + 1)
        let u2 = self.knots_u[self.knots_u.len() - 1 - self.degree_u]; // Value(Upper() - myUDeg)
        let v1 = self.knots_v[self.degree_v]; // Value(myVDeg + 1)
        let v2 = self.knots_v[self.knots_v.len() - 1 - self.degree_v]; // Value(Upper() - myVDeg)
        (u1, u2, v1, v2)
    }

    /// OCCT Geom_BSplineSurface::UKnot(UIndex)
    /// (Geom_BSplineSurface_1.cxx L686-690).
    pub fn u_knot(&self, the_i: i32) -> f64 {
        let (uknots, _) = knots_mults_of_flat(&self.knots_u);
        at(&uknots, the_i)
    }

    /// OCCT Geom_BSplineSurface::VKnot(VIndex)
    /// (Geom_BSplineSurface_1.cxx L694-698).
    pub fn v_knot(&self, the_i: i32) -> f64 {
        let (vknots, _) = knots_mults_of_flat(&self.knots_v);
        at(&vknots, the_i)
    }

    /// OCCT Geom_BSplineSurface::FirstVKnotIndex()
    /// (Geom_BSplineSurface_1.cxx L1422-1432).
    pub fn first_v_knot_index(&self) -> i32 {
        if self.is_periodic_v {
            1
        } else {
            let (_, vmults) = knots_mults_of_flat(&self.knots_v);
            first_uknot_index_mults(self.degree_v, &vmults)
        }
    }

    /// OCCT Geom_BSplineSurface::LastVKnotIndex()
    /// (Geom_BSplineSurface_1.cxx L1450-1460).
    pub fn last_v_knot_index(&self) -> i32 {
        if self.is_periodic_v {
            self.v_knot_count()
        } else {
            let (_, vmults) = knots_mults_of_flat(&self.knots_v);
            last_uknot_index_mults(self.degree_v, &vmults)
        }
    }

    /// The myVKnots.Length() of the stored V direction (the compressed knot
    /// count; the periodic arm of LastVKnotIndex,
    /// Geom_BSplineSurface_1.cxx L1452-1455).
    fn v_knot_count(&self) -> i32 {
        let (vknots, _) = knots_mults_of_flat(&self.knots_v);
        vknots.len() as i32
    }

    /// OCCT Geom_BSplineSurface::LocateU(U, ParametricTolerance, I1, I2)
    /// (Geom_BSplineSurface_1.cxx L1464-1514) — the hxx default form
    /// (WithKnotRepetition = false, Geom_BSplineSurface.hxx L679-684).
    pub fn locate_u(&self, the_u: f64, the_parametric_tolerance: f64, the_i1: &mut i32, the_i2: &mut i32) {
        self.locate_u_opt(the_u, the_parametric_tolerance, the_i1, the_i2, false);
    }

    /// OCCT Geom_BSplineSurface::LocateU(U, ParametricTolerance, I1, I2,
    /// WithKnotRepetition) (Geom_BSplineSurface_1.cxx L1464-1514).
    pub fn locate_u_opt(
        &self,
        the_u: f64,
        the_parametric_tolerance: f64,
        the_i1: &mut i32,
        the_i2: &mut i32,
        with_knot_repetition: bool,
    ) {
        // OCCT L1470: NewU = U, vbid = myVKnots.Value(1).
        let mut new_u = the_u;
        let (vknots, _) = knots_mults_of_flat(&self.knots_v);
        let mut vbid = at(&vknots, 1);
        // OCCT L1471: TheKnots = WithKnotRepetition ? myUFlatKnots : myUKnots.
        let uknots_flat: Vec<f64>;
        let the_knots: &[f64] = if with_knot_repetition {
            &self.knots_u
        } else {
            let (k, _) = knots_mults_of_flat(&self.knots_u);
            uknots_flat = k;
            &uknots_flat
        };

        // OCCT L1473: PeriodicNormalization(NewU, vbid).
        self.periodic_normalization(&mut new_u, &mut vbid);

        let knot_len = the_knots.len() as i32;
        let u_first = at(the_knots, 1);
        let u_last = at(the_knots, knot_len);
        let p_parametric_tolerance = the_parametric_tolerance.abs();
        if (new_u - u_first).abs() <= p_parametric_tolerance {
            // OCCT L1478-1481.
            *the_i1 = 1;
            *the_i2 = 1;
        } else if (new_u - u_last).abs() <= p_parametric_tolerance {
            // OCCT L1482-1485.
            *the_i1 = knot_len;
            *the_i2 = knot_len;
        } else if new_u < u_first {
            // OCCT L1486-1490.
            *the_i2 = 1;
            *the_i1 = 0;
        } else if new_u > u_last {
            // OCCT L1491-1495.
            *the_i1 = knot_len;
            *the_i2 = *the_i1 + 1;
        } else {
            // OCCT L1496-1513.
            *the_i1 = 1;
            hunt(the_knots, new_u, the_i1);
            *the_i1 = (*the_i1).min(knot_len).max(1);
            while *the_i1 + 1 <= knot_len && (at(the_knots, *the_i1 + 1) - new_u).abs() <= p_parametric_tolerance {
                *the_i1 += 1;
            }
            if (at(the_knots, *the_i1) - new_u).abs() <= p_parametric_tolerance {
                *the_i2 = *the_i1;
            } else {
                *the_i2 = *the_i1 + 1;
            }
        }
    }

    /// OCCT Geom_BSplineSurface::PeriodicNormalization(Uparameter,
    /// Vparameter) (Geom_BSplineSurface.cxx L1244-1296).
    fn periodic_normalization(&self, u_parameter: &mut f64, v_parameter: &mut f64) {
        if self.is_periodic_u {
            // OCCT L1250-1251: the flat-knot range bounds.
            let a_max_val = self.knots_u[self.knots_u.len() - 1 - self.degree_u];
            let a_min_val = self.knots_u[self.degree_u];
            let eps = epsilon_of(*u_parameter).abs();
            let period = a_max_val - a_min_val;

            if period <= eps {
                panic!("Standard_OutOfRange: PeriodicNormalization: Uparameter is too great number");
            }

            let is_less = a_min_val - *u_parameter > 0.0;
            let is_greater = *u_parameter - a_max_val > 0.0;
            if is_less || is_greater {
                // OCCT L1266-1269: modf(aDPar / Period, &aNbPer) — only the
                // integral part is consumed.
                let a_d_par = if is_less {
                    a_max_val - *u_parameter
                } else {
                    a_min_val - *u_parameter
                };
                let a_nb_per = (a_d_par / period).trunc();
                *u_parameter += a_nb_per * period;
            }
        }
        if self.is_periodic_v {
            // OCCT L1274-1275.
            let a_max_val = self.knots_v[self.knots_v.len() - 1 - self.degree_v];
            let a_min_val = self.knots_v[self.degree_v];
            let eps = epsilon_of(*v_parameter).abs();
            let period = a_max_val - a_min_val;

            if period <= eps {
                panic!("Standard_OutOfRange: PeriodicNormalization: Vparameter is too great number");
            }

            let is_less = a_min_val - *v_parameter > 0.0;
            let is_greater = *v_parameter - a_max_val > 0.0;
            if is_less || is_greater {
                // OCCT L1290-1293.
                let a_d_par = if is_less {
                    a_max_val - *v_parameter
                } else {
                    a_min_val - *v_parameter
                };
                let a_nb_per = (a_d_par / period).trunc();
                *v_parameter += a_nb_per * period;
            }
        }
    }

    /// OCCT Geom_BSplineSurface::UIso(U) (Geom_BSplineSurface_1.cxx
    /// L598-635) — the BSplSLib::Iso U arm over the flat U knots, the result
    /// curve built over the V direction data.  The rcad curve stores the
    /// expanded knot vector == the stored knots_v (myVFlatKnots).
    pub fn u_iso(&self, the_u: f64) -> BSplineCurve3 {
        let rational = self.is_rational_u() || self.is_rational_v();
        let weights = if rational { Some(&self.weights) } else { None };
        let (cpoles, cweights) = bspl_slib_iso(
            the_u,
            true,
            self.degree_u,
            &self.knots_u,
            &self.control_points,
            weights.map(|w| &w[..]),
            self.is_periodic_u,
        );
        BSplineCurve3 {
            degree: self.degree_v,
            knots: self.knots_v.clone(),
            control_points: cpoles,
            weights: cweights,
            is_periodic: self.is_periodic_v,
        }
    }

    /// OCCT Geom_BSplineSurface::VIso(V) (Geom_BSplineSurface_1.cxx
    /// L775-812) — the BSplSLib::Iso V arm over the flat V knots, the result
    /// curve built over the U direction data.
    pub fn v_iso(&self, the_v: f64) -> BSplineCurve3 {
        let rational = self.is_rational_u() || self.is_rational_v();
        let weights = if rational { Some(&self.weights) } else { None };
        let (cpoles, cweights) = bspl_slib_iso(
            the_v,
            false,
            self.degree_v,
            &self.knots_v,
            &self.control_points,
            weights.map(|w| &w[..]),
            self.is_periodic_v,
        );
        BSplineCurve3 {
            degree: self.degree_u,
            knots: self.knots_u.clone(),
            control_points: cpoles,
            weights: cweights,
            is_periodic: self.is_periodic_u,
        }
    }

    /// OCCT Geom_BSplineSurface::Segment(U1, U2, V1, V2) — the hxx default
    /// form (Geom_BSplineSurface.hxx L583-588: theUTolerance = theVTolerance
    /// = Precision::PConfusion()) over the 6-arg overload
    /// (Geom_BSplineSurface.cxx L857-876).
    pub fn segment(&mut self, the_u1: f64, the_u2: f64, the_v1: f64, the_v2: f64) {
        self.segment_tol(the_u1, the_u2, the_v1, the_v2, PCONFUSION, PCONFUSION);
    }

    /// OCCT Geom_BSplineSurface::Segment(U1, U2, V1, V2, theUTolerance,
    /// theVTolerance) (Geom_BSplineSurface.cxx L857-876).
    pub fn segment_tol(
        &mut self,
        the_u1: f64,
        the_u2: f64,
        the_v1: f64,
        the_v2: f64,
        the_u_tolerance: f64,
        the_v_tolerance: f64,
    ) {
        if (the_u2 < the_u1) || (the_v2 < the_v1) {
            // OCCT L864-867.
            panic!("Standard_DomainError: Geom_BSplineSurface::Segment");
        }

        // OCCT L869-870.
        let a_max_u = the_u2.abs().max(the_u1.abs());
        let eps_u = epsilon_of(a_max_u).max(the_u_tolerance);

        // OCCT L872-873.
        let a_max_v = the_v2.abs().max(the_v1.abs());
        let eps_v = epsilon_of(a_max_v).max(the_v_tolerance);

        // OCCT L875: segment(U1, U2, V1, V2, EpsU, EpsV, true, true).
        self.segment_eps(the_u1, the_u2, the_v1, the_v2, eps_u, eps_v, true, true);
    }

    /// OCCT Geom_BSplineSurface::segment(U1, U2, V1, V2, EpsU, EpsV,
    /// SegmentInU, SegmentInV) (Geom_BSplineSurface.cxx L548-853).
    #[allow(clippy::too_many_arguments)]
    fn segment_eps(
        &mut self,
        the_u1: f64,
        the_u2: f64,
        the_v1: f64,
        the_v2: f64,
        eps_u: f64,
        eps_v: f64,
        segment_in_u: bool,
        segment_in_v: bool,
    ) {
        // OCCT L557: ClearEvalRepresentation() — the rcad surface has no
        // eval representation (architecture difference; no-op).

        // OCCT L558-570: the deltaU period clamp.
        let mut su = SurfDir::of(&self.knots_u, self.degree_u, self.is_periodic_u);
        let mut delta_u = the_u2 - the_u1;
        if su.periodic {
            let a_u_period = at(&su.knots, su.knots.len() as i32) - at(&su.knots, 1);
            if delta_u - a_u_period > PCONFUSION {
                panic!("Standard_DomainError: Geom_BSplineSurface::Segment");
            }
            if delta_u > a_u_period {
                delta_u = a_u_period;
            }
        }

        // OCCT L572-584: the deltaV period clamp.
        let mut sv = SurfDir::of(&self.knots_v, self.degree_v, self.is_periodic_v);
        let mut delta_v = the_v2 - the_v1;
        if sv.periodic {
            let a_v_period = at(&sv.knots, sv.knots.len() as i32) - at(&sv.knots, 1);
            if delta_v - a_v_period > PCONFUSION {
                panic!("Standard_DomainError: Geom_BSplineSurface::Segment");
            }
            if delta_v > a_v_period {
                delta_v = a_v_period;
            }
        }

        let rational = self.is_rational_u() || self.is_rational_v();
        let mut poles = self.control_points.clone();
        let mut weights = self.weights.clone();

        // OCCT L586-588: NewU1, NewU2, NewV1, NewV2; U, V; indexU, indexV.
        let mut new_u1;
        let mut new_u2;
        let mut new_v1;
        let mut new_v2;
        let mut u;
        let mut v;
        let mut index_u;

        // OCCT L590-599: LocateParameter on U1.
        index_u = 0;
        new_u1 = 0.0;
        locate_parameter_knots_mults(
            su.degree,
            &su.knots,
            &su.mults,
            the_u1,
            su.periodic,
            1,
            su.knots.len() as i32,
            &mut index_u,
            &mut new_u1,
        );
        // OCCT L600-609: LocateParameter on U2.
        index_u = 0;
        new_u2 = 0.0;
        locate_parameter_knots_mults(
            su.degree,
            &su.knots,
            &su.mults,
            the_u2,
            su.periodic,
            1,
            su.knots.len() as i32,
            &mut index_u,
            &mut new_u2,
        );
        // OCCT L610-620: inserting the UKnots.
        if segment_in_u {
            let u_knots = [new_u1.min(new_u2), new_u1.max(new_u2)];
            let u_mults = [su.degree as i32, su.degree as i32];
            insert_dir_knots(
                true,
                &mut su,
                &mut poles,
                &mut weights,
                rational,
                &u_knots,
                &u_mults,
                eps_u,
            );
        }

        // OCCT L622-641: LocateParameter on V1/V2.
        let mut index_v = 0;
        new_v1 = 0.0;
        locate_parameter_knots_mults(
            sv.degree,
            &sv.knots,
            &sv.mults,
            the_v1,
            sv.periodic,
            1,
            sv.knots.len() as i32,
            &mut index_v,
            &mut new_v1,
        );
        index_v = 0;
        new_v2 = 0.0;
        locate_parameter_knots_mults(
            sv.degree,
            &sv.knots,
            &sv.mults,
            the_v2,
            sv.periodic,
            1,
            sv.knots.len() as i32,
            &mut index_v,
            &mut new_v2,
        );
        // OCCT L642-652: inserting the VKnots.
        if segment_in_v {
            let v_knots = [new_v1.min(new_v2), new_v1.max(new_v2)];
            let v_mults = [sv.degree as i32, sv.degree as i32];
            insert_dir_knots(
                false,
                &mut sv,
                &mut poles,
                &mut weights,
                rational,
                &v_knots,
                &v_mults,
                eps_v,
            );
        }

        // OCCT L654-672: set the origine at NewU1 (periodic U arm).
        if su.periodic && segment_in_u {
            let mut index = 0i32;
            u = 0.0;
            locate_parameter_knots_mults(
                su.degree,
                &su.knots,
                &su.mults,
                the_u1,
                su.periodic,
                1,
                su.knots.len() as i32,
                &mut index,
                &mut u,
            );
            if (at(&su.knots, index + 1) - u).abs() <= eps_u {
                index += 1;
            }
            set_dir_origin(true, &mut su, &mut poles, &mut weights, rational, index);
            set_dir_not_periodic(true, &mut su, &mut poles, &mut weights, rational);
        }

        // OCCT L674-677: compute index1 and index2 to set the new knots and mults.
        let mut index1_u = 0i32;
        let mut index2_u = 0i32;
        let from_u1 = 1i32;
        let to_u2 = su.knots.len() as i32;
        // OCCT L678-686.
        u = 0.0;
        locate_parameter_knots_mults(
            su.degree,
            &su.knots,
            &su.mults,
            new_u1,
            su.periodic,
            from_u1,
            to_u2,
            &mut index1_u,
            &mut u,
        );
        // OCCT L687-690.
        if (at(&su.knots, index1_u + 1) - u).abs() <= eps_u {
            index1_u += 1;
        }
        // OCCT L691-699.
        locate_parameter_knots_mults(
            su.degree,
            &su.knots,
            &su.mults,
            new_u1 + delta_u,
            su.periodic,
            from_u1,
            to_u2,
            &mut index2_u,
            &mut u,
        );
        // OCCT L700-703.
        if (at(&su.knots, index2_u + 1) - u).abs() <= eps_u || index2_u == index1_u {
            index2_u += 1;
        }

        // OCCT L705-716: the new U knots and mults window.
        let nbuknots = (index2_u - index1_u + 1) as usize;
        let mut nuknots = vec![0.0f64; nbuknots];
        let mut numults = vec![0i32; nbuknots];
        let mut k = 1i32;
        for i in index1_u..=index2_u {
            set_at_knot(&mut nuknots, k, at(&su.knots, i));
            set_at_mult(&mut numults, k, ati(&su.mults, i));
            k += 1;
        }
        // OCCT L717-721.
        if segment_in_u {
            set_at_mult(&mut numults, 1, su.degree as i32 + 1);
            set_at_mult(&mut numults, nbuknots as i32, su.degree as i32 + 1);
        }

        // OCCT L723-741: set the origine at NewV1 (periodic V arm).
        if sv.periodic && segment_in_v {
            let mut index = 0i32;
            v = 0.0;
            locate_parameter_knots_mults(
                sv.degree,
                &sv.knots,
                &sv.mults,
                the_v1,
                sv.periodic,
                1,
                sv.knots.len() as i32,
                &mut index,
                &mut v,
            );
            if (at(&sv.knots, index + 1) - v).abs() <= eps_v {
                index += 1;
            }
            set_dir_origin(false, &mut sv, &mut poles, &mut weights, rational, index);
            set_dir_not_periodic(false, &mut sv, &mut poles, &mut weights, rational);
        }

        // OCCT L743-746: compute index1 and index2 to set the new knots and mults.
        let mut index1_v = 0i32;
        let mut index2_v = 0i32;
        let from_v1 = 1i32;
        let to_v2 = sv.knots.len() as i32;
        // OCCT L747-755.
        v = 0.0;
        locate_parameter_knots_mults(
            sv.degree,
            &sv.knots,
            &sv.mults,
            new_v1,
            sv.periodic,
            from_v1,
            to_v2,
            &mut index1_v,
            &mut v,
        );
        // OCCT L756-759.
        if (at(&sv.knots, index1_v + 1) - v).abs() <= eps_v {
            index1_v += 1;
        }
        // OCCT L760-768.
        locate_parameter_knots_mults(
            sv.degree,
            &sv.knots,
            &sv.mults,
            new_v1 + delta_v,
            sv.periodic,
            from_v1,
            to_v2,
            &mut index2_v,
            &mut v,
        );
        // OCCT L769-772.
        if (at(&sv.knots, index2_v + 1) - v).abs() <= eps_v || index2_v == index1_v {
            index2_v += 1;
        }

        // OCCT L774-785: the new V knots and mults window.
        let nbvknots = (index2_v - index1_v + 1) as usize;
        let mut nvknots = vec![0.0f64; nbvknots];
        let mut nvmults = vec![0i32; nbvknots];
        let mut k = 1i32;
        for i in index1_v..=index2_v {
            set_at_knot(&mut nvknots, k, at(&sv.knots, i));
            set_at_mult(&mut nvmults, k, ati(&sv.mults, i));
            k += 1;
        }
        // OCCT L786-790.
        if segment_in_v {
            set_at_mult(&mut nvmults, 1, sv.degree as i32 + 1);
            set_at_mult(&mut nvmults, nbvknots as i32, sv.degree as i32 + 1);
        }

        // OCCT L792-799: compute index1 and index2 to set the new poles and
        // weights (U direction).
        let mut pindex1_u = pole_index(su.degree, index1_u, su.periodic, &su.mults);
        let mut pindex2_u = pole_index(su.degree, index2_u, su.periodic, &su.mults);

        pindex1_u += 1;
        let nb_u = poles.len() as i32; // myPoles.ColLength()
        pindex2_u = (pindex2_u + 1).min(nb_u);

        let nbupoles = (pindex2_u - pindex1_u + 1) as usize;

        // OCCT L801-808: the V direction pole window.
        let mut pindex1_v = pole_index(sv.degree, index1_v, sv.periodic, &sv.mults);
        let mut pindex2_v = pole_index(sv.degree, index2_v, sv.periodic, &sv.mults);

        pindex1_v += 1;
        let nb_v = poles.first().map(|r| r.len()).unwrap_or(0) as i32; // myPoles.RowLength()
        pindex2_v = (pindex2_v + 1).min(nb_v);

        let nbvpoles = (pindex2_v - pindex1_v + 1) as usize;

        // OCCT L810-842: the pole (and weight) window copy.  The OCCT
        // pindex values are 1-based pole indices; the rcad grid is 0-based.
        let mut npoles = vec![vec![DVec3::ZERO; nbvpoles]; nbupoles];
        let mut nweights = vec![vec![0.0f64; nbvpoles]; nbupoles];
        let mut k = 0usize;
        if rational {
            for i in (pindex1_u as usize)..=(pindex2_u as usize) {
                let mut l = 0usize;
                for j in (pindex1_v as usize)..=(pindex2_v as usize) {
                    npoles[k][l] = poles[i - 1][j - 1];
                    nweights[k][l] = weights[i - 1][j - 1];
                    l += 1;
                }
                k += 1;
            }
        } else {
            for i in (pindex1_u as usize)..=(pindex2_u as usize) {
                let mut l = 0usize;
                for j in (pindex1_v as usize)..=(pindex2_v as usize) {
                    npoles[k][l] = poles[i - 1][j - 1];
                    l += 1;
                }
                k += 1;
            }
            // OCCT L841: myWeights = BSplSLib::UnitWeights(nbupoles, nbvpoles).
            for row in nweights.iter_mut() {
                for w in row.iter_mut() {
                    *w = 1.0;
                }
            }
        }

        // OCCT L844-848: commit the new data.
        su.knots = nuknots;
        su.mults = numults;
        sv.knots = nvknots;
        sv.mults = nvmults;
        poles = npoles;
        weights = nweights;

        // OCCT L850: myMaxDerivInvOk = false — no rcad deriv cache
        // (architecture difference; no-op).

        // OCCT L851-852: updateUKnots(); updateVKnots() — the flat mirror
        // rebuild over the stored rcad fields.
        self.knots_u = update_knots_flat(&su);
        self.knots_v = update_knots_flat(&sv);
        self.degree_u = su.degree;
        self.degree_v = sv.degree;
        self.is_periodic_u = su.periodic;
        self.is_periodic_v = sv.periodic;
        self.control_points = poles;
        self.weights = weights;
    }
}
