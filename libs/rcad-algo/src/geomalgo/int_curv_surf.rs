//! IntCurveSurface polygon / polyhedron sampling classes plus
//! IntPatch_Polyhedron.
//!
//! 1:1 translations:
//!   - IntCurveSurface_ThePolygonOfHInter.hxx/.cxx — the curve sampling
//!     polygon (Closed flag, deflection-over-estimation, ApproxParamOnCurve).
//!   - IntCurveSurface_PolygonUtils.pxx — InitUniform / InitWithParams /
//!     ApproxParamOnCurve.
//!   - IntCurveSurface_ThePolyhedronOfHInter.hxx/.cxx — the surface sampling
//!     polyhedron (singularity flags, PlaneEquation, TriConnex, bounding).
//!   - IntCurveSurface_PolyhedronUtils.pxx — the grid sampling, TriConnex,
//!     PlaneEquation, ComputeMaxDeflection / ComputeMaxBorderDeflection.
//!   - IntPatch_Polyhedron.hxx/.cxx + IntPatch_HInterTool (NbSamplesU/V,
//!     SingularOn*) — the auto-subdivision polyhedron.
//!
//! Points and triangles use the OCCT 1-based indexing convention.

use glam::DVec3;
use rcad_kernel::geom::{CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::precision::INFINITE_VALUE;

/// OCCT IntCurveSurface_PolyhedronUtils THE_MIN_EDGE_LENGTH_SQUARED.
const MIN_EDGE_LENGTH_SQUARED: f64 = 1e-15;

/// OCCT IntPatch_Polyhedron.cxx LONGUEUR_MINI_EDGE_TRIANGLE.
const MIN_EDGE_TRIANGLE: f64 = 1e-14;

/// OCCT IntPatch_Polyhedron.cxx DEFLECTION_COEFF.
const DEFLECTION_COEFF: f64 = 1.1;

/// OCCT IntPatch_Polyhedron.cxx NBMAXUV.
const NB_MAX_UV: usize = 30;

/// OCCT Epsilon(100.) — the ULP gap at 100.0 (IntCurveSurface_ThePolyhedronOfHInter
/// constructor initial TheDeflection).
fn epsilon_100() -> f64 {
    rcad_kernel::math::direct_polynomial_roots::epsilon(100.0)
}

/// OCCT gp::Resolution() = 1e-15.
const GP_RESOLUTION: f64 = 1e-15;

/// OCCT gp_Lin(P1, Dir(P2-P1)).Distance(Pm) — distance from a point to the
/// line through P1 and P2 (the polygon deflection estimate).
fn line_distance(p1: DVec3, p2: DVec3, pm: DVec3) -> f64 {
    let d = p2 - p1;
    let len = d.length();
    if len < GP_RESOLUTION {
        return 0.0;
    }
    (d.cross(pm - p1)).length() / len
}

// ============================================================================
// IntCurveSurface_ThePolygonOfHInter
// ============================================================================

/// OCCT IntCurveSurface_ThePolygonOfHInter — a discretized polygon of a curve.
#[derive(Debug, Clone)]
pub struct ThePolygonOfHInter {
    bnd: BndBox,
    deflection: f64,
    nb_pnt_in: usize,
    points: Vec<DVec3>,
    closed_polygon: bool,
    binf: f64,
    bsup: f64,
    params: Option<Vec<f64>>,
}

impl ThePolygonOfHInter {
    /// OCCT IntCurveSurface_ThePolygonOfHInter(Curve, NbPnt) — the parameter
    /// range is the curve's natural domain (the OCCT Geom_Line endpoints are
    /// ±Precision::Infinite()).
    pub fn new(curve: &dyn CurveEval, nb_pnt: usize) -> Self {
        let nb = if nb_pnt < 5 { 5 } else { nb_pnt };
        let [a, b] = curve.default_domain();
        let binf = if a.is_infinite() { -INFINITE_VALUE } else { a };
        let bsup = if b.is_infinite() { INFINITE_VALUE } else { b };
        let mut poly = ThePolygonOfHInter {
            bnd: BndBox::new(),
            deflection: 0.0,
            nb_pnt_in: nb,
            points: Vec::with_capacity(nb),
            closed_polygon: false,
            binf,
            bsup,
            params: None,
        };
        poly.init_uniform(curve);
        poly
    }

    /// OCCT IntCurveSurface_ThePolygonOfHInter(Curve, U1, U2, NbPnt).
    pub fn new_range(curve: &dyn CurveEval, u1: f64, u2: f64, nb_pnt: usize) -> Self {
        let nb = if nb_pnt < 5 { 5 } else { nb_pnt };
        let mut poly = ThePolygonOfHInter {
            bnd: BndBox::new(),
            deflection: 0.0,
            nb_pnt_in: nb,
            points: Vec::with_capacity(nb),
            closed_polygon: false,
            binf: u1,
            bsup: u2,
            params: None,
        };
        poly.init_uniform(curve);
        poly
    }

    /// OCCT IntCurveSurface_ThePolygonOfHInter(Curve, Upars) — explicit
    /// parameter sampling.
    pub fn new_params(curve: &dyn CurveEval, u_pars: &[f64]) -> Self {
        let nb = u_pars.len();
        let mut poly = ThePolygonOfHInter {
            bnd: BndBox::new(),
            deflection: 0.0,
            nb_pnt_in: nb,
            points: Vec::with_capacity(nb),
            closed_polygon: false,
            binf: u_pars[0],
            bsup: u_pars[nb - 1],
            params: None,
        };
        poly.init_with_params(curve, u_pars);
        poly
    }

    /// OCCT ThePnts(i) / BeginOfSeg(i).
    pub fn begin_of_seg(&self, index: usize) -> DVec3 {
        self.points[index - 1]
    }

    /// OCCT EndOfSeg(i) — ThePnts(i + 1).
    pub fn end_of_seg(&self, index: usize) -> DVec3 {
        self.points[index]
    }

    /// OCCT Closed().
    pub fn closed(&self) -> bool {
        self.closed_polygon
    }

    /// OCCT Closed(flag).
    pub fn set_closed(&mut self, flag: bool) {
        self.closed_polygon = flag;
    }

    /// OCCT NbSegments() — NbPntIn - 1.
    pub fn nb_segments(&self) -> usize {
        self.nb_pnt_in - 1
    }

    /// OCCT InfParameter().
    pub fn inf_parameter(&self) -> f64 {
        self.binf
    }

    /// OCCT SupParameter().
    pub fn sup_parameter(&self) -> f64 {
        self.bsup
    }

    /// OCCT DeflectionOverEstimation().
    pub fn deflection_over_estimation(&self) -> f64 {
        self.deflection
    }

    /// OCCT ApproxParamOnCurve(Index, ParamOnLine) (IntCurveSurface_PolygonUtils.pxx
    /// L155-201).
    pub fn approx_param_on_curve(&self, index: usize, param_on_line: f64) -> f64 {
        if param_on_line < 0.0 || param_on_line > 1.0 {
            return self.binf + (param_on_line * (self.bsup - self.binf)) / (self.nb_pnt_in - 1) as f64;
        }
        let mut index = index;
        let mut param_on_line = param_on_line;
        if index == self.nb_pnt_in && param_on_line == 0.0 {
            index -= 1;
            param_on_line = 1.0;
        }
        let (du, u) = match &self.params {
            None => {
                let du = (self.bsup - self.binf) / (self.nb_pnt_in - 1) as f64;
                (du, self.binf + du * (index - 1) as f64)
            }
            Some(pars) => (pars[index] - pars[index - 1], pars[index - 1]),
        };
        u + du * param_on_line
    }

    /// OCCT Init(Curve) — uniform sampling (IntCurveSurface_PolygonUtils.pxx
    /// L42-87).
    fn init_uniform(&mut self, curve: &dyn CurveEval) {
        let du = (self.bsup - self.binf) / (self.nb_pnt_in - 1) as f64;
        let mut u = self.binf;
        for _ in 0..self.nb_pnt_in {
            let p = curve.point_at(u);
            self.bnd.add_point(p);
            self.points.push(p);
            u += du;
        }

        // Deflection estimate: distance from the midpoints to the chords.
        self.deflection = 0.0;
        if self.nb_pnt_in > 3 {
            u = self.binf + du * 0.5;
            for i in 0..self.nb_pnt_in - 1 {
                let pm = curve.point_at(u);
                let t = line_distance(self.points[i], self.points[i + 1], pm);
                if t > self.deflection {
                    self.deflection = t;
                }
                u += du;
            }
            self.bnd.enlarge(1.5 * self.deflection);
        } else {
            self.bnd.enlarge(1e-10);
        }
        self.closed_polygon = false;
    }

    /// OCCT Init(Curve, Upars) — explicit parameter sampling
    /// (IntCurveSurface_PolygonUtils.pxx L100-144).
    fn init_with_params(&mut self, curve: &dyn CurveEval, u_pars: &[f64]) {
        self.params = Some(u_pars.to_vec());
        for i in 0..self.nb_pnt_in {
            let u = u_pars[i];
            let p = curve.point_at(u);
            self.bnd.add_point(p);
            self.points.push(p);
        }

        self.deflection = 0.0;
        if self.nb_pnt_in > 3 {
            for i in 0..self.nb_pnt_in - 1 {
                let u = 0.5 * (u_pars[i] + u_pars[i + 1]);
                let pm = curve.point_at(u);
                let t = line_distance(self.points[i], self.points[i + 1], pm);
                if t > self.deflection {
                    self.deflection = t;
                }
            }
            self.bnd.enlarge(1.5 * self.deflection);
        } else {
            self.bnd.enlarge(1e-10);
        }
        self.closed_polygon = false;
    }
}

// ============================================================================
// Shared polyhedron helpers (IntCurveSurface_PolyhedronUtils.pxx)
// ============================================================================

/// OCCT IntCurveSurface_PolyhedronUtils::Triangle (L252-261) — the three
/// 1-based point indices of a 1-based triangle.
pub(crate) fn triangle_indices(index: usize, nb_delta_v: usize) -> (usize, usize, usize) {
    let idx = index - 1;
    let line = 1 + idx / (nb_delta_v * 2);
    let colon = 1 + idx % (nb_delta_v * 2);
    let colpnt = (colon + 1) / 2;
    let p1 = (line - 1) * (nb_delta_v + 1) + colpnt;
    let p2 = line * (nb_delta_v + 1) + colpnt + ((colon - 1) % 2);
    let p3 = (line - 1 + (colon % 2)) * (nb_delta_v + 1) + colpnt + 1;
    (p1, p2, p3)
}

/// OCCT IntCurveSurface_PolyhedronUtils::TriConnex (L273-501) — navigate to
/// the triangle connected through the pivot/edge.  Returns (TriCon, OtherP).
pub(crate) fn tri_connex_core(
    triang: usize,
    pivot: usize,
    pedge: usize,
    nb_delta_u: usize,
    nb_delta_v: usize,
) -> (i32, i32) {
    let pivot_m1 = pivot as i64 - 1;
    let nb_delta_v_p1 = (nb_delta_v + 1) as i64;
    let nb_delta_v_m2 = (nb_delta_v + nb_delta_v) as i64;

    let lig_p = pivot_m1 / nb_delta_v_p1;
    let col_p = pivot_m1 - lig_p * nb_delta_v_p1;

    let mut lig_e: i64 = 0;
    let mut col_e: i64 = 0;
    let mut typ_e: i64 = 0;
    if pedge != 0 {
        lig_e = (pedge as i64 - 1) / nb_delta_v_p1;
        col_e = (pedge as i64 - 1) - lig_e * nb_delta_v_p1;
        if lig_p == lig_e {
            typ_e = 1;
        } else if col_p == col_e {
            typ_e = 2;
        } else {
            typ_e = 3;
        }
    }

    let mut lin_t: i64 = 0;
    let mut col_t: i64 = 0;
    let mut lin_o: i64 = 0;
    let mut col_o: i64 = 0;

    if triang != 0 {
        let t = (triang as i64 - 1) / nb_delta_v_m2;
        let tt = (triang as i64 - 1) - t * nb_delta_v_m2;
        lin_t = 1 + t;
        col_t = 1 + tt;
        if typ_e == 0 {
            if lig_p == lin_t {
                lig_e = lig_p - 1;
                col_e = col_p - 1;
                typ_e = 3;
            } else if col_t == lig_p + lig_p {
                lig_e = lig_p;
                col_e = col_p - 1;
                typ_e = 1;
            } else {
                lig_e = lig_p + 1;
                col_e = col_p + 1;
                typ_e = 3;
            }
        }
        match typ_e {
            1 => {
                if lin_t == lig_p {
                    lin_t += 1;
                    lin_o = lig_p + 1;
                    col_o = if col_p > col_e { col_p } else { col_e };
                } else {
                    lin_t -= 1;
                    lin_o = lig_p - 1;
                    col_o = if col_p < col_e { col_p } else { col_e };
                }
            }
            2 => {
                if col_t == col_p + col_p {
                    col_t += 1;
                    lin_o = if lig_p > lig_e { lig_p } else { lig_e };
                    col_o = col_p + 1;
                } else {
                    col_t -= 1;
                    lin_o = if lig_p < lig_e { lig_p } else { lig_e };
                    col_o = col_p - 1;
                }
            }
            _ => {
                if (col_t & 1) == 0 {
                    col_t -= 1;
                    lin_o = if lig_p > lig_e { lig_p } else { lig_e };
                    col_o = if col_p < col_e { col_p } else { col_e };
                } else {
                    col_t += 1;
                    lin_o = if lig_p < lig_e { lig_p } else { lig_e };
                    col_o = if col_p > col_e { col_p } else { col_e };
                }
            }
        }
    } else if pedge == 0 {
        // Unknown triangle and edge.
        lin_t = if 1 > lig_p { 1 } else { lig_p };
        col_t = if 1 > col_p + col_p { 1 } else { col_p + col_p };
        if lig_p == 0 {
            lin_o = lig_p + 1;
        } else {
            lin_o = lig_p - 1;
        }
        col_o = col_p;
    } else {
        // Unknown triangle, known edge — take the left/down connectivity.
        match typ_e {
            1 => {
                lin_t = lig_p + 1;
                col_t = if col_p > col_e { col_p } else { col_e };
                col_t += col_t;
                lin_o = lig_p + 1;
                col_o = if col_p > col_e { col_p } else { col_e };
            }
            2 => {
                lin_t = if lig_p > lig_e { lig_p } else { lig_e };
                col_t = col_p + col_p;
                lin_o = if lig_p < lig_e { lig_p } else { lig_e };
                col_o = col_p - 1;
            }
            _ => {
                lin_t = if lig_p > lig_e { lig_p } else { lig_e };
                col_t = col_p + col_e;
                lin_o = if lig_p > lig_e { lig_p } else { lig_e };
                col_o = if col_p < col_e { col_p } else { col_e };
            }
        }
    }

    let mut tri_con = ((lin_t - 1) * nb_delta_v_m2 + col_t) as i32;

    if lin_t < 1 {
        lin_o = 0;
        col_o = col_p + col_p - col_e;
        if col_o < 0 {
            col_o = 0;
            lin_o = 1;
        } else if col_o > nb_delta_v as i64 {
            col_o = nb_delta_v as i64;
            lin_o = 1;
        }
        tri_con = 0;
    } else if lin_t > nb_delta_u as i64 {
        lin_o = nb_delta_u as i64;
        col_o = col_p + col_p - col_e;
        if col_o < 0 {
            col_o = 0;
            lin_o = nb_delta_u as i64 - 1;
        } else if col_o > nb_delta_v as i64 {
            col_o = nb_delta_v as i64;
            lin_o = nb_delta_u as i64 - 1;
        }
        tri_con = 0;
    }

    if col_t < 1 {
        col_o = 0;
        lin_o = lig_p + lig_p - lig_e;
        if lin_o < 0 {
            lin_o = 0;
            col_o = 1;
        } else if lin_o > nb_delta_u as i64 {
            lin_o = nb_delta_u as i64;
            col_o = 1;
        }
        tri_con = 0;
    } else if col_t > nb_delta_v as i64 {
        col_o = nb_delta_v as i64;
        lin_o = lig_p + lig_p - lig_e;
        if lin_o < 0 {
            lin_o = 0;
            col_o = nb_delta_v as i64 - 1;
        } else if lin_o > nb_delta_u as i64 {
            lin_o = nb_delta_u as i64;
            col_o = nb_delta_v as i64 - 1;
        }
        tri_con = 0;
    }

    let other_p = (lin_o * nb_delta_v_p1 + col_o + 1) as i32;
    (tri_con, other_p)
}

/// OCCT IntCurveSurface_PolyhedronUtils::PlaneEquation (L509-549).
pub(crate) fn plane_equation_of(p1: DVec3, p2: DVec3, p3: DVec3) -> (DVec3, f64) {
    let v1 = p2 - p1;
    let v2 = p3 - p2;
    let v3 = p1 - p3;

    if v1.length_squared() <= MIN_EDGE_LENGTH_SQUARED {
        return (DVec3::X, 0.0);
    }
    if v2.length_squared() <= MIN_EDGE_LENGTH_SQUARED {
        return (DVec3::X, 0.0);
    }
    if v3.length_squared() <= MIN_EDGE_LENGTH_SQUARED {
        return (DVec3::X, 0.0);
    }

    let mut normal = v1.cross(v2) + v2.cross(v3) + v3.cross(v1);
    let norm_len = normal.length();
    if norm_len < GP_RESOLUTION {
        (normal, 0.0)
    } else {
        normal /= norm_len;
        (normal, normal.dot(p1))
    }
}

/// OCCT IntCurveSurface_PolyhedronUtils::Contain (L557-567).
pub(crate) fn contain_in(p1: DVec3, p2: DVec3, p3: DVec3, test_pnt: DVec3) -> bool {
    let v1 = (p2 - p1).cross(test_pnt - p1);
    let v2 = (p3 - p2).cross(test_pnt - p2);
    let v3 = (p1 - p3).cross(test_pnt - p3);
    v1.dot(v2) >= 0.0 && v2.dot(v3) >= 0.0 && v3.dot(v1) >= 0.0
}

/// OCCT IntCurveSurface_PolyhedronUtils::ComputeDeflectionWithCenter
/// (L616-652).
pub(crate) fn deflection_with_center(p1: DVec3, p2: DVec3, p3: DVec3, center: DVec3) -> f64 {
    if p1.distance_squared(p2) <= MIN_EDGE_LENGTH_SQUARED {
        return 0.0;
    }
    if p1.distance_squared(p3) <= MIN_EDGE_LENGTH_SQUARED {
        return 0.0;
    }
    if p2.distance_squared(p3) <= MIN_EDGE_LENGTH_SQUARED {
        return 0.0;
    }
    let xyz1 = p2 - p1;
    let xyz2 = p3 - p2;
    let xyz3 = p1 - p3;
    let mut normal = (xyz1.cross(xyz2)) + (xyz2.cross(xyz3)) + (xyz3.cross(xyz1));
    let norm_len = normal.length();
    if norm_len < GP_RESOLUTION {
        return 0.0;
    }
    normal /= norm_len;
    (p1 - center).dot(normal).abs()
}

/// OCCT IntCurveSurface_PolyhedronUtils::ComputeBorderDeflection (L155-226) —
/// the deflection of the boundary isoline (U or V fixed) sampled at
/// 2·theNbSamples+1 points.
fn compute_border_deflection(
    surface: &dyn SurfaceEval,
    the_parameter: f64,
    p_min: f64,
    p_max: f64,
    is_u_iso: bool,
    nb_samples: usize,
) -> f64 {
    if nb_samples <= 0 {
        return 0.0;
    }
    let a_delta = (p_max - p_min) / nb_samples as f64;
    let nb_pnts = 2 * nb_samples + 1;
    let mut varying = Vec::with_capacity(nb_pnts);
    for i in 0..=nb_samples {
        varying.push(p_min + i as f64 * a_delta);
        if i < nb_samples {
            varying.push(p_min + (i as f64 + 0.5) * a_delta);
        }
    }

    let mut a_deflection = f64::MIN;
    for i in 0..nb_samples {
        // Boundary points at 1-based 2i+1 / 2i+3, midpoint at 2i+2.
        let (a_p1, a_p2, a_par_mid) = if is_u_iso {
            (
                surface.point_at(the_parameter, varying[2 * i]),
                surface.point_at(the_parameter, varying[2 * i + 2]),
                surface.point_at(the_parameter, varying[2 * i + 1]),
            )
        } else {
            (
                surface.point_at(varying[2 * i], the_parameter),
                surface.point_at(varying[2 * i + 2], the_parameter),
                surface.point_at(varying[2 * i + 1], the_parameter),
            )
        };
        let a_p_mid = (a_p2 + a_p1) * 0.5;
        let a_dist = (a_p_mid - a_par_mid).length();
        if a_dist > a_deflection {
            a_deflection = a_dist;
        }
    }
    a_deflection
}

/// OCCT IntCurveSurface_PolyhedronUtils::ComputeMaxBorderDeflection
/// (L706-745) — the four boundary isolines.
fn compute_max_border_deflection(
    surface: &dyn SurfaceEval,
    u0: f64,
    v0: f64,
    u1: f64,
    v1: f64,
    nb_delta_u: usize,
    nb_delta_v: usize,
) -> f64 {
    let mut max_deflection = f64::MIN;
    let d = compute_border_deflection(surface, u0, v0, v1, true, nb_delta_v);
    if d > max_deflection {
        max_deflection = d;
    }
    let d = compute_border_deflection(surface, u1, v0, v1, true, nb_delta_v);
    if d > max_deflection {
        max_deflection = d;
    }
    let d = compute_border_deflection(surface, v0, u0, u1, false, nb_delta_u);
    if d > max_deflection {
        max_deflection = d;
    }
    let d = compute_border_deflection(surface, v1, u0, u1, false, nb_delta_u);
    if d > max_deflection {
        max_deflection = d;
    }
    max_deflection
}

/// OCCT IntCurveSurface_PolyhedronUtils::SetDeflectionOverEstimation
/// (L834-849) — minimum deflection 0.0001.
pub(crate) fn set_deflection_over_estimation(deflection: &mut f64, bnd: &mut BndBox, flec: f64) {
    const MIN_DEFLECTION: f64 = 0.0001;
    if flec < MIN_DEFLECTION {
        *deflection = MIN_DEFLECTION;
        bnd.enlarge(MIN_DEFLECTION);
    } else {
        *deflection = flec;
        bnd.enlarge(flec);
    }
}

/// OCCT IntCurveSurface_PolyhedronUtils::ComputeMaxDeflection (L661-695).
fn compute_max_deflection<S: SurfaceEval + ?Sized>(
    surface: &S,
    nb_delta_u: usize,
    nb_delta_v: usize,
    points: &[DVec3],
    us: &[f64],
    vs: &[f64],
) -> f64 {
    let nb_triangles = nb_delta_u * nb_delta_v * 2;
    if nb_triangles <= 0 {
        return 0.0;
    }
    let mut tol = 0.0;
    for i in 1..=nb_triangles {
        let (i1, i2, i3) = triangle_indices(i, nb_delta_v);
        let (p1, u1, v1) = (points[i1 - 1], us[i1 - 1], vs[i1 - 1]);
        let (p2, u2, v2) = (points[i2 - 1], us[i2 - 1], vs[i2 - 1]);
        let (p3, u3, v3) = (points[i3 - 1], us[i3 - 1], vs[i3 - 1]);
        let u_center = (u1 + u2 + u3) / 3.0;
        let v_center = (v1 + v2 + v3) / 3.0;
        let a_center = surface.point_at(u_center, v_center);
        let tol1 = deflection_with_center(p1, p2, p3, a_center);
        if tol1 > tol {
            tol = tol1;
        }
    }
    tol
}

// ============================================================================
// IntCurveSurface_ThePolyhedronOfHInter
// ============================================================================

/// OCCT IntCurveSurface_ThePolyhedronOfHInter — a triangulated sampling of a
/// surface on a uniform (or explicit) UV grid.
#[derive(Debug, Clone)]
pub struct ThePolyhedronOfHInter {
    nb_delta_u: usize,
    nb_delta_v: usize,
    bnd: BndBox,
    components_bnd: Vec<BndBox>,
    deflection: f64,
    points: Vec<DVec3>,
    us: Vec<f64>,
    vs: Vec<f64>,
    is_on_bounds: Vec<bool>,
    u_min_singular: bool,
    u_max_singular: bool,
    v_min_singular: bool,
    v_max_singular: bool,
    border_deflection: f64,
}

impl ThePolyhedronOfHInter {
    /// OCCT IntCurveSurface_ThePolyhedronOfHInter(Surface, nbdU, nbdV, U1,
    /// V1, U2, V2) — nbdU/nbdV clamped to a minimum of 3.
    pub fn new(
        surface: &dyn SurfaceEval,
        nbdu: usize,
        nbdv: usize,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) -> Self {
        let nbdu = if nbdu < 3 { 3 } else { nbdu };
        let nbdv = if nbdv < 3 { 3 } else { nbdv };
        let mut poly = ThePolyhedronOfHInter {
            nb_delta_u: nbdu,
            nb_delta_v: nbdv,
            bnd: BndBox::new(),
            components_bnd: Vec::new(),
            deflection: epsilon_100(),
            points: Vec::new(),
            us: Vec::new(),
            vs: Vec::new(),
            is_on_bounds: Vec::new(),
            u_min_singular: false,
            u_max_singular: false,
            v_min_singular: false,
            v_max_singular: false,
            border_deflection: 0.0,
        };
        poly.init_uniform(surface, u1, v1, u2, v2);
        poly
    }

    /// OCCT IntCurveSurface_ThePolyhedronOfHInter(Surface, Upars, Vpars) —
    /// explicit parameter arrays.  Panics (Standard_OutOfRange) when either
    /// array holds fewer than two values.
    pub fn new_params(surface: &dyn SurfaceEval, u_pars: &[f64], v_pars: &[f64]) -> Self {
        assert!(
            u_pars.len() >= 2 && v_pars.len() >= 2,
            "IntCurveSurface_ThePolyhedronOfHInter() - parameter arrays must contain at least two values"
        );
        let nbdu = u_pars.len() - 1;
        let nbdv = v_pars.len() - 1;
        let mut poly = ThePolyhedronOfHInter {
            nb_delta_u: nbdu,
            nb_delta_v: nbdv,
            bnd: BndBox::new(),
            components_bnd: Vec::new(),
            deflection: epsilon_100(),
            points: Vec::new(),
            us: Vec::new(),
            vs: Vec::new(),
            is_on_bounds: Vec::new(),
            u_min_singular: false,
            u_max_singular: false,
            v_min_singular: false,
            v_max_singular: false,
            border_deflection: 0.0,
        };
        poly.init_with_params(surface, u_pars, v_pars);
        poly
    }

    /// OCCT Size(nbdu, nbdv).
    pub fn size(&self) -> (usize, usize) {
        (self.nb_delta_u, self.nb_delta_v)
    }

    /// OCCT NbTriangles().
    pub fn nb_triangles(&self) -> usize {
        self.nb_delta_u * self.nb_delta_v * 2
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        (self.nb_delta_u + 1) * (self.nb_delta_v + 1)
    }

    /// OCCT Triangle(Index, P1, P2, P3) — 1-based.
    pub fn triangle(&self, index: usize) -> (usize, usize, usize) {
        triangle_indices(index, self.nb_delta_v)
    }

    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> DVec3 {
        self.points[index - 1]
    }

    /// OCCT Parameters(Index, U, V).
    pub fn parameters(&self, index: usize) -> (f64, f64) {
        (self.us[index - 1], self.vs[index - 1])
    }

    /// OCCT Point(Index, U, V) — the point and its parameters.
    pub fn point_with_parameters(&self, index: usize) -> (DVec3, f64, f64) {
        (self.points[index - 1], self.us[index - 1], self.vs[index - 1])
    }

    /// OCCT TriConnex(Triang, Pivot, Pedge, TriCon, OtherP) — returns
    /// (TriCon, OtherP).
    pub fn tri_connex(&self, triang: usize, pivot: usize, pedge: usize) -> (i32, i32) {
        tri_connex_core(triang, pivot, pedge, self.nb_delta_u, self.nb_delta_v)
    }

    /// OCCT PlaneEquation(Triang, NormalVector, PolarDistance).
    pub fn plane_equation(&self, triang: usize) -> (DVec3, f64) {
        let (i1, i2, i3) = self.triangle(triang);
        plane_equation_of(self.point(i1), self.point(i2), self.point(i3))
    }

    /// OCCT Contain(Triang, ThePnt).
    pub fn contain(&self, triang: usize, pnt: DVec3) -> bool {
        let (i1, i2, i3) = self.triangle(triang);
        contain_in(self.point(i1), self.point(i2), self.point(i3), pnt)
    }

    /// OCCT Bounding().
    pub fn bounding(&self) -> &BndBox {
        &self.bnd
    }

    /// OCCT DeflectionOverEstimation().
    pub fn deflection_over_estimation(&self) -> f64 {
        self.deflection
    }

    /// OCCT HasUMinSingularity() / setters.
    pub fn has_u_min_singularity(&self) -> bool {
        self.u_min_singular
    }
    pub fn has_u_max_singularity(&self) -> bool {
        self.u_max_singular
    }
    pub fn has_v_min_singularity(&self) -> bool {
        self.v_min_singular
    }
    pub fn has_v_max_singularity(&self) -> bool {
        self.v_max_singular
    }
    pub fn set_u_min_singularity(&mut self, sing: bool) {
        self.u_min_singular = sing;
    }
    pub fn set_u_max_singularity(&mut self, sing: bool) {
        self.u_max_singular = sing;
    }
    pub fn set_v_min_singularity(&mut self, sing: bool) {
        self.v_min_singular = sing;
    }
    pub fn set_v_max_singularity(&mut self, sing: bool) {
        self.v_max_singular = sing;
    }

    /// OCCT IsOnBound(Index1, Index2) (PolyhedronUtils.pxx L754-781).
    pub fn is_on_bound(&self, index1: usize, index2: usize) -> bool {
        let diff = (index1 as i64 - index2 as i64).abs();
        if diff != 1 && diff != (self.nb_delta_v + 1) as i64 {
            return false;
        }
        for i in 0..=self.nb_delta_u {
            if index1 == 1 + i * (self.nb_delta_v + 1) && index2 == index1 - 1 {
                return false;
            }
            if index1 == (1 + i) * (self.nb_delta_v + 1) && index2 == index1 + 1 {
                return false;
            }
        }
        self.is_on_bounds[index1 - 1] && self.is_on_bounds[index2 - 1]
    }

    /// OCCT GetBorderDeflection().
    pub fn get_border_deflection(&self) -> f64 {
        self.border_deflection
    }

    /// OCCT Init(Surface, U0, V0, U1, V1) (IntCurveSurface_ThePolyhedronOfHInter.cxx
    /// L107-135) — uniform grid + deflection + bounding.
    fn init_uniform(&mut self, surface: &dyn SurfaceEval, u0: f64, v0: f64, u1: f64, v1: f64) {
        let du = (u1 - u0) / self.nb_delta_u as f64;
        let dv = (v1 - v0) / self.nb_delta_v as f64;
        let nb_u = self.nb_delta_u + 1;
        let nb_v = self.nb_delta_v + 1;

        let mut u_params = Vec::with_capacity(nb_u);
        let mut v_params = Vec::with_capacity(nb_v);
        for i in 0..nb_u {
            u_params.push(u0 + i as f64 * du);
        }
        for j in 0..nb_v {
            v_params.push(v0 + j as f64 * dv);
        }

        for i1 in 0..nb_u {
            for i2 in 0..nb_v {
                let p = surface.point_at(u_params[i1], v_params[i2]);
                self.points.push(p);
                self.us.push(u_params[i1]);
                self.vs.push(v_params[i2]);
                self.is_on_bounds.push(
                    i1 == 0 || i1 == self.nb_delta_u || i2 == 0 || i2 == self.nb_delta_v,
                );
                self.bnd.add_point(p);
            }
        }

        let tol = compute_max_deflection(
            surface,
            self.nb_delta_u,
            self.nb_delta_v,
            &self.points,
            &self.us,
            &self.vs,
        );
        set_deflection_over_estimation(&mut self.deflection, &mut self.bnd, tol * 1.2);
        self.fill_bounding();
        self.border_deflection = compute_max_border_deflection(
            surface, u0, v0, u1, v1, self.nb_delta_u, self.nb_delta_v,
        );
    }

    /// OCCT Init(Surface, Upars, Vpars) (L139-168).
    fn init_with_params(&mut self, surface: &dyn SurfaceEval, u_pars: &[f64], v_pars: &[f64]) {
        for i1 in 0..=self.nb_delta_u {
            for i2 in 0..=self.nb_delta_v {
                let p = surface.point_at(u_pars[i1], v_pars[i2]);
                self.points.push(p);
                self.us.push(u_pars[i1]);
                self.vs.push(v_pars[i2]);
                self.is_on_bounds.push(
                    i1 == 0 || i1 == self.nb_delta_u || i2 == 0 || i2 == self.nb_delta_v,
                );
                self.bnd.add_point(p);
            }
        }

        let tol = compute_max_deflection(
            surface,
            self.nb_delta_u,
            self.nb_delta_v,
            &self.points,
            &self.us,
            &self.vs,
        );
        set_deflection_over_estimation(&mut self.deflection, &mut self.bnd, tol * 1.2);
        self.fill_bounding();
        self.border_deflection = compute_max_border_deflection(
            surface,
            u_pars[0],
            v_pars[0],
            u_pars[u_pars.len() - 1],
            v_pars[v_pars.len() - 1],
            self.nb_delta_u,
            self.nb_delta_v,
        );
    }

    /// OCCT FillBounding() (PolyhedronUtils.pxx L575-608).
    fn fill_bounding(&mut self) {
        self.components_bnd = Vec::with_capacity(self.nb_triangles());
        for i_tri in 1..=self.nb_triangles() {
            let (np1, np2, np3) = self.triangle(i_tri);
            let (p1, p2, p3) = (self.point(np1), self.point(np2), self.point(np3));
            let mut boite = BndBox::new();
            if p1.distance_squared(p2) > MIN_EDGE_LENGTH_SQUARED {
                if p1.distance_squared(p3) > MIN_EDGE_LENGTH_SQUARED {
                    if p2.distance_squared(p3) > MIN_EDGE_LENGTH_SQUARED {
                        boite.add_point(p1);
                        boite.add_point(p2);
                        boite.add_point(p3);
                        boite.enlarge(self.deflection);
                    }
                }
            }
            boite.enlarge(self.deflection);
            self.components_bnd.push(boite);
        }
    }
}

// ============================================================================
// IntPatch_HInterTool sampling (IntPatch_HInterTool.cxx L37-112) + IntPatch_Polyhedron
// ============================================================================

/// Count of distinct values in an expanded knot vector (OCCT NbUKnots).
fn nb_knots(knots: &[f64]) -> usize {
    let mut n = 0usize;
    let mut last = f64::NAN;
    for &k in knots {
        if n == 0 || (k - last).abs() > 1e-15 {
            n += 1;
            last = k;
        }
    }
    n
}

/// OCCT IntPatch_HInterTool::NbSamplesU (IntPatch_HInterTool.cxx L76-112).
fn nb_samples_u(surface: &Surface3) -> usize {
    let n = match surface {
        Surface3::Plane(_) => 2,
        Surface3::Bezier(b) => 3 + b.control_points.len(),
        Surface3::BSpline(b) => {
            let mut nbs = nb_knots(&b.knots_u) * b.degree_u;
            let rational = b.weights.iter().flatten().any(|w| (*w - 1.0).abs() > 1e-12);
            if !rational {
                nbs *= 2;
            }
            if nbs < 4 {
                nbs = 4;
            }
            nbs
        }
        Surface3::Torus(_) => 20,
        Surface3::Cylinder(_)
        | Surface3::Cone(_)
        | Surface3::Sphere(_)
        | Surface3::Revolution(_)
        | Surface3::LinearExtrusion(_)
        | Surface3::Offset(_)
        | Surface3::Ellipsoid(_)
        | Surface3::Helicoid(_)
        | Surface3::Pipe(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::TriBezier(_)
        | Surface3::Trimmed(_) => 10,
    };
    if n > NB_MAX_UV {
        NB_MAX_UV
    } else {
        n
    }
}

/// OCCT IntPatch_HInterTool::NbSamplesV (L37-74).
fn nb_samples_v(surface: &Surface3) -> usize {
    let n = match surface {
        Surface3::Plane(_) => 2,
        Surface3::Bezier(b) => 3 + b.control_points.first().map_or(0, |r| r.len()),
        Surface3::BSpline(b) => {
            let mut nbs = nb_knots(&b.knots_v) * b.degree_v;
            let rational = b.weights.iter().flatten().any(|w| (*w - 1.0).abs() > 1e-12);
            if !rational {
                nbs *= 2;
            }
            if nbs < 4 {
                nbs = 4;
            }
            nbs
        }
        Surface3::Cylinder(_)
        | Surface3::Cone(_)
        | Surface3::Sphere(_)
        | Surface3::Torus(_)
        | Surface3::Revolution(_)
        | Surface3::LinearExtrusion(_)
        | Surface3::Offset(_)
        | Surface3::Ellipsoid(_)
        | Surface3::Helicoid(_)
        | Surface3::Pipe(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::TriBezier(_)
        | Surface3::Trimmed(_) => 15,
        _ => 10,
    };
    if n > NB_MAX_UV {
        NB_MAX_UV
    } else {
        n
    }
}

/// OCCT IntPatch_Polyhedron — the polyhedron used by IntPatch, with the
/// auto-computed subdivision from IntPatch_HInterTool.
#[derive(Debug, Clone)]
pub struct IntPatchPolyhedron {
    bnd: BndBox,
    components_bnd: Vec<BndBox>,
    deflection: f64,
    nb_delta_u: usize,
    nb_delta_v: usize,
    points: Vec<DVec3>,
    us: Vec<f64>,
    vs: Vec<f64>,
    u_min_singular: bool,
    u_max_singular: bool,
    v_min_singular: bool,
    v_max_singular: bool,
}

impl IntPatchPolyhedron {
    /// OCCT IntPatch_Polyhedron(Surface) — auto-computed subdivision.
    pub fn new(surface: &Surface3) -> Self {
        let nbdu = nb_samples_u(surface);
        let nbdv = nb_samples_v(surface);
        Self::new_sub(surface, nbdu, nbdv)
    }

    /// OCCT IntPatch_Polyhedron(Surface, nbu, nbv) — nbu/nbv clamped to a
    /// minimum of 1 (bug #51 fix).
    pub fn new_sub(surface: &Surface3, nbu: usize, nbv: usize) -> Self {
        let nbdu = if nbu < 1 { 1 } else { nbu };
        let nbdv = if nbv < 1 { 1 } else { nbv };
        let [u0, u1, v0, v1] = surface.default_domain();

        let mut poly = IntPatchPolyhedron {
            bnd: BndBox::new(),
            components_bnd: Vec::new(),
            deflection: epsilon_100(),
            nb_delta_u: nbdu,
            nb_delta_v: nbdv,
            points: Vec::new(),
            us: Vec::new(),
            vs: Vec::new(),
            // IntPatch_HInterTool::SingularOnUMin/... — always false.
            u_min_singular: false,
            u_max_singular: false,
            v_min_singular: false,
            v_max_singular: false,
        };

        // Build the UV parameter arrays and evaluate the grid.
        let u_delta = (u1 - u0) / nbdu as f64;
        let v_delta = (v1 - v0) / nbdv as f64;
        for i in 0..=nbdu {
            for j in 0..=nbdv {
                let u = u0 + i as f64 * u_delta;
                let v = v0 + j as f64 * v_delta;
                let p = surface.point_at(u, v);
                poly.points.push(p);
                poly.us.push(u);
                poly.vs.push(v);
                poly.bnd.add_point(p);
            }
        }

        // Compute max deflection via the surface adaptor, scale, bound.
        let mut tol = compute_max_deflection(surface, nbdu, nbdv, &poly.points, &poly.us, &poly.vs);
        tol *= DEFLECTION_COEFF;
        set_deflection_over_estimation(&mut poly.deflection, &mut poly.bnd, tol);
        poly.fill_bounding();
        poly
    }

    /// OCCT Size(nbdu, nbdv).
    pub fn size(&self) -> (usize, usize) {
        (self.nb_delta_u, self.nb_delta_v)
    }

    /// OCCT NbTriangles().
    pub fn nb_triangles(&self) -> usize {
        self.nb_delta_u * self.nb_delta_v * 2
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        (self.nb_delta_u + 1) * (self.nb_delta_v + 1)
    }

    /// OCCT Triangle(Index, P1, P2, P3).
    pub fn triangle(&self, index: usize) -> (usize, usize, usize) {
        triangle_indices(index, self.nb_delta_v)
    }

    /// OCCT Point(Index) — 1-based.
    pub fn point(&self, index: usize) -> DVec3 {
        self.points[index - 1]
    }

    /// OCCT Parameters(Index, U, V).
    pub fn parameters(&self, index: usize) -> (f64, f64) {
        (self.us[index - 1], self.vs[index - 1])
    }

    /// OCCT TriConnex(Triang, Pivot, Pedge, TriCon, OtherP) (IntPatch_Polyhedron.cxx
    /// L299-568) — the polyhedron connectivity with the degenerate-point
    /// guards.  Returns (TriCon, OtherP).
    pub fn tri_connex(&self, triang: usize, pivot: usize, pedge: usize) -> (i32, i32) {
        let (mut tri_con, mut other_p) = tri_connex_core(triang, pivot, pedge, self.nb_delta_u, self.nb_delta_v);

        // OCCT L546-566: Pivot == Pedge (degenerate) — return the triangle
        // unchanged with OtherP = 0.
        if pedge != 0 && self.point(pivot).distance_squared(self.point(pedge)) <= MIN_EDGE_TRIANGLE {
            other_p = 0;
            tri_con = triang as i32;
            return (tri_con, other_p);
        }
        // OtherP == Pedge (degenerate) — return 0 (known uncorrected bug).
        if pedge != 0 && self.point(other_p as usize).distance_squared(self.point(pedge)) <= MIN_EDGE_TRIANGLE {
            return (0, other_p);
        }
        (tri_con, other_p)
    }

    /// OCCT PlaneEquation(Triang, NormalVector, PolarDistance) — uses the
    /// LONGUEUR_MINI_EDGE_TRIANGLE threshold.
    pub fn plane_equation(&self, triang: usize) -> (DVec3, f64) {
        let (i1, i2, i3) = self.triangle(triang);
        let (p1, p2, p3) = (self.point(i1), self.point(i2), self.point(i3));
        let v1 = p2 - p1;
        let v2 = p3 - p2;
        let v3 = p1 - p3;

        if v1.length_squared() <= MIN_EDGE_TRIANGLE {
            return (DVec3::X, 0.0);
        }
        if v2.length_squared() <= MIN_EDGE_TRIANGLE {
            return (DVec3::X, 0.0);
        }
        if v3.length_squared() <= MIN_EDGE_TRIANGLE {
            return (DVec3::X, 0.0);
        }

        let mut normal = v1.cross(v2) + v2.cross(v3) + v3.cross(v1);
        let norm_len = normal.length();
        if norm_len < GP_RESOLUTION {
            (normal, 0.0)
        } else {
            normal /= norm_len;
            (normal, normal.dot(p1))
        }
    }

    /// OCCT Contain(Triang, ThePnt).
    pub fn contain(&self, triang: usize, pnt: DVec3) -> bool {
        let (i1, i2, i3) = self.triangle(triang);
        contain_in(self.point(i1), self.point(i2), self.point(i3), pnt)
    }

    /// OCCT Bounding().
    pub fn bounding(&self) -> &BndBox {
        &self.bnd
    }

    /// OCCT DeflectionOverEstimation().
    pub fn deflection_over_estimation(&self) -> f64 {
        self.deflection
    }

    /// OCCT HasUMinSingularity() etc.
    pub fn has_u_min_singularity(&self) -> bool {
        self.u_min_singular
    }
    pub fn has_u_max_singularity(&self) -> bool {
        self.u_max_singular
    }
    pub fn has_v_min_singularity(&self) -> bool {
        self.v_min_singular
    }
    pub fn has_v_max_singularity(&self) -> bool {
        self.v_max_singular
    }

    /// OCCT FillBounding() (IntPatch_Polyhedron.cxx L246-274).
    fn fill_bounding(&mut self) {
        self.components_bnd = Vec::with_capacity(self.nb_triangles());
        for i_tri in 1..=self.nb_triangles() {
            let (p1, p2, p3) = {
                let (np1, np2, np3) = self.triangle(i_tri);
                (self.point(np1), self.point(np2), self.point(np3))
            };
            let mut boite = BndBox::new();
            if p1.distance_squared(p2) > MIN_EDGE_TRIANGLE {
                if p1.distance_squared(p3) > MIN_EDGE_TRIANGLE {
                    if p2.distance_squared(p3) > MIN_EDGE_TRIANGLE {
                        boite.add_point(p1);
                        boite.add_point(p2);
                        boite.add_point(p3);
                    }
                }
            }
            boite.enlarge(self.deflection);
            self.components_bnd.push(boite);
        }
    }
}
