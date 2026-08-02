// OCCT IntTools_FaceFace — face-face intersection.
//
// Computes intersection curves between two surfaces.
// Handles analytic cases: Plane-Plane (PerformPlanes via IntAna_QuadQuadGeo
// and ClassifyLin2d), Plane-Sphere, with the general IntPatch path to follow.

use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, Line2d, Line3, Plane, Surface3, CurveEval, SurfaceEval,
};
use glam::{DVec2, DVec3};

use super::face_make_curve;
use crate::geomalgo::int_patch::IntPatchLine;

/// OCCT IntTools_FaceFace::CorrectPlaneBoundaries (L3126-3144) +
/// CorrectSurfaceBoundaries (L2050-2150): the FF domain of a face is not its
/// exact UV rectangle but a slightly enlarged one — a plane is expanded by 10%
/// of each parameter range, a quadric by the tolerance in the non-periodic
/// directions (periodic directions are clamped to the natural domain).  OCCT
/// uses this enlarged domain for BOTH the IntPatch boundary-vertex collection
/// (PutPointsOnLine walks the TopolTool boundary = the enlarged rectangle) and
/// the LineConstructor domain classification.
pub fn correct_ff_uv(surf: &Surface3, uv: [f64; 4], tol: f64) -> [f64; 4] {
    let mut out = uv;
    match surf {
        Surface3::Plane(_) => {
            // OCCT CorrectPlaneBoundaries: dU = 0.1 * (aUmax - aUmin).
            if out[0].is_finite() && out[1].is_finite() {
                let du = 0.1 * (out[1] - out[0]);
                out[0] -= du;
                out[1] += du;
            }
            if out[2].is_finite() && out[3].is_finite() {
                let dv = 0.1 * (out[3] - out[2]);
                out[2] -= dv;
                out[3] += dv;
            }
        }
        _ => {
            use rcad_kernel::geom::SurfaceEval;
            let d = surf.default_domain();
            let isuperiodic = surf.is_u_periodic();
            let isvperiodic = surf.is_v_periodic();
            // OCCT: enlarge for Bezier/BSpline/Extrusion/Revolution/Cylinder.
            let enlarge = matches!(
                surf,
                Surface3::Cylinder(_)
                    | Surface3::BSpline(_)
                    | Surface3::Bezier(_)
                    | Surface3::LinearExtrusion(_)
                    | Surface3::Revolution(_)
            );
            let snap = |cur: f64, lo: f64, hi: f64| {
                if cur.is_finite() && (cur - lo) > tol {
                    cur - tol
                } else {
                    lo
                }
            };
            let snap_hi = |cur: f64, lo: f64, hi: f64| {
                if cur.is_finite() && (hi - cur) > tol {
                    cur + tol
                } else {
                    hi
                }
            };
            if !isuperiodic && enlarge {
                out[0] = snap(out[0], d[0], d[1]);
                out[1] = snap_hi(out[1], d[0], d[1]);
            }
            if !isvperiodic && enlarge {
                out[2] = snap(out[2], d[2], d[3]);
                out[3] = snap_hi(out[3], d[2], d[3]);
            }
            // Periodic directions are clamped to the natural domain.
            if isuperiodic {
                out[0] = d[0];
                out[1] = d[1];
            }
            if isvperiodic {
                out[2] = d[2];
                out[3] = d[3];
            }
        }
    }
    out
}

#[derive(Debug, Clone)]
pub struct IntersectionCurve {
    pub curve: Curve3,
    pub t_range: [f64; 2],
    pub pcurve1: Option<Curve2d>,
    pub pcurve2: Option<Curve2d>,
    pub tolerance: f64,
    pub tang_tolerance: f64,
}

/// OCCT IntAna_ResultType (IntAna_QuadQuadGeo.hxx L31-42).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntAnaResultType {
    Empty,
    Line,
    Circle,
    Point,
    Same,
    Parallel,
    NoGeometricSolution,
}

/// OCCT IntAna_QuadQuadGeo result for the plane-plane case: type + line.
struct QuadQuadGeoPlnPln {
    typeres: IntAnaResultType,
    /// Intersection line (origin, direction) — valid for IntAna_Line.
    line_origin: DVec3,
    line_dir: DVec3,
}

/// OCCT IntAna_QuadQuadGeo::Perform(Plane, Plane, TolAng, Tol)
/// (IntAna_QuadQuadGeo.cxx L389-512).
fn quad_quad_geo_pln_pln(p1: &Plane, p2: &Plane, tol_ang: f64, tol: f64) -> QuadQuadGeoPlnPln {
    // P1.Coefficients(A1, B1, C1, D1); P2.Coefficients(A2, B2, C2, D2)
    let a1 = p1.normal.x;
    let b1 = p1.normal.y;
    let c1 = p1.normal.z;
    let d1 = -p1.normal.dot(p1.origin);
    let a2 = p2.normal.x;
    let b2 = p2.normal.y;
    let c2 = p2.normal.z;
    let d2 = -p2.normal.dot(p2.origin);

    let a_vn1 = DVec3::new(a1, b1, c1);
    let a_vn2 = DVec3::new(a2, b2, c2);
    let vd = a_vn1.cross(a_vn2);

    let a_loc_p1 = p1.origin;
    let a_loc_p2 = p2.origin;

    let dist1 = a2 * a_loc_p1.x + b2 * a_loc_p1.y + c2 * a_loc_p1.z + d2;
    let dist2 = a1 * a_loc_p2.x + b1 * a_loc_p2.y + c1 * a_loc_p2.z + d1;

    let a_mvd = vd.length();
    if a_mvd <= tol_ang {
        // normals are collinear - planes are same or parallel
        let typeres = if dist1.abs() <= tol && dist2.abs() <= tol {
            IntAnaResultType::Same
        } else {
            IntAnaResultType::Empty
        };
        return QuadQuadGeoPlnPln {
            typeres,
            line_origin: DVec3::ZERO,
            line_dir: DVec3::X,
        };
    }

    let a_eps = 1e-16;
    let mut denom = a1 * a2 + b1 * b2 + c1 * c2;
    let denom2 = denom * denom;
    let ddenom = 1.0 - denom2;
    denom = if ddenom.abs() <= a_eps { a_eps } else { ddenom };

    let par1 = dist1 / denom;
    let par2 = -dist2 / denom;

    let inter1 = a_vn1.cross(vd);
    let inter2 = a_vn2.cross(vd);

    let x1 = a_loc_p1.x + par1 * inter1.x;
    let y1 = a_loc_p1.y + par1 * inter1.y;
    let z1 = a_loc_p1.z + par1 * inter1.z;
    let x2 = a_loc_p2.x + par2 * inter2.x;
    let y2 = a_loc_p2.y + par2 * inter2.y;
    let z2 = a_loc_p2.z + par2 * inter2.z;

    let mut pt1 = DVec3::new((x1 + x2) * 0.5, (y1 + y2) * 0.5, (z1 + z2) * 0.5);
    let dir1 = vd.normalize();

    // OCCT L448-509: when the angle between the planes is small, the origin
    // of the intersection line is computed with error and should be refined.
    let a_tresh_ang = 2e-6; // 1.e-4 deg
    let a_tresh_dist = 1e-12;
    if a_mvd < a_tresh_ang {
        let a_dist1 = a1 * pt1.x + b1 * pt1.y + c1 * pt1.z + d1;
        let a_dist2 = a2 * pt1.x + b2 * pt1.y + c2 * pt1.z + d2;
        if a_dist1.abs() > a_tresh_dist || a_dist2.abs() > a_tresh_dist {
            // 1. IntAna_IntConicQuad(aL1, P1): line through pt1 along n1 × plane 1
            let a_dn1 = a_vn1.normalize();
            let a_l1_origin = pt1;
            let a_l1_dir = a_dn1;
            let (icq1_done, icq1_pnt) = int_conic_quad_line_pln(
                a_l1_origin, a_l1_dir, p1, tol_ang, tol);
            if !icq1_done {
                return QuadQuadGeoPlnPln {
                    typeres: IntAnaResultType::Empty,
                    line_origin: DVec3::ZERO,
                    line_dir: DVec3::X,
                };
            }
            // 2. IntAna_IntConicQuad(aL2, P2): line through the point along dir1 × n1
            let a_dl2 = dir1.cross(a_vn1).normalize();
            let a_l2_origin = icq1_pnt;
            let (icq2_done, icq2_pnt, icq2_parallel) = int_conic_quad_line_pln_par(
                a_l2_origin, a_dl2, p2, tol_ang, tol);
            if !icq2_done {
                return QuadQuadGeoPlnPln {
                    typeres: IntAnaResultType::Empty,
                    line_origin: DVec3::ZERO,
                    line_dir: DVec3::X,
                };
            }
            if icq2_parallel {
                return QuadQuadGeoPlnPln {
                    typeres: IntAnaResultType::Empty,
                    line_origin: DVec3::ZERO,
                    line_dir: DVec3::X,
                };
            }
            pt1 = icq2_pnt;
        }
    }

    QuadQuadGeoPlnPln {
        typeres: IntAnaResultType::Line,
        line_origin: pt1,
        line_dir: dir1,
    }
}

/// OCCT IntAna_IntConicQuad::Perform(Line, Plane, Tolang, Tol)
/// (IntAna_IntConicQuad.cxx L436-492). Returns (done, intersection point).
fn int_conic_quad_line_pln(
    orig: DVec3, dir: DVec3, p: &Plane, tol_ang: f64, tol: f64,
) -> (bool, DVec3) {
    let (done, pnt, _parallel) =
        int_conic_quad_line_pln_par(orig, dir, p, tol_ang, tol);
    (done, pnt)
}

/// OCCT IntAna_IntConicQuad::Perform(Line, Plane, Tolang, Tol) — full result.
fn int_conic_quad_line_pln_par(
    orig: DVec3, dir: DVec3, p: &Plane, tol_ang: f64, tol: f64,
) -> (bool, DVec3, bool) {
    let a = p.normal.x;
    let b = p.normal.y;
    let c = p.normal.z;
    let d = -p.normal.dot(p.origin);
    let (al, bl, cl) = (dir.x, dir.y, dir.z);

    let direc = a * al + b * bl + c * cl;
    let dis = a * orig.x + b * orig.y + c * orig.z + d;

    let mut parallel = false;
    if direc.abs() < tol_ang {
        parallel = true;
        // OCCT L464-475: when Len == 0 this block is skipped.
        // IntAna_QuadQuadGeo refinement calls with default Len (0), so skip.
        let _ = tol;
    }
    if parallel {
        // inquadric = |Dis| < Tolang (not used by the caller beyond done/point)
        return (true, orig, true);
    }
    // single intersection point
    let par = -dis / direc;
    (
        true,
        DVec3::new(orig.x + par * al, orig.y + par * bl, orig.z + par * cl),
        false,
    )
}

/// OCCT ElCLib::Parameter(gp_Lin2d, gp_Pnt2d) — parameter of a point on a 2D line.
fn elclib_lin2d_parameter(line: &Line2d, p: DVec2) -> f64 {
    (p - line.origin).dot(line.direction)
}

/// OCCT INTER (IntTools_FaceFace.cxx L2563-2569) — the 1:1 translation.
fn inter(d1: f64, d2: f64, tol: f64) -> bool {
    (d1 > tol && d2 < -tol)
        || (d1 < -tol && d2 > tol)
        || ((d1 <= tol && d1 >= -tol) && (d2 > tol || d2 < -tol))
        || ((d2 <= tol && d2 >= -tol) && (d1 > tol || d1 < -tol))
}

/// OCCT COINC (IntTools_FaceFace.cxx L2570-2571).
fn coinc(d1: f64, d2: f64, tol: f64) -> bool {
    (d1 >= -tol && d1 <= tol) && (d2 >= -tol && d2 <= tol)
}

/// OCCT ClassifyLin2d (IntTools_FaceFace.cxx L2575-2740): clip the 2D line
/// against the surface UV rectangle [xmin,xmax]x[ymin,ymax], returning the
/// parameter interval [P1,P2] where the line is inside. Returns false when
/// the line does not cross the rectangle (or the common part is degenerate).
fn classify_lin2d(
    uv: [f64; 4], line: &Line2d, the_tol: f64,
) -> (bool, f64, f64) {
    let (xmin, xmax, ymin, ymax) = (uv[0], uv[1], uv[2], uv[3]);
    let mut d1;
    let mut d2;

    // gp_Lin2d coefficients A*u + B*v + C = 0 from origin+unit dir:
    //   A = -dir.y, B = dir.x, C = dir.y*ox - dir.x*oy
    let a = -line.direction.y;
    let b = line.direction.x;
    let c = -(a * line.origin.x + b * line.origin.y);

    let mut par = [0.0f64; 2];
    let mut nbi = 0usize;

    // xmin, ymin <-> xmin, ymax
    d1 = a * xmin + b * ymin + c;
    d2 = a * xmin + b * ymax + c;
    if inter(d1, d2, the_tol) {
        // intersection with boundary
        let y = -(c + a * xmin) / b;
        par[nbi] = elclib_lin2d_parameter(line, DVec2::new(xmin, y));
        nbi += 1;
    } else if coinc(d1, d2, the_tol) {
        // coincidence with boundary
        par[0] = elclib_lin2d_parameter(line, DVec2::new(xmin, ymin));
        par[1] = elclib_lin2d_parameter(line, DVec2::new(xmin, ymax));
        nbi = 2;
    }

    if nbi == 2 {
        if (par[0] - par[1]).abs() > the_tol {
            let (p1, p2) = if par[0] < par[1] { (par[0], par[1]) } else { (par[1], par[0]) };
            return (true, p1, p2);
        } else {
            return (false, 0.0, 0.0);
        }
    }

    // xmin, ymax <-> xmax, ymax
    d1 = d2;
    d2 = a * xmax + b * ymax + c;
    if d1 > the_tol || d1 < -the_tol {
        // to avoid checking coincidence with the same point
        if inter(d1, d2, the_tol) {
            let x = -(c + b * ymax) / a;
            par[nbi] = elclib_lin2d_parameter(line, DVec2::new(x, ymax));
            nbi += 1;
        } else if coinc(d1, d2, the_tol) {
            par[0] = elclib_lin2d_parameter(line, DVec2::new(xmin, ymax));
            par[1] = elclib_lin2d_parameter(line, DVec2::new(xmax, ymax));
            nbi = 2;
        }
    }

    if nbi == 2 {
        if (par[0] - par[1]).abs() > the_tol {
            let (p1, p2) = if par[0] < par[1] { (par[0], par[1]) } else { (par[1], par[0]) };
            return (true, p1, p2);
        } else {
            return (false, 0.0, 0.0);
        }
    }

    // xmax, ymax <-> xmax, ymin
    d1 = d2;
    d2 = a * xmax + b * ymin + c;
    if d1 > the_tol || d1 < -the_tol {
        if inter(d1, d2, the_tol) {
            let y = -(c + a * xmax) / b;
            par[nbi] = elclib_lin2d_parameter(line, DVec2::new(xmax, y));
            nbi += 1;
        } else if coinc(d1, d2, the_tol) {
            par[0] = elclib_lin2d_parameter(line, DVec2::new(xmax, ymax));
            par[1] = elclib_lin2d_parameter(line, DVec2::new(xmax, ymin));
            nbi = 2;
        }
    }

    if nbi == 2 {
        if (par[0] - par[1]).abs() > the_tol {
            let (p1, p2) = if par[0] < par[1] { (par[0], par[1]) } else { (par[1], par[0]) };
            return (true, p1, p2);
        } else {
            return (false, 0.0, 0.0);
        }
    }

    // xmax, ymin <-> xmin, ymin
    d1 = d2;
    d2 = a * xmin + b * ymin + c;
    if d1 > the_tol || d1 < -the_tol {
        if inter(d1, d2, the_tol) {
            let x = -(c + b * ymin) / a;
            par[nbi] = elclib_lin2d_parameter(line, DVec2::new(x, ymin));
            nbi += 1;
        } else if coinc(d1, d2, the_tol) {
            par[0] = elclib_lin2d_parameter(line, DVec2::new(xmax, ymin));
            par[1] = elclib_lin2d_parameter(line, DVec2::new(xmin, ymin));
            nbi = 2;
        }
    }

    if nbi == 2 {
        if (par[0] - par[1]).abs() > the_tol {
            let (p1, p2) = if par[0] < par[1] { (par[0], par[1]) } else { (par[1], par[0]) };
            return (true, p1, p2);
        } else {
            return (false, 0.0, 0.0);
        }
    }

    (false, 0.0, 0.0)
}

/// OCCT IntTools_Tools::ComputeIntRange (IntTools_Tools.cxx L783-804).
fn compute_int_range(tol1: f64, tol2: f64, angle: f64) -> f64 {
    const PI: f64 = std::f64::consts::PI;
    const ANGULAR: f64 = 1e-12; // Precision::Angular()
    if (PI / 2.0 - angle).abs() < ANGULAR {
        tol2
    } else {
        let an_angle = if angle > PI / 2.0 { PI - angle } else { angle };
        let a1 = tol1 * (PI / 2.0 - an_angle).tan();
        let a2 = tol2 / an_angle.sin();
        a1 + a2
    }
}

/// OCCT PerformPlanes (IntTools_FaceFace.cxx L2427-2559): plane-plane
/// intersection producing a trimmed line + 2D pcurves on both planes.
fn perform_planes(
    p1: &Plane, uv1: [f64; 4],
    p2: &Plane, uv2: [f64; 4],
    tol_f1: f64, tol_f2: f64,
    tol_tang: f64,
) -> (bool, bool, Vec<IntersectionCurve>) {
    // IntAna_QuadQuadGeo aPlnInter(aPln1, aPln2, TolAng, TolTang)
    let tol_ang = 1e-8;
    let res = quad_quad_geo_pln_pln(p1, p2, tol_ang, tol_tang);

    if res.typeres == IntAnaResultType::Same {
        return (true, true, Vec::new()); // tangent faces
    }
    if res.typeres != IntAnaResultType::Line {
        return (true, false, Vec::new());
    }

    // Project the 3D intersection line onto each plane's (u,v) frame.
    // The intersection line lies in both planes, so the projected 2D
    // direction keeps unit length and the parameters match the 3D line.
    let proj_to_2d = |plane: &Plane| {
        let o = res.line_origin;
        let dd = res.line_dir;
        let o2 = DVec2::new(
            (o - plane.origin).dot(plane.u_dir),
            (o - plane.origin).dot(plane.v_dir),
        );
        let d2 = DVec2::new(dd.dot(plane.u_dir), dd.dot(plane.v_dir));
        Line2d::new(o2, d2)
    };
    let lin2d1 = proj_to_2d(p1);
    let lin2d2 = proj_to_2d(p2);

    // classify line2d1 relatively first plane
    let (crossed1, p11, p12) = classify_lin2d(uv1, &lin2d1, tol_tang);
    if !crossed1 {
        return (true, false, Vec::new());
    }
    // classify line2d2 relatively second plane
    let (crossed2, p21, p22) = classify_lin2d(uv2, &lin2d2, tol_tang);
    if !crossed2 {
        return (true, false, Vec::new());
    }

    // Analysis of parametric intervals: must have common part
    if p21 >= p12 {
        return (true, false, Vec::new());
    }
    if p22 <= p11 {
        return (true, false, Vec::new());
    }

    let pmin = p11.max(p21);
    let pmax = p12.min(p22);

    if pmax - pmin <= tol_tang {
        return (true, false, Vec::new());
    }

    let mut a_curve = IntersectionCurve {
        curve: Curve3::Line(Line3 {
            origin: res.line_origin,
            direction: res.line_dir,
        }),
        t_range: [pmin, pmax],
        pcurve1: Some(Curve2d::Line(lin2d1)),
        pcurve2: Some(Curve2d::Line(lin2d2)),
        tolerance: tol_f1.max(tol_f2),
        tang_tolerance: 0.0,
    };

    // Computation of the tangential tolerance
    let a_d1 = p1.normal;
    let a_d2 = p2.normal;
    let an_angle = a_d1.angle_between(a_d2);
    let a_dt = compute_int_range(tol_f1, tol_f2, an_angle);
    let a_tang_tol = (a_dt * a_dt + tol_f1 * tol_f1).sqrt();
    a_curve.tang_tolerance = a_tang_tol;

    (true, false, vec![a_curve])
}

pub struct FaceFace {
    surf1: Surface3,
    surf2: Surface3,
    uv1: [f64; 4],
    uv2: [f64; 4],
    tol1: f64,
    tol2: f64,
    fuzzy: f64,
    curves: Vec<IntersectionCurve>,
    /// Raw IntPatch lines before MakeCurve domain clipping (OCCT myIntersector lines).
    lines: Vec<IntPatchLine>,
    /// Face indices in the DS (used for FClass2d domain classification in MakeCurve).
    face1: usize,
    face2: usize,
    done: bool,
    tangent_faces: bool,
}

impl FaceFace {
    pub fn new() -> Self {
        FaceFace {
            surf1: Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            }),
            surf2: Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            }),
            uv1: [0.0, 0.0, 0.0, 0.0],
            uv2: [0.0, 0.0, 0.0, 0.0],
            tol1: 1e-7,
            tol2: 1e-7,
            fuzzy: rcad_kernel::CONFUSION,
            curves: Vec::new(),
            lines: Vec::new(),
            face1: 0,
            face2: 0,
            done: false,
            tangent_faces: false,
        }
    }

    pub fn set_surfaces(&mut self, s1: Surface3, s2: Surface3) {
        self.surf1 = s1;
        self.surf2 = s2;
    }

    pub fn set_uv_bounds(&mut self, uv1: [f64; 4], uv2: [f64; 4]) {
        self.uv1 = uv1;
        self.uv2 = uv2;
    }

    pub fn set_tolerances(&mut self, t1: f64, t2: f64) {
        self.tol1 = t1.max(1e-7);
        self.tol2 = t2.max(1e-7);
    }

    pub fn set_fuzzy_value(&mut self, fuzz: f64) {
        self.fuzzy = fuzz.max(rcad_kernel::CONFUSION);
    }

    /// Face indices in the DS, used for FClass2d domain classification in MakeCurve.
    pub fn set_face_indices(&mut self, f1: usize, f2: usize) {
        self.face1 = f1;
        self.face2 = f2;
    }

    pub fn is_done(&self) -> bool {
        self.done
    }
    pub fn has_intersection(&self) -> bool {
        !self.curves.is_empty()
    }
    pub fn tangent_faces(&self) -> bool {
        self.tangent_faces
    }
    pub fn make_curves(&self) -> Vec<IntersectionCurve> {
        self.curves.clone()
    }
    pub fn points(&self) -> Vec<crate::bop::int_tools::pnt_on_2_faces::PntOn2Faces> {
        Vec::new()
    }

    /// OCCT IntTools_FaceFace::Perform — compute intersection.
    pub fn perform(&mut self, ds: &crate::bop::ds::DS) {
        self.curves.clear();
        self.lines.clear();
        self.tangent_faces = false;
        let s1 = self.surf1.clone();
        let s2 = self.surf2.clone();
        let a_fuzz = self.fuzzy / 2.0;
        let tol_f1 = self.tol1 + a_fuzz;
        let tol_f2 = self.tol2 + a_fuzz;
        let tol = tol_f1 + tol_f2;
        let tol_tang = tol;

        // OCCT loads the surface adaptors with corrected UV bounds ONLY for the
        // non-plane-plane branches (plane-plane keeps the raw UV rectangle):
        //   plane x quadric / quadric x plane: CorrectPlaneBoundaries on the
        //     plane (10% expansion), CorrectSurfaceBoundaries on the quadric.
        //   quadric x quadric: CorrectSurfaceBoundaries on both.
        let is_pln_pln = matches!((&s1, &s2), (Surface3::Plane(_), Surface3::Plane(_)));
        let (uv1, uv2) = if is_pln_pln {
            (self.uv1, self.uv2)
        } else {
            let c1 = correct_ff_uv(&s1, self.uv1, tol * 2.0);
            let c2 = correct_ff_uv(&s2, self.uv2, tol * 2.0);
            (c1, c2)
        };

        match (&s1, &s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                let (done, tangent, curves) =
                    perform_planes(p1, uv1, p2, uv2, tol_f1, tol_f2, tol_tang);
                if !done {
                    self.done = false;
                    return;
                }
                self.tangent_faces = tangent;
                self.curves = curves;
                self.done = true;
            }
            _ => {
                // OCCT IntTools_FaceFace::Perform L441-474: for the remaining
                // analytic pairs route through IntPatch_Intersection
                // (-> IntPatch_ImpImpIntersection -> IntAna_QuadQuadGeo).
                self.intersect_int_patch(&s1, &s2, uv1, uv2, tol);
                // OCCT IntTools_FaceFace::MakeCurve (L695-1846): clip the raw
                // analytic lines to the two faces' domains.
                self.curves = face_make_curve::make_curves(
                    ds,
                    self.face1,
                    self.face2,
                    &s1,
                    uv1,
                    &s2,
                    uv2,
                    tol,
                    &self.lines,
                );
                self.done = true;
            }
        }
    }

    /// OCCT IntTools_FaceFace::Perform L523-531: myIntersector.Perform(...).
    ///
    /// Runs IntPatch_Intersection on the two surfaces and converts the
    /// produced IntPatch lines into rcad IntersectionCurves. The lines are
    /// the full analytic intersection curves (untrimmed); domain clipping
    /// happens in MakeCurve (IntTools_FaceFace.cxx L695-1846), which is the
    /// PaveFiller's responsibility downstream.
    fn intersect_int_patch(
        &mut self,
        s1: &Surface3,
        s2: &Surface3,
        uv1: [f64; 4],
        uv2: [f64; 4],
        tol: f64,
    ) {
        let mut inter =
            crate::geomalgo::int_patch::IntPatchIntersection::new();
        inter.perform(s1, s2, uv1, uv2, tol, tol);
        if inter.tangent_faces() {
            self.tangent_faces = true;
            self.curves.clear();
            self.lines.clear();
            return;
        }
        if inter.is_empty() {
            self.curves.clear();
            self.lines.clear();
            return;
        }
        // Keep the raw IntPatch lines; MakeCurve (face_make_curve::make_curves)
        // clips them to the face UV domains.
        self.lines = inter.sequence_of_line().to_vec();
    }
}

impl Default for FaceFace {
    fn default() -> Self {
        Self::new()
    }
}
