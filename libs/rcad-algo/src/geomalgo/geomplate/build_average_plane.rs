//! OCCT GeomPlate_BuildAveragePlane (TKGeomAlgo/GeomPlate) — 1:1 port of
//! GeomPlate_BuildAveragePlane.cxx (whole file L46-750) for the point path.
//!
//! Mapping: gp_Pnt/gp_Vec -> DVec3, Geom_Plane -> rcad Plane, Geom_Line ->
//! rcad Line3, NCollection_HArray1/Sequence -> Vec, math_Jacobi ->
//! rcad MathJacobi, GeomLib::AxeOfInertia -> rcad base::geom_lib.
//! Architecture notes:
//! - ctor2 leaves myTol/myNbBoundPoints uninitialized in OCCT; Rust requires
//!   an initializer, so they are zeroed (never read on the ctor2 path).
//! - OCCT NCollection_Array1 is 1-based; the Vec helpers below index with an
//!   explicit `- 1` offset inside `1..=n` loops to preserve OCCT form.

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::FRAC_PI_3;
use std::f64::consts::PI;

use glam::DVec3;

use rcad_kernel::base::geom_lib;
use rcad_kernel::geom::{Line3, Plane};
use rcad_kernel::math::gp::Ax2;
use rcad_kernel::math::math_jacobi::MathJacobi;
use rcad_kernel::math::MatD;

// OCCT Precision::SquareConfusion() = Confusion * Confusion.
const SQUARE_CONFUSION: f64 = 1e-14;
// OCCT gp::Resolution() = gp::RealSmall() (gp.hxx L60), the smallest
// normalized positive double.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

// ============================================================================
// GeomPlate_Aij
// ============================================================================

/// OCCT GeomPlate_Aij (GeomPlate_Aij.hxx L27-48) — a couple of normal indices
/// with their cross vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aij {
    pub ind1: usize,
    pub ind2: usize,
    pub vec: DVec3,
}

impl Aij {
    /// OCCT GeomPlate_Aij(anInd1, anInd2, aVec).
    pub fn new(ind1: usize, ind2: usize, vec: DVec3) -> Self {
        Aij { ind1, ind2, vec }
    }
}

// ============================================================================
// gp / ElSLib helpers (private)
// ============================================================================

/// OCCT gp_Trsf::SetRotation (gp_Trsf.cxx L90-99) with the axis through the
/// origin, applied to a vector — gp_Mat::SetRotation (gp_Mat.cxx L122-158).
/// This is gp_Vec::Rotated(Ax1(gp_Pnt(0,0,0), Axis), Ang) (gp_Vec.cxx).
fn gp_vec_rotated(axis_dir: DVec3, ang: f64, v: DVec3) -> DVec3 {
    let a_v = axis_dir / axis_dir.length();
    let a = a_v.x;
    let b = a_v.y;
    let c = a_v.z;

    let a_cos = ang.cos();
    let a_sin = ang.sin();
    let a_om_cos = 1.0 - a_cos;

    let a2 = a * a;
    let b2 = b * b;
    let c2 = c * c;
    let ab = a * b;
    let ac = a * c;
    let bc = b * c;

    // gp_Mat rows (gp_Mat.cxx L148-157).
    let m00 = 1.0 + a_om_cos * (-(b2 + c2));
    let m01 = a_om_cos * ab - a_sin * c;
    let m02 = a_om_cos * ac + a_sin * b;
    let m10 = a_om_cos * ab + a_sin * c;
    let m11 = 1.0 + a_om_cos * (-(a2 + c2));
    let m12 = a_om_cos * bc - a_sin * a;
    let m20 = a_om_cos * ac - a_sin * b;
    let m21 = a_om_cos * bc + a_sin * a;
    let m22 = 1.0 + a_om_cos * (-(a2 + b2));

    DVec3::new(
        m00 * v.x + m01 * v.y + m02 * v.z,
        m10 * v.x + m11 * v.y + m12 * v.z,
        m20 * v.x + m21 * v.y + m22 * v.z,
    )
}

/// OCCT gp_Vec::Angle (gp_Vec.hxx L488-495) -> gp_Dir::Angle (gp_Dir.cxx
/// L27-53): acos between 45 and 135 degrees, asin elsewhere.
fn gp_vec_angle(a: DVec3, b: DVec3) -> f64 {
    assert!(
        a.length() > GP_RESOLUTION && b.length() > GP_RESOLUTION,
        "gp_VectorWithNullMagnitude"
    );
    // gp_Vec::Angle converts both operands to gp_Dir (normalized).
    let da = a.normalize_or_zero();
    let db = b.normalize_or_zero();
    let cosinus = da.dot(db);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = da.cross(db).length();
        if cosinus < 0.0 {
            PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT gp_Vec::IsOpposite (gp_Vec.hxx L130-134):
/// PI - Angle(theOther) <= theAngularTolerance.
fn gp_vec_is_opposite(a: DVec3, b: DVec3, angular_tolerance: f64) -> bool {
    let an_ang = PI - gp_vec_angle(a, b);
    an_ang <= angular_tolerance
}

/// OCCT ElSLib::Parameters(Pln, P, U, V) (ElSLib.hxx L507) dispatching to
/// ElSLib::PlaneParameters (ElSLib.cxx L1547-1555): the plane frame is
/// orthonormal, so the local X/Y coordinates equal the dot products of
/// (P - Location) with XDirection / YDirection.
fn elslib_plane_parameters(plane: &Plane, p: DVec3) -> (f64, f64) {
    let vp = p - plane.origin;
    (vp.dot(plane.u_dir), vp.dot(plane.v_dir))
}

// ============================================================================
// GeomPlate_BuildAveragePlane
// ============================================================================

/// Computes an average inertial plane from an array of points.
///
/// OCCT: `GeomPlate_BuildAveragePlane`.
#[derive(Debug, Clone)]
pub struct BuildAveragePlane {
    my_pts: Vec<DVec3>,
    my_umax: f64,
    my_vmax: f64,
    my_vmin: f64,
    my_umin: f64,
    my_plane: Option<Plane>,
    my_tol: f64,
    my_line: Option<Line3>,
    my_ox: DVec3,
    my_oy: DVec3,
    my_g: DVec3,
    my_nb_bound_points: usize,
}

impl BuildAveragePlane {
    /// OCCT ctor 1 (GeomPlate_BuildAveragePlane.cxx L46-101).
    ///
    /// Tol differentiates the plane result from the line result.
    /// POption = 1: automatic parametrisation; POption = 2: parametrisation
    /// by eigen vectors.  NOption = 1: the average plane is the inertial
    /// plane; NOption = 2: the average plane is the plane of max. flux.
    pub fn new(
        pts: Vec<DVec3>,
        nb_bound_points: usize,
        tol: f64,
        p_option: i32,
        n_option: i32,
    ) -> Self {
        let mut bap = BuildAveragePlane {
            my_pts: pts,
            my_umax: 0.0,
            my_vmax: 0.0,
            my_vmin: 0.0,
            my_umin: 0.0,
            my_plane: None,
            my_tol: tol,
            my_line: None,
            my_ox: DVec3::ZERO,
            my_oy: DVec3::ZERO,
            my_g: DVec3::ZERO,
            my_nb_bound_points: nb_bound_points,
        };

        let oz = bap.def_plan(n_option);

        if oz.length_squared() > 0.0 {
            if p_option == 1 {
                // myPlane = new Geom_Plane(myG, OZ).
                let plane = Plane::new(bap.my_g, oz);
                bap.my_ox = plane.u_dir;
                bap.my_oy = plane.v_dir;
                bap.my_plane = Some(plane);
            } else {
                bap.base_plan(oz);
                // gp_Dir NDir(myOX ^ myOY); gp_Dir UDir(myOX);
                // gp_Ax3 triedre(myG, NDir, UDir); Geom_Plane(triedre).
                let n_dir = bap.my_ox.cross(bap.my_oy).normalize_or_zero();
                let u_dir = bap.my_ox.normalize_or_zero();
                let plane = Plane::with_axes(bap.my_g, n_dir, u_dir);
                bap.my_plane = Some(plane);
            }
            let pln = bap.my_plane.as_ref().unwrap();
            // ElSLib::Parameters(P, myG, myUmax, myVmax) then scan all points.
            let (mut umax, mut vmax) = elslib_plane_parameters(pln, bap.my_g);
            let mut umin = umax;
            let mut vmin = vmax;
            for i in 1..=bap.my_pts.len() {
                let (u, v) = elslib_plane_parameters(pln, bap.my_pts[i - 1]);
                if umax < u {
                    umax = u;
                }
                if umin > u {
                    umin = u;
                }
                if vmax < v {
                    vmax = v;
                }
                if vmin > v {
                    vmin = v;
                }
            }
            bap.my_umax = umax;
            bap.my_umin = umin;
            bap.my_vmax = vmax;
            bap.my_vmin = vmin;
        }

        if bap.is_line() {
            // myLine = new Geom_Line(myG, myOX).
            bap.my_line = Some(Line3::new(bap.my_g, bap.my_ox));
        }
        bap
    }

    /// OCCT ctor 2 (GeomPlate_BuildAveragePlane.cxx L109-258) — creates the
    /// plane from the "best vector" of a normal sequence.
    pub fn from_normals(normals: &[DVec3], pts: Vec<DVec3>) -> Self {
        let mut bap = BuildAveragePlane {
            my_pts: pts,
            my_umax: 0.0,
            my_vmax: 0.0,
            my_vmin: 0.0,
            my_umin: 0.0,
            my_plane: None,
            // OCCT leaves myTol/myNbBoundPoints uninitialized on this path;
            // zeroed (never read here).
            my_tol: 0.0,
            my_line: None,
            my_ox: DVec3::ZERO,
            my_oy: DVec3::ZERO,
            my_g: DVec3::ZERO,
            my_nb_bound_points: 0,
        };

        let mut best_vec;
        let nn = normals.len();

        if nn == 1 {
            best_vec = normals[0];
        } else if nn == 2 {
            best_vec = normals[0] + normals[1];
            let a_sq_magn = best_vec.length_squared();
            if a_sq_magn < SQUARE_CONFUSION {
                let a_sq1 = normals[0].length_squared();
                let a_sq2 = normals[1].length_squared();
                if a_sq1 > a_sq2 {
                    best_vec = normals[0].normalize_or_zero();
                } else {
                    best_vec = normals[1].normalize_or_zero();
                }
            } else {
                best_vec /= a_sq_magn.sqrt();
            }
        } else {
            // The common case (L152-234).
            let mut max_angle = 0.0f64;
            for i in 1..=nn - 1 {
                for j in (i + 1)..=nn {
                    let angle = gp_vec_angle(normals[i - 1], normals[j - 1]);
                    if angle > max_angle {
                        max_angle = angle;
                    }
                }
            }
            max_angle *= 1.2;
            max_angle /= 2.0;
            let n_int = 50;

            // NCollection_Array1 OptVec/OptScal(1, NN*(NN-1)/2) — 1-based,
            // stored with a -1 offset.
            let n_pairs = nn * (nn - 1) / 2;
            let mut opt_vec = vec![DVec3::ZERO; n_pairs];
            let mut opt_scal = vec![0.0f64; n_pairs];

            let mut k = 1usize;
            for i in 1..=nn - 1 {
                for j in (i + 1)..=nn {
                    opt_scal[k - 1] = f64::MIN;

                    let step = max_angle / n_int as f64;
                    let mut vec = normals[i - 1] + normals[j - 1];

                    let a_sq_magn = vec.length_squared();
                    if a_sq_magn < SQUARE_CONFUSION {
                        // OCCT `continue` still runs `j++, k++`.
                        k += 1;
                        continue;
                    }

                    vec /= a_sq_magn.sqrt();

                    let cross1 = normals[i - 1].cross(normals[j - 1]);
                    let cross2 = vec.cross(cross1);
                    // gp_Ax1 Axe(gp_Pnt(0, 0, 0), Cross2) — direction is
                    // normalized by gp_Ax1.
                    let axe_dir = cross2.normalize_or_zero();

                    let mut vec1 = gp_vec_rotated(axe_dir, -max_angle, vec);

                    for _n in 0..=(2 * n_int) {
                        vec1 = gp_vec_rotated(axe_dir, step, vec1);
                        let mut min_scal = f64::MAX;
                        for m in 1..=nn {
                            let scal = vec1.dot(normals[m - 1]);
                            if scal < min_scal {
                                min_scal = scal;
                            }
                        }
                        if min_scal > opt_scal[k - 1] {
                            opt_scal[k - 1] = min_scal;
                            opt_vec[k - 1] = vec1;
                        }
                    }
                    k += 1;
                } // for i, for j
            }
            // Find maximum among all maximums.
            let mut best_scal = f64::MIN;
            let mut index = 0usize;
            for k in 1..=opt_scal.len() {
                if opt_scal[k - 1] > best_scal {
                    best_scal = opt_scal[k - 1];
                    index = k;
                }
            }
            best_vec = opt_vec[index - 1];
        }

        // Making the plane myPlane (L237-258).
        // OCCT `gp_Ax2 Axe;` default is the OXYZ frame; AxeOfInertia
        // overwrites it below.
        let mut axe = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
        let mut is_singular = false;
        // OCCT copies myPts into a local Array1 for the out-parameter call;
        // the rcad port takes a const slice.
        geom_lib::axe_of_inertia(&bap.my_pts, &mut axe, &mut is_singular, 1.0e-7);
        let best_dir = best_vec.normalize_or_zero();
        let mut x_dir = best_dir.cross(axe.x_direction).normalize_or_zero();
        x_dir = x_dir.cross(best_dir).normalize_or_zero();

        let plane = Plane::with_axes(axe.location, best_dir, x_dir);

        // Initializing myUmin, myVmin, myUmax, myVmax (L247-277).
        let (mut umax, mut vmax) = elslib_plane_parameters(&plane, axe.location);
        let mut umin = umax;
        let mut vmin = vmax;
        for i in 1..=bap.my_pts.len() {
            // gp_Vec aVec(Pln.Location(), myPts(i))  ==  myPts(i) - Location.
            let a_vec = bap.my_pts[i - 1] - plane.origin;
            let norm_vec = plane.normal;
            // (aVec * NormVec) * NormVec  ==  (aVec . NormVec) * NormVec.
            let norm_vec = norm_vec * a_vec.dot(norm_vec);

            let (u, v) = elslib_plane_parameters(&plane, bap.my_pts[i - 1] - norm_vec);
            if u > umax {
                umax = u;
            }
            if u < umin {
                umin = u;
            }
            if v > vmax {
                vmax = v;
            }
            if v < vmin {
                vmin = v;
            }
        }
        bap.my_umax = umax;
        bap.my_umin = umin;
        bap.my_vmax = vmax;
        bap.my_vmin = vmin;
        // Initializing myOX, myOY.
        bap.my_ox = plane.u_dir;
        bap.my_oy = plane.v_dir;
        bap.my_plane = Some(plane);
        bap
    }

    /// OCCT Plane() (L282-289) — returns the average plane.
    ///
    /// Raises (panics) if the object is a line.
    pub fn plane(&self) -> Option<&Plane> {
        assert!(
            !self.is_line(),
            "Cannot use the function 'GeomPlate_BuildAveragePlane::Plane()', the Object is a 'Geom_Line'"
        );
        self.my_plane.as_ref()
    }

    /// OCCT Line() (L491-496) — returns the line when 2 eigenvalues are null.
    ///
    /// Raises (panics) if the object is a plane.
    pub fn line(&self) -> Option<&Line3> {
        assert!(
            !self.is_plane(),
            "Cannot use the function 'GeomPlate_BuildAveragePlane::Line()', the Object is a 'Geom_Plane'"
        );
        self.my_line.as_ref()
    }

    /// OCCT IsPlane() (L501-504).
    pub fn is_plane(&self) -> bool {
        let oz = self.my_ox.cross(self.my_oy);
        oz.length_squared() != 0.0
    }

    /// OCCT IsLine() (L509-512).
    pub fn is_line(&self) -> bool {
        let oz = self.my_ox.cross(self.my_oy);
        oz.length_squared() == 0.0
    }

    /// OCCT MinMaxBox() (L292-301) — the minimal box including all normal
    /// projections of the initial array on the plane.
    pub fn min_max_box(&self, umin: &mut f64, umax: &mut f64, vmin: &mut f64, vmax: &mut f64) {
        *umax = self.my_umax;
        *umin = self.my_umin;
        *vmax = self.my_vmax;
        *vmin = self.my_vmin;
    }

    /// OCCT DefPlan (L305-382) — defines the average plane direction.
    /// NOption = 1: the inertial plane; NOption = 2: the plane of max. flux.
    fn def_plan(&mut self, n_option: i32) -> DVec3 {
        let mut oz = DVec3::ZERO;
        // Barycenter myG (L311-323).
        let mut gb = DVec3::ZERO;
        let nb = self.my_pts.len();
        for point in &self.my_pts {
            gb += *point;
        }
        self.my_g = gb / nb as f64;

        if n_option == 1 {
            let mut axe = Ax2::new(DVec3::ZERO, DVec3::Z, DVec3::X);
            let mut is_singular = false;
            geom_lib::axe_of_inertia(&self.my_pts, &mut axe, &mut is_singular, self.my_tol);

            self.my_ox = axe.x_direction;
            self.my_oy = axe.y_direction;

            oz = axe.direction;

            if self.my_nb_bound_points != 0 && self.my_pts.len() != self.my_nb_bound_points {
                // Boundary cross-product consistency check (L340-367).
                let mut a = DVec3::ZERO;
                for i in 3..=self.my_nb_bound_points {
                    let b = self.my_pts[i - 2] - self.my_pts[0];
                    let c = self.my_pts[i - 1] - self.my_pts[0];
                    let d = b.cross(c);
                    a += d;
                }
                let oz1 = a;
                let mut the_angle = gp_vec_angle(oz, oz1);
                if the_angle > FRAC_PI_2 {
                    the_angle = PI - the_angle;
                }
                if the_angle > FRAC_PI_3 {
                    oz = oz1;
                }
            }
        } else if n_option == 2 {
            // Pure boundary cross-product sum (L369-381).
            let mut a = DVec3::ZERO;
            for i in 3..=self.my_nb_bound_points {
                let b = self.my_pts[i - 2] - self.my_pts[0];
                let c = self.my_pts[i - 1] - self.my_pts[0];
                let d = b.cross(c);
                a += d;
            }
            oz = a;
        }
        oz
    }

    /// OCCT BasePlan (L386-467) — computes a base of the average plane
    /// defined by (myG, N) using eigen vectors.
    fn base_plan(&mut self, oz: DVec3) {
        let mut m = MatD::new(3, 3);
        let nb = self.my_pts.len();

        for i in 1..=nb {
            let point = self.my_pts[i - 1];
            // Proj = P - myG projected onto the OZ plane.
            let mut proj = point - self.my_g;
            let mut scal = proj.dot(oz);
            scal /= oz.dot(oz);
            proj -= oz * scal;
            m.set(1, 1, m.get(1, 1) + proj.x * proj.x);
            m.set(2, 2, m.get(2, 2) + proj.y * proj.y);
            m.set(3, 3, m.get(3, 3) + proj.z * proj.z);
            m.set(1, 2, m.get(1, 2) + proj.x * proj.y);
            m.set(1, 3, m.get(1, 3) + proj.x * proj.z);
            m.set(2, 3, m.get(2, 3) + proj.y * proj.z);
        }
        m.set(2, 1, m.get(1, 2));
        m.set(3, 1, m.get(1, 3));
        m.set(3, 2, m.get(2, 3));

        let jacobi = MathJacobi::new(&m);
        let n1 = jacobi.value(1);
        let n2 = jacobi.value(2);
        let n3 = jacobi.value(3);

        let r1 = n1.min(n2.min(n3));
        let r2: f64;
        let m1: usize;
        let m2: usize;
        let m3: usize;
        if r1 == n1 {
            m1 = 1;
            r2 = n2.min(n3);
            if r2 == n2 {
                m2 = 2;
                m3 = 3;
            } else {
                m2 = 3;
                m3 = 2;
            }
        } else if r1 == n2 {
            m1 = 2;
            r2 = n1.min(n3);
            if r2 == n1 {
                m2 = 1;
                m3 = 3;
            } else {
                m2 = 3;
                m3 = 1;
            }
        } else {
            m1 = 3;
            r2 = n1.min(n2);
            if r2 == n1 {
                m2 = 1;
                m3 = 2;
            } else {
                m2 = 2;
                m3 = 1;
            }
        }
        // J.Vector(m1, V1) — V1 is fetched but never used by OCCT either.
        let _v1 = jacobi.vector(m1);
        let v2 = jacobi.vector(m2);
        let v3 = jacobi.vector(m3);

        if ((n1.abs() <= self.my_tol) && (n2.abs() <= self.my_tol))
            || ((n2.abs() <= self.my_tol) && (n3.abs() <= self.my_tol))
            || ((n1.abs() <= self.my_tol) && (n3.abs() <= self.my_tol))
        {
            self.my_ox = DVec3::new(v3.get(1), v3.get(2), v3.get(3));
            self.my_oy = oz.cross(self.my_ox);
        } else {
            self.my_ox = DVec3::new(v3.get(1), v3.get(2), v3.get(3));
            self.my_oy = DVec3::new(v2.get(1), v2.get(2), v2.get(3));
        }
    }

    /// OCCT HalfSpace (L515-750) — static half-space normal filtering.
    pub fn half_space(
        new_normals: &[DVec3],
        normals: &mut Vec<DVec3>,
        bset: &mut Vec<Aij>,
        lin_tol: f64,
        ang_tol: f64,
    ) -> bool {
        let square_tol = lin_tol * lin_tol;

        // 1 (L530-534)
        let save_normals = normals.clone();
        let save_bset = bset.clone();

        let null_vec = DVec3::ZERO;
        let mut b1set: Vec<Aij> = Vec::new();
        let mut b2set: Vec<Aij> = Vec::new();

        let mut i = 1usize;
        if normals.is_empty() {
            if new_normals.len() == 1 {
                normals.push(*new_normals.last().unwrap());
                return true;
            }
            // 2 (L545-558)
            let mut cross = new_normals[0].cross(new_normals[1]);
            if cross.length_squared() <= square_tol {
                return false;
            }

            cross = cross.normalize_or_zero();
            bset.push(Aij::new(1, 2, cross));
            bset.push(Aij::new(2, 1, -cross));
            normals.push(new_normals[0]);
            normals.push(new_normals[1]);
            i = 3;
        }

        while i <= new_normals.len() {
            // 3 (L562-571)
            for j in 1..=bset.len() {
                let scal = bset[j - 1].vec.dot(new_normals[i - 1]);
                if scal >= -lin_tol {
                    b2set.push(bset[j - 1]);
                }
            }

            let ii = normals.len() + 1;
            for j in 1..=ii - 1 {
                if normals[j - 1].length_squared() == 0.0 {
                    continue;
                }
                // 4 (L575-637)
                let mut cross = new_normals[i - 1].cross(normals[j - 1]);
                if cross.length_squared() <= square_tol {
                    *normals = save_normals.clone();
                    *bset = save_bset.clone();
                    return false;
                }
                cross = cross.normalize_or_zero();
                let mut is_new = true;
                for k in 1..=b2set.len() {
                    if gp_vec_is_opposite(b2set[k - 1].vec, -cross, ang_tol) {
                        let ind1 = b2set[k - 1].ind1;
                        let ind2 = b2set[k - 1].ind2;
                        if ind1 == ii || ind2 == ii {
                            is_new = false;
                            break;
                        }
                        let cross1 = normals[ind1 - 1].cross(new_normals[i - 1]);
                        let cross2 = normals[ind2 - 1].cross(new_normals[i - 1]);
                        if cross1.length_squared() <= square_tol
                            || cross2.length_squared() <= square_tol
                        {
                            *normals = save_normals.clone();
                            *bset = save_bset.clone();
                            return false;
                        }
                        if gp_vec_is_opposite(cross1, cross2, ang_tol) {
                            let cross2 = normals[ind1 - 1].cross(normals[ind2 - 1]);
                            if gp_vec_is_opposite(cross1, cross2, ang_tol) {
                                *normals = save_normals.clone();
                                *bset = save_bset.clone();
                                return false;
                            }
                        } else if gp_vec_angle(new_normals[i - 1], normals[ind1 - 1])
                            > gp_vec_angle(new_normals[i - 1], normals[ind2 - 1])
                        {
                            b2set[k - 1].ind2 = ind1;
                            b2set[k - 1].ind1 = ii;
                        } else {
                            b2set[k - 1].ind1 = ii;
                        }
                        is_new = false;
                        break;
                    }
                }
                if is_new {
                    b1set.push(Aij::new(ii, j, cross));
                }

                // Cross.Reverse();
                cross = -cross;
                let mut is_new = true;
                for k in 1..=b2set.len() {
                    if gp_vec_is_opposite(b2set[k - 1].vec, -cross, ang_tol) {
                        let ind1 = b2set[k - 1].ind1;
                        let ind2 = b2set[k - 1].ind2;
                        if ind1 == ii || ind2 == ii {
                            is_new = false;
                            break;
                        }
                        let cross1 = normals[ind1 - 1].cross(new_normals[i - 1]);
                        let cross2 = normals[ind2 - 1].cross(new_normals[i - 1]);
                        if cross1.length_squared() <= square_tol
                            || cross2.length_squared() <= square_tol
                        {
                            *normals = save_normals.clone();
                            *bset = save_bset.clone();
                            return false;
                        }
                        if gp_vec_is_opposite(cross1, cross2, ang_tol) {
                            let cross2 = normals[ind1 - 1].cross(normals[ind2 - 1]);
                            if gp_vec_is_opposite(cross1, cross2, ang_tol) {
                                *normals = save_normals.clone();
                                *bset = save_bset.clone();
                                return false;
                            }
                        } else if gp_vec_angle(new_normals[i - 1], normals[ind1 - 1])
                            > gp_vec_angle(new_normals[i - 1], normals[ind2 - 1])
                        {
                            b2set[k - 1].ind2 = ind1;
                            b2set[k - 1].ind1 = ii;
                        } else {
                            b2set[k - 1].ind1 = ii;
                        }
                        is_new = false;
                        break;
                    }
                }
                if is_new {
                    b1set.push(Aij::new(ii, j, cross));
                }
            }

            // 5 (L640-656)
            for j in 1..=b1set.len() {
                let mut is_ge_null = true;
                for k in 1..=normals.len() {
                    if normals[k - 1].length_squared() == 0.0 {
                        continue;
                    }
                    if b1set[j - 1].vec.dot(normals[k - 1]) < -lin_tol {
                        is_ge_null = false;
                        break;
                    }
                }
                if is_ge_null {
                    b2set.push(b1set[j - 1]);
                }
            }

            // 6 (L659-666)
            if b2set.is_empty() {
                *normals = save_normals.clone();
                *bset = save_bset.clone();
                return false;
            }

            // 7 (L669-673)
            *bset = b2set.clone();
            b2set.clear();
            b1set.clear();
            normals.push(new_normals[i - 1]);

            // 8 (L676-692)
            for j in 1..=normals.len() {
                if normals[j - 1].length_squared() == 0.0 {
                    continue;
                }
                let mut is_found = false;
                for k in 1..=bset.len() {
                    if j == bset[k - 1].ind1 || j == bset[k - 1].ind2 {
                        is_found = true;
                        break;
                    }
                }
                if !is_found {
                    normals[j - 1] = null_vec;
                }
            }

            i += 1;
        }

        true
    }
}
