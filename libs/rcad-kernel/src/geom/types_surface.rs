//! Surface payload types: plane and the quadrics, helicoid, pipe, BSpline /
//! Bezier / triangular Bezier, offset, trimmed, extrusion, revolution, ruled
//! and Coons surfaces, plus the [`Surface3`] and [`PrimitiveSolid`] enums and
//! the `bspline_is_planar` / `bspline_to_plane` helpers.
//!
//! Extracted verbatim from `geom/mod.rs` (project Rule 5 file-size split);
//! all public paths stay stable through the `pub use` re-exports in `mod.rs`.
use glam::{DVec2, DVec3};
use serde::{Deserialize, Serialize};

use super::{Curve3, Point3, Vec3, any_perpendicular};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Plane {
    pub origin: Point3,
    /// Unit normal vector (OCCT gp_Dir invariant). Use Plane::new() for normalized construction.
    pub normal: Vec3,
    /// U-axis direction (X direction of OCCT gp_Ax3). Orthogonal to normal.
    /// Determines the U=0 direction in the plane's (u, v) parameterization.
    pub u_dir: Vec3,
    /// V-axis direction (Y direction of OCCT gp_Ax3, = normal × u_dir).
    /// Orthogonal to both normal and u_dir.
    pub v_dir: Vec3,
}

impl Plane {
    /// OCCT-aligned: construct from origin and normal.
    /// Equivalent to `gp_Pln(gp_Pnt, gp_Dir)` which internally creates
    /// a `gp_Ax3(P, V)` whose X direction is the unit vector perpendicular
    /// to V having a zero in the coordinate of the smallest |component| of V
    /// (gp_Ax3.cxx L29-80):
    ///   1. A,B,C = V.X,V.Y,V.Z; Aabs,Babs,Cabs = |A|,|B|,|C|
    ///   2. If |B| is smallest:  D = |A|>|C| ? (-C,0,A) : (C,0,-A)
    ///   3. elif |A| is smallest: D = |B|>|C| ? (0,-C,B) : (0,C,-B)
    ///   4. else:                 D = |A|>|B| ? (-B,A,0) : (B,-A,0)
    ///   5. u_dir = D.normalize(); v_dir = N × u_dir
    pub fn new(origin: DVec3, normal: DVec3) -> Self {
        let normal = normal.normalize_or_zero();
        let (a, b, c) = (normal.x, normal.y, normal.z);
        let (aabs, babs, cabs) = (a.abs(), b.abs(), c.abs());
        let d = if babs <= aabs && babs <= cabs {
            if aabs > cabs {
                DVec3::new(-c, 0.0, a)
            } else {
                DVec3::new(c, 0.0, -a)
            }
        } else if aabs <= babs && aabs <= cabs {
            if babs > cabs {
                DVec3::new(0.0, -c, b)
            } else {
                DVec3::new(0.0, c, -b)
            }
        } else if aabs > babs {
            DVec3::new(-b, a, 0.0)
        } else {
            DVec3::new(b, -a, 0.0)
        };
        let u_dir = d.normalize_or_zero();
        let v_dir = normal.cross(u_dir).normalize_or_zero();
        Plane {
            origin,
            normal,
            u_dir,
            v_dir,
        }
    }

    /// OCCT-aligned: construct from origin, normal, and explicit u_dir.
    /// Equivalent to `gp_Pln(gp_Ax3(origin, normal, u_dir))`.
    /// v_dir is computed as `normal × u_dir` (right-handed orthonormal frame).
    pub fn with_axes(origin: DVec3, normal: DVec3, u_dir: DVec3) -> Self {
        let normal = normal.normalize_or_zero();
        let u_dir = u_dir.normalize();
        let v_dir = normal.cross(u_dir).normalize_or_zero();
        Plane {
            origin,
            normal,
            u_dir,
            v_dir,
        }
    }

    /// OCCT gp_Pln::Transform(theT) — gp_Pln.hxx L274:
    /// `myPosition.Transform(theT)` = gp_Ax3::Transform (gp_Ax3.hxx L306-311)
    /// — the Location point, the X direction and the Y direction are all
    /// transformed (the main direction follows as the axis transform).
    pub fn transform(&mut self, trsf: &crate::math::gp::Trsf) {
        self.origin = trsf.apply(self.origin); // gp_Ax3: axis location.
        self.normal = trsf.transform_dir(self.normal); // gp_Ax3: axis direction.
        self.u_dir = trsf.transform_dir(self.u_dir); // gp_Ax3: vxdir.
        self.v_dir = trsf.transform_dir(self.v_dir); // gp_Ax3: vydir.
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CylindricalSurface {
    pub origin: Point3,
    pub axis: Vec3,
    pub radius: f64,
    /// Reference direction for u=0 (perpendicular to axis).
    /// Preserved through rotation so UV mapping stays consistent.
    pub ref_dir: Vec3,
    /// Optional explicit Y direction of the surface frame
    /// (P(u,v) = O + R*(cos u*X + sin u*Y) + v*axis, u = atan2(P·Y, P·X)).
    /// OCCT gp_Ax3 keeps the X/Y directions explicitly, and a swept lateral
    /// face (BRepSweep_Translation::MakeEmptyFace + GeomAdaptor::
    /// SurfaceOfLinearExtrusion::Cylinder with gp_Ax3::ZReverse — which
    /// reverses ONLY the axis, gp_Ax3.hxx L131) carries the generating
    /// circle's Y even when the axis is reversed — a left-handed frame whose
    /// u equals the circle's parameter keeps the sweep pcurves consistent.
    /// When None the right-handed Y = axis × ref_dir is used (the default for
    /// all other constructions).
    #[serde(default)]
    pub y_dir: Option<Vec3>,
}

impl CylindricalSurface {
    /// Create a cylinder with [`any_perpendicular(axis)`](any_perpendicular) as the reference direction.
    pub fn new(origin: Point3, axis: Vec3, radius: f64) -> Self {
        Self {
            origin,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            ref_dir: any_perpendicular(axis),
            y_dir: None,
        }
    }

    /// Create a cylinder with an explicit reference direction for u=0.
    pub fn new_with_ref_dir(origin: Point3, axis: Vec3, radius: f64, ref_dir: Vec3) -> Self {
        Self {
            origin,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            ref_dir: ref_dir.normalize_or_zero(),
            y_dir: None,
        }
    }

    /// Effective Y direction of the surface frame: the explicit `y_dir` when
    /// set (the left-handed swept-lateral frame), else axis × ref_dir.
    pub fn y_axis(&self) -> Vec3 {
        self.y_dir
            .unwrap_or_else(|| self.axis.cross(self.ref_dir))
            .normalize_or_zero()
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SphericalSurface {
    pub center: Point3,
    pub axis: Vec3,
    pub radius: f64,
    /// Reference direction for u=0 (perpendicular to axis).
    /// Preserved through rotation so UV mapping stays consistent.
    pub ref_dir: Vec3,
}

impl SphericalSurface {
    /// Direction perpendicular to ref_dir in the equatorial plane (axis × ref_dir).
    /// Used for UV mapping: U = atan2(dot(ref_dir_perp), dot(ref_dir)).
    pub fn ref_dir_perp(&self) -> Vec3 {
        self.axis.cross(self.ref_dir).normalize_or_zero()
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConicalSurface {
    /// Point on the cone axis where the surface radius equals `radius`.
    ///
    /// Historically this field was used as an apex for zero-radius primitive
    /// cones. For general conical surfaces, the true apex is derived from this
    /// reference point, `radius`, and `half_angle_rad`.
    pub apex: Point3,
    pub axis: Vec3,
    /// Radius of the reference circle at `apex`.
    pub radius: f64,
    pub half_angle_rad: f64,
    /// Reference direction for u=0 (perpendicular to axis).
    ///
    /// OCCT `gp_Cone` carries a full `gp_Ax3` (Position + XDirection); the
    /// XDirection defines the u=0 generatrix.  Preserved through rotation so UV
    /// mapping stays consistent — same convention as [`CylindricalSurface`] and
    /// [`SphericalSurface`].
    #[serde(default = "default_cone_ref_dir")]
    pub ref_dir: Vec3,
}

/// Serde fallback for the cone `ref_dir` field when reading legacy data that
/// predates the field.  New data always carries the explicit reference.
fn default_cone_ref_dir() -> DVec3 {
    DVec3::X
}

impl ConicalSurface {
    /// Create a cone with the default reference direction for u=0
    /// (`any_perpendicular(axis)`, matching OCCT `gp_Ax2` default).
    pub fn new(apex: Point3, axis: Vec3, radius: f64, half_angle_rad: f64) -> Self {
        Self {
            apex,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            half_angle_rad,
            ref_dir: any_perpendicular(axis),
        }
    }

    /// Create a cone with an explicit reference direction for u=0.
    pub fn new_with_ref_dir(apex: Point3, axis: Vec3, radius: f64, half_angle_rad: f64, ref_dir: Vec3) -> Self {
        Self {
            apex,
            axis: axis.normalize_or_zero(),
            radius: radius.abs(),
            half_angle_rad,
            ref_dir: ref_dir.normalize_or_zero(),
        }
    }

    pub fn axis_dir(&self) -> DVec3 {
        self.axis.normalize_or_zero()
    }

    pub fn apex_point(&self) -> DVec3 {
        let tan_half = self.half_angle_rad.tan();
        if tan_half.abs() < 1e-12 {
            self.apex
        } else {
            self.apex - self.axis_dir() * (self.radius / tan_half)
        }
    }

    pub fn axial_from_slant(&self, slant: f64) -> f64 {
        slant * self.half_angle_rad.cos()
    }

    pub fn slant_from_axial(&self, axial: f64) -> f64 {
        let cos_half = self.half_angle_rad.cos();
        if cos_half.abs() < 1e-12 {
            0.0
        } else {
            axial / cos_half
        }
    }

    pub fn radius_at_slant(&self, slant: f64) -> f64 {
        self.radius + slant * self.half_angle_rad.sin()
    }

    pub fn radius_at_axial(&self, axial: f64) -> f64 {
        self.radius + axial * self.half_angle_rad.tan()
    }

    /// UV coordinates of world point `p` relative to this conical surface.
    ///
    /// `u` = azimuth (−π, π], `v` = slant distance from the reference circle at
    /// `self.apex`, matching [`SurfaceEval::point_at`].  When `p` is off the
    /// surface the returned `(u, v)` corresponds to the closest point on the cone.
    pub fn world_to_uv(self, p: DVec3) -> DVec2 {
        let axis = self.axis_dir();
        let x_ax = self.ref_dir.normalize_or_zero();
        let y_ax = axis.cross(x_ax).normalize();
        let local = p - self.apex;
        let along = local.dot(axis);
        let perp = local - axis * along;
        let radial = perp.length();

        let u = if radial < 1e-15 {
            0.0
        } else {
            // OCCT ElSLib::ConeParameters: atan2(P·YDirection, P·XDirection) on
            // the raw offset vector — NOT on the normalized radial. Normalizing
            // re-introduces rounding that flips the sign of a ~0 Y component
            // (atan2(-1e-16, X) -> 2PI), projecting a u=0 generatrix vertex to
            // u=2PI (same rationale as closest_point_on_surface Cone branch).
            local.dot(y_ax).atan2(local.dot(x_ax))
        };
        // OCCT gp_Cone / ElSLib::ConeD0 parameterization: u in [0, 2*PI].
        // atan2 returns (-PI, PI]; normalize so section-edge pcurves and the
        // FClass2d boundary sampling share the [0, 2*PI] domain with the
        // cone's natural pcurves (make_cone lateral/edge pcurves).
        let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };

        let cos_half = self.half_angle_rad.cos();
        let sin_half = self.half_angle_rad.sin();
        let v = along * cos_half + (radial - self.radius) * sin_half;

        DVec2::new(u, v)
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ToroidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    #[serde(default = "torus_ref_dir_default")]
    pub ref_dir: Vec3,
    pub major_radius: f64,
    pub minor_radius: f64,
}

fn torus_ref_dir_default() -> Vec3 {
    DVec3::X
}

/// An ellipsoidal surface aligned to a local orthonormal frame.
///
/// Parameterization matches sphere-like angles:
/// - `u` = longitude `[0, 2π]`
/// - `v` = colatitude `[0, π]` (0 at +axis pole)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EllipsoidalSurface {
    pub center: Point3,
    pub axis: Vec3,
    /// Reference direction used to derive the local X axis.
    pub ref_dir: Vec3,
    pub radius_x: f64,
    pub radius_y: f64,
    pub radius_z: f64,
}

/// A classical helicoid surface around an axis.
///
/// Parameterization:
/// `S(u, v) = origin + v * (cos(u) * x_axis + sin(u) * y_axis) + (pitch/(2*pi))*u * axis`
///
/// `u` is the azimuth / screw parameter and `v` is the signed radial distance.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct HelicoidSurface {
    pub origin: Point3,
    pub axis: Vec3,
    /// Reference direction used to derive the local X axis.
    pub ref_dir: Vec3,
    /// Axial advance per full revolution.
    pub pitch: f64,
}

/// A circular pipe/tube surface around a spine curve.
///
/// `u` is the azimuth angle around the local section frame and `v` follows the
/// natural parameter of the spine curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipeSurface {
    pub spine: Box<Curve3>,
    /// Initial/reference direction projected onto the normal plane of the
    /// spine tangent at evaluation time.
    pub ref_dir: Vec3,
    pub radius: f64,
}

/// A non-uniform rational B-spline surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BSplineSurface {
    pub degree_u: usize,
    pub degree_v: usize,
    /// Full knot vector for u (with multiplicities expanded).
    pub knots_u: Vec<f64>,
    /// Full knot vector for v (with multiplicities expanded).
    pub knots_v: Vec<f64>,
    /// Control point grid [u_index][v_index].
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight grid [u_index][v_index]; 1.0 for non-rational.
    pub weights: Vec<Vec<f64>>,
    /// OCCT Geom_BSplineSurface::IsUPeriodic() — the caller-declared
    /// constructor argument myUPeriodic (Geom_BSplineSurface.cxx ctor
    /// L230-249, CheckSurfaceData).
    #[serde(default)]
    pub is_periodic_u: bool,
    /// OCCT Geom_BSplineSurface::IsVPeriodic() — myVPeriodic.
    #[serde(default)]
    pub is_periodic_v: bool,
}

impl BSplineSurface {
    /// OCCT static Rational (Geom_BSplineSurface.cxx L110-138) — the V
    /// direction flag: set when two vertically adjacent weights differ by
    /// more than Epsilon(x) (the nextafter ULP, Standard_Real.hxx L242-246).
    /// OCCT derives myVRational at construction from the weight grid; the
    /// rcad surface is an immutable value, so the same derivation evaluates
    /// on read.
    pub fn is_rational_v(&self) -> bool {
        let w = &self.weights;
        for j in 0..w.first().map(|r| r.len()).unwrap_or(0) {
            for i in 0..w.len().saturating_sub(1) {
                let a = w[i][j];
                let b = w.get(i + 1).and_then(|r| r.get(j)).copied().unwrap_or(a);
                if (a - b).abs() > standard_epsilon(a) {
                    return true;
                }
            }
        }
        false
    }

    /// OCCT static Rational (Geom_BSplineSurface.cxx L126-137) — the U
    /// direction counterpart over horizontally adjacent weights.
    pub fn is_rational_u(&self) -> bool {
        let w = &self.weights;
        for (i, row) in w.iter().enumerate() {
            for j in 0..row.len().saturating_sub(1) {
                let a = row[j];
                let b = row.get(j + 1).copied().unwrap_or(a);
                if (a - b).abs() > standard_epsilon(a) {
                    return true;
                }
            }
        }
        false
    }
}

// OCCT Epsilon(theValue) (Standard_Real.hxx L242-248) — the ULP of `x`
// via nextafter toward the infinity of the same sign.  `epsilon_of` in
// `base::extrema_ext_elc` is the single canonical definition; this
// re-export keeps the historical local name for the call sites.
use crate::base::extrema_ext_elc::epsilon_of as standard_epsilon;

/// Returns `true` if the BSpline surface is planar (degree ≤ 1 in both directions
/// and all control points lie within `tol` of a single plane).
pub fn bspline_is_planar(bsp: &BSplineSurface, tol: f64) -> bool {
    if bsp.control_points.is_empty() {
        return false;
    }

    // Collect all unique control points
    let pts: Vec<DVec3> = bsp
        .control_points
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .collect();
    if pts.len() < 3 {
        return true; // trivially planar
    }

    // Find the first 3 non-collinear points to define the plane
    let origin = pts[0];
    let mut normal = DVec3::ZERO;
    for i in 1..pts.len() - 1 {
        let d1 = pts[i] - origin;
        let d2 = pts[i + 1] - origin;
        let n = d1.cross(d2);
        if n.length_squared() > tol * tol {
            normal = n.normalize();
            break;
        }
    }
    if normal.length_squared() < 0.5 {
        // All points are collinear — any plane containing the line works
        return true;
    }

    // Check all points lie within tol of the plane
    pts.iter().all(|&p| {
        let d = (p - origin).dot(normal);
        d.abs() <= tol
    })
}

/// Convert a planar BSpline surface to the best-fit `Plane`.
/// The BSpline must satisfy `bspline_is_planar`. Returns the plane from the
/// first 3 non-collinear control points (or any plane if all points collinear).
pub fn bspline_to_plane(bsp: &BSplineSurface) -> Plane {
    let pts: Vec<DVec3> = bsp
        .control_points
        .iter()
        .flat_map(|row| row.iter())
        .copied()
        .collect();

    let origin = pts[0];
    let mut normal = DVec3::Z;
    for i in 1..pts.len() - 1 {
        let d1 = pts[i] - origin;
        let d2 = pts[i + 1] - origin;
        let n = d1.cross(d2);
        if n.length_squared() > 1e-30 {
            normal = n.normalize();
            break;
        }
    }

    Plane::new(origin, normal)
}

/// A rational or non-rational Bezier surface (tensor-product bicubic patch).
///
/// Evaluated by applying de Casteljau in u, then in v. Domain is `[0, 1] × [0, 1]`.
/// Analogous to OCCT `Geom_BezierSurface`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BezierSurface {
    /// Control point grid [u_count][v_count].
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight grid [u_count][v_count]; 1.0 for non-rational.
    pub weights: Vec<Vec<f64>>,
}

/// A triangular rational Bezier surface using barycentric coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TriBezierSurface {
    /// Triangular control net rows. Row `i` has `degree + 1 - i` points.
    pub control_points: Vec<Vec<DVec3>>,
    /// Weight rows with the same triangular layout as `control_points`.
    pub weights: Vec<Vec<f64>>,
}

/// A surface offset from a base surface by a fixed distance along the normal.
///
/// `S(u,v) = basis.point_at(u,v) + offset_distance * basis.normal_at(u,v)`
///
/// The offset normal is the same as the basis normal. Domain equals the basis domain.
/// Analogous to OCCT `Geom_OffsetSurface`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OffsetSurface {
    pub basis: Box<Surface3>,
    /// Offset distance along the outward normal (positive = outward).
    pub offset_distance: f64,
}

/// A rectangular trimmed surface: a base surface restricted to the UV box
/// `[u1, u2] × [v1, v2]`.
///
/// Evaluation delegates fully to the basis surface; only the reported domain
/// changes. Analogous to OCCT `Geom_RectangularTrimmedSurface`.
///
/// Appears in STEP as `RECTANGULAR_TRIMMED_SURFACE`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmedSurface {
    /// The underlying surface being trimmed.
    pub basis: Box<Surface3>,
    /// Trim bounds `[u1, u2, v1, v2]`.
    pub trim: [f64; 4],
}

impl TrimmedSurface {
    pub fn new(basis: Surface3, u1: f64, u2: f64, v1: f64, v2: f64) -> Self {
        Self {
            basis: Box::new(basis),
            trim: [u1, u2, v1, v2],
        }
    }
}

/// Surface formed by translating a 3D profile curve along a direction.
/// S(u,v) = profile.point_at(u) + v * direction
/// Analogous to OCCT Geom_SurfaceOfLinearExtrusion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LinearExtrusionSurface {
    pub profile: Box<Curve3>,
    /// Normalized extrusion direction.
    pub direction: Vec3,
}

/// Surface formed by rotating a 3D profile curve around an axis.
/// S(u,v) = rotate(profile.point_at(v), axis_origin, axis_dir, angle=u)
/// Analogous to OCCT Geom_SurfaceOfRevolution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RevolutionSurface {
    pub profile: Box<Curve3>,
    pub axis_origin: Point3,
    /// Normalized rotation axis direction.
    pub axis_dir: Vec3,
}

/// Surface linearly interpolating between two 3D curves with a shared parameter domain.
/// S(u,v) = lerp(start.point_at(u), end.point_at(u), v)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RuledSurface {
    pub start: Box<Curve3>,
    pub end: Box<Curve3>,
}

/// A Coons patch blending four boundary curves over `[0,1] x [0,1]`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoonsSurface {
    /// Boundary curve at `v = 0`, parameterized by `u`.
    pub south: Box<Curve3>,
    /// Boundary curve at `v = 1`, parameterized by `u`.
    pub north: Box<Curve3>,
    /// Boundary curve at `u = 0`, parameterized by `v`.
    pub west: Box<Curve3>,
    /// Boundary curve at `u = 1`, parameterized by `v`.
    pub east: Box<Curve3>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Surface3 {
    Plane(Plane),
    Cylinder(CylindricalSurface),
    Sphere(SphericalSurface),
    Cone(ConicalSurface),
    Torus(ToroidalSurface),
    Ellipsoid(EllipsoidalSurface),
    Helicoid(HelicoidSurface),
    Pipe(PipeSurface),
    BSpline(BSplineSurface),
    LinearExtrusion(LinearExtrusionSurface),
    Revolution(RevolutionSurface),
    Ruled(RuledSurface),
    Coons(CoonsSurface),
    Bezier(BezierSurface),
    TriBezier(TriBezierSurface),
    Offset(OffsetSurface),
    Trimmed(TrimmedSurface),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum PrimitiveSolid {
    Box {
        width: f64,
        height: f64,
        depth: f64,
    },
    Sphere {
        radius: f64,
    },
    Cylinder {
        radius: f64,
        height: f64,
    },
    Cone {
        base_radius: f64,
        height: f64,
    },
    Torus {
        major_radius: f64,
        minor_radius: f64,
    },
}
