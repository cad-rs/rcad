//! OCCT GProp_GProps (TKGeomBase/GProp) — the global properties accumulator.
//!
//! 1:1 translation of:
//! - `GProp_GProps.hxx` L96-243 / `GProp_GProps.cxx` L26-197 (the
//!   accumulator: mass, centre of mass, matrix of inertia, static moments,
//!   moment about an axis, radius of gyration, principal properties, and the
//!   Add-form combination with the Huygens transfer),
//! - `GProp.cxx` L33-54 (`GProp::HOperator`, the Huygens operator),
//! - `GProp_PrincipalProps.hxx/.cxx` (the principal-properties presentation),
//! plus the BRepGProp drivers that fill the accumulator from the shape:
//! - `BRepGProp::LinearProperties` (BRepGProp.cxx L121-165) with the
//!   `BRepGProp_Cinert` edge integrator (BRepGProp_Cinert.cxx L30-172) and
//!   `BRepGProp_EdgeTool` (BRepGProp_EdgeTool.cxx L29-62 IntegrationOrder,
//!   L64-81 NbIntervals/Intervals -> GeomAdaptor_Curve::NbIntervals
//!   GeomAbs_CN -> BSplCLib::Intervals),
//! - `BRepGProp::SurfaceProperties` (BRepGProp.cxx L267-279 -> the static
//!   surfaceProperties L172-266) with the `BRepGProp_Sinert` FIXED-ORDER
//!   Gauss engine (BRepGProp_Gauss.cxx L1145-1211 domain Compute,
//!   L1306-1393 natural Compute, computeSInertiaOfElementaryPart L391-424,
//!   convert L473-491),
//! - `BRepGProp::VolumeProperties` (BRepGProp.cxx L555-586) aggregating the
//!   existing per-face Vinert components of `base::gprop::volume`
//!   (shape_vinert: the BRepGProp_Gauss.cxx L1215-1302 accumulation +
//!   computeVInertiaOfElementaryPart isByPoint L306-339) with the
//!   GProp_GProps::Add L48-64 combination and the Vinert convert L493-532.
//!
//! The drivers are the POOL-FREE (shape-tree) forms: the offset-engine
//! consumers carry `Shape` handles that live outside a `BRep` pool (the
//! `curve_pool_free` / `curve_on_surface_pool_free` convention in
//! topo::topods).  Limitations inherited from that convention (documented at
//! each site): TopLoc_Location is treated as identity, the MeshCinert
//! (triangulation) branches are not taken, and shapes stored in a pool
//! (Solid/Shell children as ShapeRef indices) are not reachable from a bare
//! `Shape`.

use std::collections::HashMap;

use glam::DVec3;

use crate::base::bnd_lib::curve2d_bounding_box;
use crate::base::gprop::surface::{
    add_inf, curve_integration_order, mult_inf, occt_gauss, surface_u_v_integration_order, EPS_DIM,
};
use crate::core::precision::{is_infinite_value, PCONFUSION};
use crate::geom::{Curve2dEval, Curve3, CurveEval, Surface3, SurfaceEval};
use crate::math::math_jacobi::MathJacobi;
use crate::math::MatD;
use crate::topo::topods::{curve_on_surface_pool_free, Orientation, ShapeType, TShape};
use crate::topo::topo_shape::Shape;

// OCCT math::GaussPointsMax() (math.cxx L25-28).
const GPM: usize = 61;

// OCCT gp::Resolution() (gp.hxx L60) = Standard::RealSmall() = DBL_MIN.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

// OCCT GProp_GProps.cxx L54 / L74: the |dim| guard of the gravity center.
const DIM_EPS: f64 = 1.0e-20;

/// OCCT Standard::Epsilon(theValue) (Standard_Real.hxx L242-248): the
/// absolute value of the difference between theValue and the nearest double
/// in the direction of infinity of the same sign; the minimal positive
/// double for 0.
fn occt_epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        f64::from_bits(the_value.to_bits() + 1) - the_value
    } else {
        the_value - f64::from_bits(the_value.to_bits() + 1)
    }
}

// ---------------------------------------------------------------------------
// gp_Mat stand-in — the member operations the translated GProp code uses.
// ---------------------------------------------------------------------------

/// OCCT gp_Mat — the 3x3 matrix with the row-major storage
/// `myMat[row][col]` of gp_Mat.hxx L80 and the member operations used by the
/// GProp translation (SetCols, Multiply, the arithmetic forms, the
/// gp_XYZ::Multiplied matrix-vector product).
#[derive(Debug, Clone, Copy)]
pub struct GpMat {
    a: [[f64; 3]; 3],
}

impl GpMat {
    /// OCCT gp_Mat() — the zero matrix (gp_Mat.hxx).
    pub const fn zero() -> Self {
        GpMat { a: [[0.0; 3]; 3] }
    }

    /// OCCT gp_Mat(theCol1, theCol2, theCol3) (gp_Mat.cxx L30-42): the three
    /// gp_XYZ arguments are the COLUMNS.
    pub fn from_cols(the_col1: DVec3, the_col2: DVec3, the_col3: DVec3) -> Self {
        let mut m = GpMat::zero();
        m.a[0][0] = the_col1.x;
        m.a[1][0] = the_col1.y;
        m.a[2][0] = the_col1.z;
        m.a[0][1] = the_col2.x;
        m.a[1][1] = the_col2.y;
        m.a[2][1] = the_col2.z;
        m.a[0][2] = the_col3.x;
        m.a[1][2] = the_col3.y;
        m.a[2][2] = the_col3.z;
        m
    }

    /// OCCT gp_Mat::SetCols(Xyz1, Xyz2, Xyz3).
    pub fn set_cols(&mut self, the_xyz1: DVec3, the_xyz2: DVec3, the_xyz3: DVec3) {
        *self = GpMat::from_cols(the_xyz1, the_xyz2, the_xyz3);
    }

    /// OCCT gp_Mat::Multiply(Scalar) — in place.
    pub fn multiply(&mut self, the_scalar: f64) {
        for r in 0..3 {
            for c in 0..3 {
                self.a[r][c] *= the_scalar;
            }
        }
    }

    /// OCCT gp_Mat::Multiplied(Scalar).
    pub fn multiplied(&self, the_scalar: f64) -> GpMat {
        let mut m = *self;
        m.multiply(the_scalar);
        m
    }

    /// OCCT operator+(gp_Mat, gp_Mat).
    pub fn added(&self, the_other: GpMat) -> GpMat {
        let mut m = GpMat::zero();
        for r in 0..3 {
            for c in 0..3 {
                m.a[r][c] = self.a[r][c] + the_other.a[r][c];
            }
        }
        m
    }

    /// OCCT operator-(gp_Mat, gp_Mat).
    pub fn subbed(&self, the_other: GpMat) -> GpMat {
        let mut m = GpMat::zero();
        for r in 0..3 {
            for c in 0..3 {
                m.a[r][c] = self.a[r][c] - the_other.a[r][c];
            }
        }
        m
    }

    /// OCCT gp_Mat::Value(i, j) — 1-based element read.
    pub fn value(&self, the_i: usize, the_j: usize) -> f64 {
        self.a[the_i - 1][the_j - 1]
    }

    /// OCCT gp_XYZ::Multiplied(theMatrix) (gp_XYZ.hxx L339-347):
    /// result(row i) = sum_j myMat[i][j] * v[j] — the M * v product.
    pub fn mul_vec(&self, the_v: DVec3) -> DVec3 {
        DVec3::new(
            self.a[0][0] * the_v.x + self.a[0][1] * the_v.y + self.a[0][2] * the_v.z,
            self.a[1][0] * the_v.x + self.a[1][1] * the_v.y + self.a[1][2] * the_v.z,
            self.a[2][0] * the_v.x + self.a[2][1] * the_v.y + self.a[2][2] * the_v.z,
        )
    }
}

// ---------------------------------------------------------------------------
// GProp::HOperator (GProp.cxx L33-54).
// ---------------------------------------------------------------------------

/// OCCT GProp::HOperator(G, Q, Mass, Operator) — the "Huyghens operator" of
/// a system of mass `the_mass` with center of mass `the_g` at the point
/// `the_q`: Inertia/Q = Inertia/G + HOperator.
pub fn h_operator(the_g: DVec3, the_q: DVec3, the_mass: f64) -> GpMat {
    // OCCT L39-52.
    let qg = the_g - the_q;
    let ixx = qg.y * qg.y + qg.z * qg.z;
    let iyy = qg.x * qg.x + qg.z * qg.z;
    let izz = qg.y * qg.y + qg.x * qg.x;
    let ixy = -qg.x * qg.y;
    let iyz = -qg.y * qg.z;
    let ixz = -qg.x * qg.z;
    let mut operator = GpMat::zero();
    operator.set_cols(
        DVec3::new(ixx, ixy, ixz),
        DVec3::new(ixy, iyy, iyz),
        DVec3::new(ixz, iyz, izz),
    );
    operator.multiply(the_mass);
    operator
}

// ---------------------------------------------------------------------------
// GProp_PrincipalProps (hxx; cxx L27-140).
// ---------------------------------------------------------------------------

/// OCCT GProp_PrincipalProps — the principal properties of inertia of a
/// system computed by a GProp_GProps object.
#[derive(Debug, Clone, Copy)]
pub struct PrincipalProps {
    /// OCCT: i1 — the first principal moment.
    pub i1: f64,
    /// OCCT: i2 — the second principal moment.
    pub i2: f64,
    /// OCCT: i3 — the third principal moment.
    pub i3: f64,
    /// OCCT: r1 — the first principal radius of gyration.
    pub r1: f64,
    /// OCCT: r2 — the second principal radius of gyration.
    pub r2: f64,
    /// OCCT: r3 — the third principal radius of gyration.
    pub r3: f64,
    /// OCCT: v1 — the first axis of inertia.
    pub v1: DVec3,
    /// OCCT: v2 — the second axis of inertia.
    pub v2: DVec3,
    /// OCCT: v3 — the third axis of inertia.
    pub v3: DVec3,
    /// OCCT: g — the center of mass of the system.
    pub g: DVec3,
}

impl PrincipalProps {
    /// OCCT GProp_PrincipalProps::HasSymmetryAxis() (cxx L57-65): the
    /// relative tolerance 1.e-10 form.
    pub fn has_symmetry_axis(&self) -> bool {
        let a_rel_tol = 1.0e-10;
        let eps1 = self.i1.abs() * a_rel_tol;
        let eps2 = self.i2.abs() * a_rel_tol;
        (self.i1 - self.i2).abs() <= eps1
            || (self.i1 - self.i3).abs() <= eps1
            || (self.i2 - self.i3).abs() <= eps2
    }

    /// OCCT GProp_PrincipalProps::HasSymmetryAxis(aTol) (cxx L67-76).
    pub fn has_symmetry_axis_tol(&self, a_tol: f64) -> bool {
        let eps1 = (self.i1 * a_tol).abs() + occt_epsilon(self.i1).abs();
        let eps2 = (self.i2 * a_tol).abs() + occt_epsilon(self.i2).abs();
        (self.i1 - self.i2).abs() <= eps1
            || (self.i1 - self.i3).abs() <= eps1
            || (self.i2 - self.i3).abs() <= eps2
    }

    /// OCCT GProp_PrincipalProps::HasSymmetryPoint() (cxx L78-85).
    pub fn has_symmetry_point(&self) -> bool {
        let a_rel_tol = 1.0e-10;
        let eps1 = self.i1.abs() * a_rel_tol;
        (self.i1 - self.i2).abs() <= eps1 && (self.i1 - self.i3).abs() <= eps1
    }

    /// OCCT GProp_PrincipalProps::HasSymmetryPoint(aTol) (cxx L87-94).
    pub fn has_symmetry_point_tol(&self, a_tol: f64) -> bool {
        let eps1 = (self.i1 * a_tol).abs() + occt_epsilon(self.i1).abs();
        (self.i1 - self.i2).abs() <= eps1 && (self.i1 - self.i3).abs() <= eps1
    }

    /// OCCT GProp_PrincipalProps::Moments(Ixx, Iyy, Izz) (cxx L96-102).
    pub fn moments(&self) -> (f64, f64, f64) {
        (self.i1, self.i2, self.i3)
    }

    /// OCCT GProp_PrincipalProps::FirstAxisOfInertia() (cxx L104-107).
    pub fn first_axis_of_inertia(&self) -> DVec3 {
        self.v1
    }

    /// OCCT GProp_PrincipalProps::SecondAxisOfInertia() (cxx L109-112).
    pub fn second_axis_of_inertia(&self) -> DVec3 {
        self.v2
    }

    /// OCCT GProp_PrincipalProps::ThirdAxisOfInertia() (cxx L114-117).
    pub fn third_axis_of_inertia(&self) -> DVec3 {
        self.v3
    }

    /// OCCT GProp_PrincipalProps::RadiusOfGyration(Rxx, Ryy, Rzz)
    /// (cxx L119-125).
    pub fn radius_of_gyration(&self) -> (f64, f64, f64) {
        (self.r1, self.r2, self.r3)
    }
}

// ---------------------------------------------------------------------------
// GProp_GProps (hxx L96-243; cxx L26-197).
// ---------------------------------------------------------------------------

/// OCCT GProp_GProps — the accumulator of the global properties of a
/// compound geometric system.  `g` is the center of mass OFFSET from `loc`
/// (the absolute-frame center of mass is `loc + g`); `inertia` is the
/// quadratic-moments matrix about `loc` in the gp_Mat layout
/// [[Ixx, -Ixy, -Ixz], [-Ixy, Iyy, -Iyz], [-Ixz, -Iyz, Izz]].
#[derive(Debug, Clone, Copy)]
pub struct GProps {
    g: DVec3,       // OCCT: g — center of mass (relative to loc)
    loc: DVec3,     // OCCT: loc — reference point for inertia accumulation
    dim: f64,       // OCCT: dim — mass / length / area / volume
    inertia: GpMat, // OCCT: inertia — quadratic moments matrix
}

impl GProps {
    /// OCCT GProp_GProps::GProp_GProps() (cxx L26-32): the origin (0, 0, 0)
    /// of the absolute Cartesian coordinate system is used.
    pub fn new() -> Self {
        GProps {
            g: DVec3::ZERO,  // gp::Origin()
            loc: DVec3::ZERO,
            dim: 0.0,
            inertia: GpMat::zero(),
        }
    }

    /// OCCT GProp_GProps::GProp_GProps(SystemLocation) (cxx L34-40).
    pub fn with_location(system_location: DVec3) -> Self {
        GProps {
            g: DVec3::ZERO,
            loc: system_location,
            dim: 0.0,
            inertia: GpMat::zero(),
        }
    }

    /// OCCT GProp_GProps::Add(Item, Density = 1.0) (cxx L42-98): brings
    /// together the properties of `the_item` with those of this framework;
    /// Huygens' theorem transfers the inertia to this reference point when
    /// the two locations differ.
    pub fn add(&mut self, the_item: &GProps, the_density: f64) {
        // OCCT L44-47.
        if the_density <= GP_RESOLUTION {
            panic!("Standard_DomainError: GProp_GProps::Add - Density <= gp::Resolution()");
        }
        // OCCT L48.
        if self.loc.distance(the_item.loc) <= GP_RESOLUTION {
            // OCCT L50-63.
            let mut gxyz = the_item.g * (the_item.dim * the_density);
            let g_scaled = self.g * self.dim;
            gxyz += g_scaled;
            self.dim += the_item.dim * the_density;
            if self.dim.abs() >= DIM_EPS {
                gxyz /= self.dim;
                self.g = gxyz;
            } else {
                self.g = DVec3::ZERO;
            }
            // OCCT L63.
            self.inertia = self.inertia.added(the_item.inertia.multiplied(the_density));
        } else {
            // OCCT L67-72.
            let itemloc = self.loc - the_item.loc;
            let itemg = the_item.loc + the_item.g;
            let mut gxyz = (the_item.g - itemloc) * (the_item.dim * the_density);
            let g_scaled = self.g * self.dim;
            gxyz += g_scaled;
            // OCCT L73-82.
            self.dim += the_item.dim * the_density;
            if self.dim.abs() >= DIM_EPS {
                gxyz /= self.dim;
                self.g = gxyz;
            } else {
                self.g = DVec3::ZERO;
            }
            // OCCT L83-96: the inertia of the Item at the location point of
            // the system via the Huyghens theorem.
            let mut item_inertia = the_item.inertia;
            if the_item.g.length() > GP_RESOLUTION {
                // Computes the inertia of Item at its dim centre.
                let h_mat = h_operator(itemg, the_item.loc, the_item.dim);
                item_inertia = item_inertia.subbed(h_mat);
            }
            // Computes the inertia of Item at the location point of the system.
            let h_mat = h_operator(itemg, self.loc, the_item.dim);
            item_inertia = item_inertia.added(h_mat);
            self.inertia = self.inertia.added(item_inertia.multiplied(the_density));
        }
    }

    /// OCCT GProp_GProps::Mass() (cxx L100-103).
    pub fn mass(&self) -> f64 {
        self.dim
    }

    /// OCCT GProp_GProps::CentreOfMass() (cxx L105-108).
    pub fn centre_of_mass(&self) -> DVec3 {
        self.loc + self.g
    }

    /// OCCT GProp_GProps::MatrixOfInertia() (cxx L110-115): the symmetric
    /// matrix of the quadratic moments about the center of mass.
    pub fn matrix_of_inertia(&self) -> GpMat {
        let h_mat = h_operator(self.g, DVec3::ZERO, self.dim);
        self.inertia.subbed(h_mat)
    }

    /// OCCT GProp_GProps::StaticMoments(Ix, Iy, Iz) (cxx L117-124): the
    /// moments of inertia about the three axes of the absolute coordinate
    /// system.
    pub fn static_moments(&self) -> (f64, f64, f64) {
        let g_abs = self.loc + self.g;
        (g_abs.x * self.dim, g_abs.y * self.dim, g_abs.z * self.dim)
    }

    /// OCCT GProp_GProps::MomentOfInertia(A) (cxx L126-146): the moment of
    /// inertia about the axis A.  `the_a_loc` is A.Location() and
    /// `the_a_dir` A.Direction() (a gp_Dir — unit by construction).
    pub fn moment_of_inertia(&self, the_a_loc: DVec3, the_a_dir: DVec3) -> f64 {
        // OCCT L134-137.
        if self.loc.distance(the_a_loc) <= GP_RESOLUTION {
            the_a_dir.dot(self.inertia.mul_vec(the_a_dir))
        } else {
            // OCCT L140-144.
            let h_mat = h_operator(self.loc + self.g, the_a_loc, self.dim);
            let axis_inertia = self.matrix_of_inertia().added(h_mat);
            the_a_dir.dot(axis_inertia.mul_vec(the_a_dir))
        }
    }

    /// OCCT GProp_GProps::RadiusOfGyration(A) (cxx L148-152).
    pub fn radius_of_gyration(&self, the_a_loc: DVec3, the_a_dir: DVec3) -> f64 {
        (self.moment_of_inertia(the_a_loc, the_a_dir) / self.dim).sqrt()
    }

    /// OCCT GProp_GProps::PrincipalProperties() (cxx L154-197): the
    /// eigenvalues and eigenvectors of the matrix of inertia (math_Jacobi).
    pub fn principal_properties(&self) -> PrincipalProps {
        // OCCT L157-165.
        let mut diag_mat = MatD::new(3, 3);
        let axis_inertia = self.matrix_of_inertia();
        for j in 1..=3 {
            for i in 1..=3 {
                diag_mat.set(i, j, axis_inertia.value(i, j));
            }
        }
        // OCCT L167-174.
        let jacobi = MathJacobi::new(&diag_mat);
        let ixx = jacobi.value(1);
        let iyy = jacobi.value(2);
        let izz = jacobi.value(3);
        let vectors = jacobi.vectors();
        let vxx = DVec3::new(vectors.get(1, 1), vectors.get(2, 1), vectors.get(3, 1));
        let vyy = DVec3::new(vectors.get(1, 2), vectors.get(2, 2), vectors.get(3, 2));
        let vzz = DVec3::new(vectors.get(1, 3), vectors.get(2, 3), vectors.get(3, 3));
        // OCCT L176-186: protection against dim == 0.
        let mut rxx = 0.0;
        let mut ryy = 0.0;
        let mut rzz = 0.0;
        if 0.0 != self.dim {
            rxx = (ixx / self.dim).abs().sqrt();
            ryy = (iyy / self.dim).abs().sqrt();
            rzz = (izz / self.dim).abs().sqrt();
        }
        // OCCT L187-196.
        PrincipalProps {
            i1: ixx,
            i2: iyy,
            i3: izz,
            r1: rxx,
            r2: ryy,
            r3: rzz,
            v1: vxx,
            v2: vyy,
            v3: vzz,
            g: self.g + self.loc,
        }
    }
}

// ---------------------------------------------------------------------------
// BRepGProp_EdgeTool (BRepGProp_EdgeTool.cxx L29-81) + the GeomAdaptor_Curve
// / BSplCLib interval machinery its NbIntervals(GeomAbs_CN) delegates to.
// ---------------------------------------------------------------------------

/// OCCT BRepAdaptor_Curve over an edge's stored 3D curve
/// (BRep_Tool::Curve(E, f, l) BRep_Tool.cxx L219-224 + GeomAdaptor_Curve::
/// Load of a Geom_TrimmedCurve: the BASIS curve with the intersected
/// parameter range).  Returns (basis curve, first, last).
fn edge_adaptor_curve(the_c: Curve3, the_first: f64, the_last: f64) -> (Curve3, f64, f64) {
    if let Curve3::Trimmed(tc) = &the_c {
        // GeomAdaptor_Curve::Load: myCurve = TC->BasisCurve();
        // myFirst = Max(TC->FirstParameter(), U1); myLast = Min(TC->LastParameter(), U2).
        let first = tc.first.max(the_first);
        let last = tc.last.min(the_last);
        ((*tc.curve).clone(), first, last)
    } else {
        (the_c, the_first, the_last)
    }
}

/// OCCT BRepGProp_EdgeTool::IntegrationOrder(BAC) (BRepGProp_EdgeTool.cxx
/// L29-62).  `the_curve` is the adaptor basis curve (the trimmed unwrap of
/// edge_adaptor_curve — BRepAdaptor_Curve::GetType exposes the basis type).
pub fn edge_tool_integration_order(the_curve: &Curve3) -> usize {
    match the_curve {
        // OCCT L33-36: GeomAbs_Line.
        Curve3::Line(_) => 2,
        // OCCT L37-40: GeomAbs_Parabola.
        Curve3::Parabola(_) => 5,
        // OCCT L41-48: GeomAbs_BezierCurve — 2 * NbPoles - 1.
        Curve3::Bezier(b) => 2 * b.control_points.len() - 1,
        // OCCT L49-56: GeomAbs_BSplineCurve — 2 * NbPoles - 1.
        Curve3::BSpline(b) => 2 * b.control_points.len() - 1,
        // OCCT L57-60: default (Circle, Ellipse, Hyperbola, Offset, ...).
        _ => 10,
    }
}

/// OCCT GeomAdaptor_Curve::NbIntervals(GeomAbs_CN) (GeomAdaptor_Curve.cxx
/// L371-410) for the C1-interval split of BRepGProp_Cinert::Perform: only a
/// BSpline curve splits (anOffsetCurve recursion into the basis with CN; a
/// Geom_OffsetCurve adaptor, not translated pool-free, behaves as 1); every
/// other curve type returns 1 (L458-461).
pub fn edge_tool_nb_intervals_cn(the_curve: &Curve3, the_first: f64, the_last: f64) -> usize {
    match the_curve {
        Curve3::BSpline(b) => {
            // OCCT L389-401: for GeomAbs_CN aCont = aDegree, so BSplCLib::
            // Intervals keeps every knot (mult > degree - degree = 0).
            let a_cont = b.degree;
            bsplib_intervals(
                &compressed_knots(&b.knots),
                b.degree,
                b.is_periodic,
                a_cont,
                the_first,
                the_last,
            )
            .0
        }
        Curve3::Offset(off) => {
            // OCCT L413-455: the offset adaptor recurses into the basis with
            // BaseS = GeomAbs_CN and counts the basis intervals inside the
            // offset range.
            let (basis, b_first, b_last) = edge_adaptor_curve((*(off.basis)).clone(), the_first, the_last);
            edge_tool_nb_intervals_cn(&basis, b_first, b_last)
        }
        // OCCT L458-461.
        _ => 1,
    }
}

/// OCCT BSplCLib::Intervals (BSplCLib.cxx L4824-4963): the interval array of
/// a BSpline for continuity `the_continuity` over [the_first, the_last].
/// Returns (NbIntervals, TI[0..=NbIntervals]) — TI is the OCCT
/// 1-based array stored 0-based.
pub fn bsplib_intervals(
    the_knots: &[f64],
    the_degree: usize,
    the_is_periodic: bool,
    the_continuity: usize,
    the_first: f64,
    the_last: f64,
) -> (usize, Vec<f64>) {
    let _ = the_degree; // LocateParameter(theDegree, ...) — unused in the no-multiplicity form
    let _ = the_continuity; // folded into the caller's knot filtering (CN keeps every knot)
    // OCCT L4834-4849: keep the knots with mult > (degree - continuity)
    // except first and last.  The rcad knot array arrives COMPRESSED (each
    // entry is a distinct knot; the multiplicities are folded away), and the
    // caller passes aCont = degree for CN so every knot qualifies.
    let mut a_new_knots: Vec<f64> = the_knots.to_vec();
    if a_new_knots.is_empty() {
        a_new_knots.push(the_first);
        a_new_knots.push(the_last);
    }
    let a_nb_new_knots = a_new_knots.len();

    // OCCT L4851-4883: the periodic range wrapping.
    let mut a_cur_first = the_first;
    let mut a_cur_last = the_last;
    let mut a_period = 0.0;
    let mut a_first_period: i64 = 0;
    let mut a_last_period: i64 = 0;
    if the_is_periodic {
        let a_lower = a_new_knots[0];
        let an_upper = a_new_knots[a_nb_new_knots - 1];
        a_period = an_upper - a_lower;
        while a_cur_first < a_lower {
            a_cur_first += a_period;
            a_first_period -= 1;
        }
        while a_cur_last < a_lower {
            a_cur_last += a_period;
            a_last_period -= 1;
        }
        while a_cur_first >= an_upper {
            a_cur_first -= a_period;
            a_first_period += 1;
        }
        while a_cur_last >= an_upper {
            a_cur_last -= a_period;
            a_last_period += 1;
        }
    }

    // OCCT L4885-4908: LocateParameter without multiplicities — the largest
    // knot index holding knots(index) <= U (1-based).
    let locate = |the_u: f64| -> usize {
        if the_u <= a_new_knots[0] {
            1
        } else if the_u >= a_new_knots[a_nb_new_knots - 1] {
            a_nb_new_knots - 1
        } else {
            let mut index = 1;
            for i in 1..=a_nb_new_knots {
                if a_new_knots[i - 1] <= the_u {
                    index = i;
                } else {
                    break;
                }
            }
            index
        }
    };
    let mut an_index1 = locate(a_cur_first);
    let mut an_index2 = locate(a_cur_last);

    // OCCT L4910-4917: the coincidence bumps at the range boundaries.
    let an_eps = PCONFUSION;
    if an_index1 < a_nb_new_knots && (a_new_knots[an_index1] - a_cur_first).abs() < an_eps {
        an_index1 += 1;
    }
    if (a_new_knots[an_index2 - 1] - a_cur_last).abs() < an_eps {
        an_index2 -= 1;
    }
    // OCCT L4918.
    let a_nb_intervals =
        (an_index2 as i64 - an_index1 as i64 + 1 + (a_last_period - a_first_period) * (a_nb_new_knots as i64 - 1)) as usize;

    // OCCT L4921-4959: the interval array.
    let mut the_intervals = vec![0.0; a_nb_intervals + 1];
    if the_is_periodic && a_last_period != a_first_period {
        let mut an_index = 1usize;
        // Part from the beginning of range to the end of the first period.
        for i in an_index1..a_nb_new_knots {
            the_intervals[an_index - 1] = a_new_knots[i - 1] + a_first_period as f64 * a_period;
            an_index += 1;
        }
        // Full periods.
        for a_period_num in (a_first_period + 1)..a_last_period {
            for i in 1..a_nb_new_knots {
                the_intervals[an_index - 1] = a_new_knots[i - 1] + a_period_num as f64 * a_period;
                an_index += 1;
            }
        }
        // Part from the beginning of the last period to the end of range.
        for i in 1..=an_index2 {
            the_intervals[an_index - 1] = a_new_knots[i - 1] + a_last_period as f64 * a_period;
            an_index += 1;
        }
    } else {
        let mut an_index = 1usize;
        for i in an_index1..=an_index2 {
            the_intervals[an_index - 1] = a_new_knots[i - 1] + a_first_period as f64 * a_period;
            an_index += 1;
        }
    }
    // OCCT L4952-4957: the range boundaries overwrite the first and last
    // slots.
    the_intervals[0] = the_first;
    the_intervals[a_nb_intervals] = the_last;
    (a_nb_intervals, the_intervals)
}

/// Geom_BSplineCurve::Knots — the COMPRESSED knot vector (the rcad knot
/// vectors are stored with multiplicities expanded).
fn compressed_knots(knots: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if out.is_empty() || (k - out[out.len() - 1]).abs() > 1e-15 {
            out.push(k);
        }
    }
    if out.is_empty() {
        out.push(0.0);
    }
    out
}

/// OCCT BRepGProp_EdgeTool::Intervals(C, T, GeomAbs_CN) — the interval
/// array TI(1, NbIntervals + 1) of BRepGProp_Cinert::Perform.
pub fn edge_tool_intervals_cn(
    the_curve: &Curve3,
    the_first: f64,
    the_last: f64,
) -> (usize, Vec<f64>) {
    match the_curve {
        Curve3::BSpline(b) => bsplib_intervals(
            &compressed_knots(&b.knots),
            b.degree,
            b.is_periodic,
            b.degree,
            the_first,
            the_last,
        ),
        Curve3::Offset(off) => {
            let (basis, b_first, b_last) = edge_adaptor_curve((*(off.basis)).clone(), the_first, the_last);
            edge_tool_intervals_cn(&basis, b_first, b_last)
        }
        _ => (1, vec![the_first, the_last]),
    }
}

// ---------------------------------------------------------------------------
// BRepGProp_Cinert (BRepGProp_Cinert.hxx; .cxx L30-172).
// ---------------------------------------------------------------------------

/// OCCT BRepGProp_Cinert — the global properties of a bounded curve
/// (Rust has no inheritance: the GProp_GProps base state is the `base`
/// field, the Stage 2e facade precedent).
#[derive(Debug, Clone, Copy)]
pub struct BRepGPropCinert {
    pub base: GProps, // OCCT: : public GProp_GProps
}

impl BRepGPropCinert {
    /// OCCT BRepGProp_Cinert::SetLocation(CLocation) (cxx L32-35).
    pub fn set_location(&mut self, the_clocation: DVec3) {
        self.base.loc = the_clocation;
    }

    /// OCCT BRepGProp_Cinert::Perform(C) (cxx L52-172) over the adaptor
    /// basis curve `the_curve` on [the_first, the_last].
    pub fn perform(&mut self, the_curve: &Curve3, the_first: f64, the_last: f64) {
        // OCCT L54-57.
        let mut i_x = 0.0;
        let mut i_y = 0.0;
        let mut i_z = 0.0;
        let mut i_xx = 0.0;
        let mut i_yy = 0.0;
        let mut i_zz = 0.0;
        let mut i_xy = 0.0;
        let mut i_xz = 0.0;
        let mut i_yz = 0.0;
        self.base.dim = 0.0;

        // OCCT L59-62: First/LastParameter + Order.
        let order = edge_tool_integration_order(the_curve).min(GPM);
        // OCCT L74-75.
        let (gauss_p, gauss_w) = match occt_gauss(order) {
            Some(v) => v,
            None => return,
        };

        // OCCT L84-101: the C1 (GeomAbs_CN) interval split.
        let nb_intervals = edge_tool_nb_intervals_cn(the_curve, the_first, the_last);
        let b_has_intervals = nb_intervals > 1;
        // OCCT L111-112.
        let uu1 = the_first.min(the_last);
        let uu2 = the_first.max(the_last);
        // OCCT L87-101: TI(1, nbIntervals + 1); the !bHasIntervals branch
        // keeps Lower = UU1 / Upper = UU2 (L106-109).
        let ti: Vec<f64> = if b_has_intervals {
            edge_tool_intervals_cn(the_curve, the_first, the_last).1
        } else {
            vec![uu1, uu2]
        };

        // OCCT L113-114: the nIndex interval loop.
        let mut p_last = DVec3::ZERO;
        for n_index in 1..=nb_intervals {
            // OCCT L116-125.
            let (lower, upper) = if b_has_intervals {
                (ti[n_index - 1].max(uu1), ti[n_index].min(uu2))
            } else {
                (uu1, uu2)
            };

            // OCCT L127-136.
            let mut dim_local = 0.0;
            let mut ix_local = 0.0;
            let mut iy_local = 0.0;
            let mut iz_local = 0.0;
            let mut ixx_local = 0.0;
            let mut iyy_local = 0.0;
            let mut izz_local = 0.0;
            let mut ixy_local = 0.0;
            let mut ixz_local = 0.0;
            let mut iyz_local = 0.0;

            // OCCT L138.
            let xloc = self.base.loc.x;
            let yloc = self.base.loc.y;
            let zloc = self.base.loc.z;

            // OCCT L146-148.
            let um = 0.5 * (upper + lower);
            let ur = 0.5 * (upper - lower);

            // OCCT L150-165: the Gauss loop.
            for i in 0..order {
                let u = um + ur * gauss_p[i];
                // OCCT L152-153: D1(C, u, P, V1).
                let p = the_curve.point_at(u);
                let v1 = the_curve.derivative_at(u);
                p_last = p;
                let mut ds = v1.length();
                let mut x = p.x;
                let mut y = p.y;
                let mut z = p.z;
                x -= xloc;
                y -= yloc;
                z -= zloc;
                // OCCT L157-165.
                ds *= gauss_w[i];
                dim_local += ds;
                ix_local += x * ds;
                iy_local += y * ds;
                iz_local += z * ds;
                ixy_local += x * y * ds;
                iyz_local += y * z * ds;
                ixz_local += x * z * ds;
                x *= x;
                y *= y;
                z *= z;
                ixx_local += (y + z) * ds;
                iyy_local += (x + z) * ds;
                izz_local += (x + y) * ds;
            }
            // OCCT L167-176: *= ur and accumulate.
            dim_local *= ur;
            ix_local *= ur;
            iy_local *= ur;
            iz_local *= ur;
            ixx_local *= ur;
            iyy_local *= ur;
            izz_local *= ur;
            ixy_local *= ur;
            ixz_local *= ur;
            iyz_local *= ur;

            self.base.dim += dim_local;
            i_x += ix_local;
            i_y += iy_local;
            i_z += iz_local;
            i_xx += ixx_local;
            i_yy += iyy_local;
            i_zz += izz_local;
            i_xy += ixy_local;
            i_xz += ixz_local;
            i_yz += iyz_local;
        }

        // OCCT L164-166.
        self.base.inertia = GpMat::from_cols(
            DVec3::new(i_xx, -i_xy, -i_xz),
            DVec3::new(-i_xy, i_yy, -i_yz),
            DVec3::new(-i_xz, -i_yz, i_zz),
        );

        // OCCT L168-173.
        if self.base.dim.abs() < GP_RESOLUTION {
            self.base.g = p_last;
        } else {
            self.base.g = DVec3::new(i_x / self.base.dim, i_y / self.base.dim, i_z / self.base.dim);
        }
    }
}

// ---------------------------------------------------------------------------
// Pool-free TopExp_Explorer carrier (the tree walk over Shape handles).
// ---------------------------------------------------------------------------

/// OCCT TopExp_Explorer(S, ToFind) — the pool-free walk over the Shape
/// tree with the cumulative orientation (TopAbs::Compose).  Architecture
/// limitation (the pool model): a Solid/Shell stores its children as
/// ShapeRef pool indices, unreachable from a bare Shape — those branches
/// yield nothing; the offset-engine shapes (wires, faces, compounds of
/// wires) are fully covered.
fn explorer_shapes(the_s: &Shape, the_to_find: ShapeType) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    explorer_walk(the_s, the_to_find, Orientation::Forward, &mut out);
    out
}

fn explorer_walk(
    the_s: &Shape,
    the_to_find: ShapeType,
    the_ori: Orientation,
    the_out: &mut Vec<Shape>,
) {
    if the_s.is_null() {
        return;
    }
    let a_type = the_s.shape_type();
    // TopAbs::Compose — the cumulated orientation of this occurrence.
    let cum = the_ori.compose(the_s.orientation);
    if a_type == the_to_find {
        // The oriented occurrence (TopExp_Explorer::Current()).
        let mut occ = the_s.clone();
        occ.orientation = cum;
        the_out.push(occ);
        // TopExp_Explorer does not descend below the requested type.
        return;
    }
    match a_type {
        ShapeType::Compound | ShapeType::CompSolid => {
            if let TShape::Compound(children) | TShape::CompSolid(children) = the_s.data.as_ref() {
                for c in children {
                    explorer_walk(c, the_to_find, cum, the_out);
                }
            }
        }
        ShapeType::Wire => {
            if let TShape::Wire(wd) = the_s.data.as_ref() {
                for e in &wd.edges {
                    explorer_walk(e, the_to_find, cum, the_out);
                }
            }
        }
        ShapeType::Face => {
            if let TShape::Face(fd) = the_s.data.as_ref() {
                // TopExp_Explorer(F, EDGE/VERTEX) traverses the wires.
                for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                    explorer_walk(w, the_to_find, cum, the_out);
                }
            }
        }
        ShapeType::Edge => {
            if let TShape::Edge(ed) = the_s.data.as_ref() {
                // TopExp_Explorer(E, VERTEX) returns the two vertices.
                for v in [&ed.first, &ed.last] {
                    explorer_walk(v, the_to_find, cum, the_out);
                }
            }
        }
        // Solid/Shell: the children are ShapeRef pool indices (architecture
        // limitation above); a Vertex matches nothing below it.
        _ => {}
    }
}

// ---------------------------------------------------------------------------
// BRepGProp::LinearProperties (BRepGProp.cxx L121-165).
// ---------------------------------------------------------------------------

/// OCCT BRepGProp::LinearProperties(S, SProps, SkipShared = false,
/// UseTriangulation = false) — the linear properties (total length, center
/// of mass, moments) of the edges of the shape, accumulated into the
/// GProp_GProps form.  The MeshCinert branch (UseTriangulation, or a
/// non-geometric edge flattened through its polygon) is not translated: an
/// edge without a 3D curve is skipped exactly as the OCCT fall-through of
/// L156-163 (theNodes is null and IsGeom is false).
pub fn linear_properties(the_s: &Shape, the_skip_shared: bool) -> GProps {
    // OCCT L122-124: the origin — gp_Pnt P(0, 0, 0); P.Transform(S.Location()).
    // The pool-free identity-location form of curve_pool_free (topods.rs)
    // applies: the location transform is unreachable without the owning
    // pool's table.
    let p = DVec3::ZERO;
    let mut s_props = GProps::with_location(p);

    // OCCT L126-128: BRepAdaptor_Curve BAC; NCollection_Map anEMap.
    let mut an_emap: HashMap<(u64, u32), ()> = HashMap::new();
    for a_e in explorer_shapes(the_s, ShapeType::Edge) {
        // OCCT L132-135.
        if the_skip_shared
            && an_emap.insert((a_e.ptr_id(), a_e.location), ()).is_some()
        {
            continue;
        }
        // OCCT L138: IsGeom = BRep_Tool::IsGeometric(aE) — the edge has a 3D
        // curve (BRep_Tool.cxx L226-231).
        let Some((a_c0, the_f, the_l)) = crate::topo::topods::curve_pool_free(&a_e) else {
            // OCCT L140: UseTriangulation || !IsGeom -> PreparePolygon
            // (BRepGProp_MeshCinert — GAP); theNodes stays null and L157-162
            // adds nothing.
            continue;
        };
        // OCCT L164-165: BAC.Initialize(aE); BRepGProp_Cinert CG(BAC, P);
        // SProps.Add(CG).
        let (a_curve, a_first, a_last) = edge_adaptor_curve(a_c0, the_f, the_l);
        let mut cg = BRepGPropCinert { base: GProps::new() };
        cg.set_location(p);
        cg.perform(&a_curve, a_first, a_last);
        s_props.add(&cg.base, 1.0);
    }
    s_props
}

// ---------------------------------------------------------------------------
// BRepGProp_Gauss Sinert engine — the Inertia accumulator, the elementary
// part compute, the fixed-order Compute forms and the convert
// (BRepGProp_Gauss.hxx L33-60; .cxx L160-191, L391-424, L418-431,
// L473-491, L1145-1211, L1306-1393).
// ---------------------------------------------------------------------------

/// OCCT BRepGProp_Gauss::Inertia (hxx L33-60) — the inertial moments of the
/// Sinert integrator.
#[derive(Debug, Clone, Copy, Default)]
pub struct Inertia {
    pub mass: f64,
    pub ix: f64,
    pub iy: f64,
    pub iz: f64,
    pub ixx: f64,
    pub iyy: f64,
    pub izz: f64,
    pub ixy: f64,
    pub ixz: f64,
    pub iyz: f64,
}

impl Inertia {
    /// OCCT multAndRestoreInertia (BRepGProp_Gauss.cxx L451-465).
    fn mult_restore(&mut self, the_value: f64) {
        self.mass *= the_value;
        self.ix *= the_value;
        self.iy *= the_value;
        self.iz *= the_value;
        self.ixx *= the_value;
        self.iyy *= the_value;
        self.izz *= the_value;
        self.ixy *= the_value;
        self.ixz *= the_value;
        self.iyz *= the_value;
    }

    /// OCCT addAndRestoreInertia (BRepGProp_Gauss.cxx L433-449).
    fn add_restore(&mut self, the_in: &Inertia) {
        self.mass += the_in.mass;
        self.ix += the_in.ix;
        self.iy += the_in.iy;
        self.iz += the_in.iz;
        self.ixx += the_in.ixx;
        self.iyy += the_in.iyy;
        self.izz += the_in.izz;
        self.ixy += the_in.ixy;
        self.ixz += the_in.ixz;
        self.iyz += the_in.iyz;
    }
}

/// OCCT BRepGProp_Gauss::computeSInertiaOfElementaryPart
/// (BRepGProp_Gauss.cxx L391-424): ds = ||n|| * weight; the static and
/// quadratic moments of the elementary surface part about theLocation.
fn compute_s_inertia_of_elementary_part(
    the_point: DVec3,
    the_normal: DVec3,
    the_location: DVec3,
    the_weight: f64,
    the_out_inertia: &mut Inertia,
) {
    // ds - Jacobien (x, y, z) -> (u, v) = ||n||.
    let ds = the_normal.length() * the_weight;
    let x = the_point.x - the_location.x;
    let y = the_point.y - the_location.y;
    let z = the_point.z - the_location.z;

    the_out_inertia.mass += ds;

    let xds = x * ds;
    let yds = y * ds;
    let zds = z * ds;

    the_out_inertia.ix += xds;
    the_out_inertia.iy += yds;
    the_out_inertia.iz += zds;
    the_out_inertia.ixy += x * yds;
    the_out_inertia.iyz += y * zds;
    the_out_inertia.ixz += x * zds;

    let xxds = x * xds;
    let yyds = y * yds;
    let zzds = z * zds;

    the_out_inertia.ixx += yyds + zzds;
    the_out_inertia.iyy += xxds + zzds;
    the_out_inertia.izz += xxds + yyds;
}

/// OCCT BRepGProp_Gauss::convert (Sinert form, BRepGProp_Gauss.cxx
/// L473-491): the gravity center relative to theLocation, the inertia
/// matrix, and the mass of the accumulated moments.
fn convert_sinert(the_inertia: &Inertia) -> (f64, DVec3, GpMat) {
    let out_mass;
    let out_g;
    if the_inertia.mass.abs() >= EPS_DIM {
        let an_inv_mass = 1.0 / the_inertia.mass;
        out_g = DVec3::new(
            the_inertia.ix * an_inv_mass,
            the_inertia.iy * an_inv_mass,
            the_inertia.iz * an_inv_mass,
        );
        out_mass = the_inertia.mass;
    } else {
        out_mass = 0.0;
        out_g = DVec3::ZERO;
    }
    let out_inertia = GpMat::from_cols(
        DVec3::new(the_inertia.ixx, -the_inertia.ixy, -the_inertia.ixz),
        DVec3::new(-the_inertia.ixy, the_inertia.iyy, -the_inertia.iyz),
        DVec3::new(-the_inertia.ixz, -the_inertia.iyz, the_inertia.izz),
    );
    (out_mass, out_g, out_inertia)
}

/// OCCT BRepGProp_Face::Normal (BRepGProp_Face.cxx L201-210): the point and
/// the D1U x D1V normal at (u, v); mySReverse is applied by the caller
/// (the face occurrence orientation).
fn face_normal(the_surf: &Surface3, u: f64, v: f64, s_reverse: bool) -> (DVec3, DVec3) {
    let (p, du, dv) = the_surf.derivatives(u, v);
    let n = du.cross(dv);
    if s_reverse {
        (p, -n)
    } else {
        (p, n)
    }
}

/// OCCT BRepGProp_Gauss::Compute(BRepGProp_Face, theLocation, out...) — the
/// NATURAL-RESTRICTION fixed-order Sinert integration
/// (BRepGProp_Gauss.cxx L1306-1393, the Sinert branch).  Returns the
/// converted (mass, gravity center relative to theLocation, inertia).
pub fn sinert_compute_natural(the_surf: &Surface3, the_location: DVec3, s_reverse: bool) -> (f64, DVec3, GpMat) {
    // OCCT L1313-1314: the surface natural parameter domain.
    let [lower_u, upper_u, lower_v, upper_v] = the_surf.default_domain();
    // OCCT L1315-1316: checkBounds — an infinite bound switches add/mult to
    // the infinite-aware forms.  For the natural path the kernel surface
    // translation (base::gprop::surface) keeps the evaluated |N| finite for
    // bounded analytic surfaces and converts an infinite domain to mass 0;
    // the same guard applies here (a natural-restriction face is a closed
    // bounded surface).
    if is_infinite_value(lower_u)
        || is_infinite_value(upper_u)
        || is_infinite_value(lower_v)
        || is_infinite_value(upper_v)
    {
        return (0.0, DVec3::ZERO, GpMat::zero());
    }
    // OCCT L1318-1319.
    let (uo0, vo0) = surface_u_v_integration_order(the_surf);
    let u_order = uo0.min(GPM);
    let v_order = vo0.min(GPM);
    // OCCT L1322-1333.
    let (gauss_pu, gauss_wu) = match occt_gauss(u_order) {
        Some(v) => v,
        None => return (0.0, DVec3::ZERO, GpMat::zero()),
    };
    let (gauss_pv, gauss_wv) = match occt_gauss(v_order) {
        Some(v) => v,
        None => return (0.0, DVec3::ZERO, GpMat::zero()),
    };
    // OCCT L1335-1338.
    let um = 0.5 * (upper_u + lower_u);
    let vm = 0.5 * (upper_v + lower_v);
    let ur = 0.5 * (upper_u - lower_u);
    let mut vr = 0.5 * (upper_v - lower_v);

    // OCCT L1343-1382: the V loop with the per-row elementary parts.
    let mut an_inertia = Inertia::default();
    for j in 0..v_order {
        let mut an_inertia_of_elementary_part = Inertia::default();
        let v = vm + vr * gauss_pv[j];
        for i in 0..u_order {
            let a_weight = gauss_wu[i];
            let u = um + ur * gauss_pu[i];
            let (a_point, a_normal) = face_normal(the_surf, u, v, s_reverse);
            // OCCT L1369-1381: the Sinert branch.
            compute_s_inertia_of_elementary_part(a_point, a_normal, the_location, a_weight, &mut an_inertia_of_elementary_part);
        }
        // OCCT L1384-1385.
        an_inertia_of_elementary_part.mult_restore(gauss_wv[j]);
        an_inertia.add_restore(&an_inertia_of_elementary_part);
    }
    // OCCT L1383-1391: vr = vr * ur; the second moments scale with vr.
    vr *= ur;
    an_inertia.ixx *= vr;
    an_inertia.iyy *= vr;
    an_inertia.izz *= vr;
    an_inertia.ixy *= vr;
    an_inertia.ixz *= vr;
    an_inertia.iyz *= vr;
    // OCCT L1384-1392: convert (the Sinert form) then theOutMass *= vr.
    let (mut out_mass, out_g, out_inertia) = convert_sinert(&an_inertia);
    out_mass *= vr;
    (out_mass, out_g, out_inertia)
}

/// OCCT BRepGProp_Gauss::Compute(BRepGProp_Face, BRepGProp_Domain,
/// theLocation, out...) — the boundary (Green-theorem) fixed-order Sinert
/// integration (BRepGProp_Gauss.cxx L1145-1211).  `the_face` is the FACE
/// occurrence whose wires carry the domain edges.  Returns None when a
/// domain-edge Load fails (OCCT L1173-1176 returns with the out-parameters
/// untouched — the caller keeps the PREVIOUS integrator state).
pub fn sinert_compute_domain(
    the_face: &Shape,
    the_surf: &Surface3,
    the_location: DVec3,
    s_reverse: bool,
) -> Option<(f64, DVec3, GpMat)> {
    // OCCT L1153-1154: theSurface.Bounds — the face UV bounds
    // (BRepTools::UVBounds through BRepAdaptor_Surface Restriction=true).
    let [bu1, bu2, bv1, bv2] = face_uv_bounds_pool_free(the_face, the_surf);
    // OCCT L1155-1156: checkBounds — the infinite-aware add/mult forms.
    let inf_bounds = is_infinite_value(bu1)
        || is_infinite_value(bu2)
        || is_infinite_value(bv1)
        || is_infinite_value(bv2);
    let add = |a: f64, b: f64| if inf_bounds { add_inf(a, b) } else { a + b };
    let mult = |a: f64, b: f64| if inf_bounds { mult_inf(a, b) } else { a * b };

    // OCCT L1158-1163.
    let (nu0, nv0) = surface_u_v_integration_order(the_surf);
    let nb_u_gaussgp_pnts = nu0.min(GPM);
    let nb_v_gaussgp_pnts = nv0.min(GPM);
    let nb_gaussgp_pnts = nb_u_gaussgp_pnts.max(nb_v_gaussgp_pnts);
    let (gauss_spv, gauss_swv) = match occt_gauss(nb_gaussgp_pnts) {
        Some(v) => v,
        None => return Some((0.0, DVec3::ZERO, GpMat::zero())),
    };

    // OCCT L1168-1211: the domain loop.
    let mut an_inertia = Inertia::default();
    for a_e in face_domain_edges(the_face) {
        // OCCT L1172-1177: Load(theDomain.Value()) — the edge pcurve; a
        // REVERSED occurrence loads the reversed pcurve and range
        // (BRepGProp_Face.cxx L164-185); a null pcurve fails the Load.
        let (a_c0, a0, b0) = match curve_on_surface_pool_free(&a_e, the_face) {
            Some(v) => v,
            None => match project_curve_on_plane_pool_free(&a_e, the_surf) {
                Some(v) => v,
                None => return None,
            },
        };
        let reversed = a_e.orientation.is_reversed();
        let (l1, l2) = if reversed {
            let x = a0;
            (a_c0.reversed_parameter(b0), a_c0.reversed_parameter(x))
        } else {
            (a0, b0)
        };
        // OCCT L1179: NbCGaussgp_Pnts = min(theSurface.IntegrationOrder(),
        // GPM) — the boundary pcurve Gauss order.
        let nb_c_gaussgp_pnts = curve_integration_order(&a_c0).min(GPM);
        // OCCT L1180: raised to at least the face surface order.
        let nb_c_gaussgp_pnts = nb_c_gaussgp_pnts.max(nb_gaussgp_pnts);
        let (gauss_cp, gauss_cw) = match occt_gauss(nb_c_gaussgp_pnts) {
            Some(v) => v,
            None => return Some((0.0, DVec3::ZERO, GpMat::zero())),
        };
        // OCCT L1183-1186.
        let lm = 0.5 * (l2 + l1);
        let lr = 0.5 * (l2 - l1);

        let mut a_c_inertia = Inertia::default();
        for i in 0..nb_c_gaussgp_pnts {
            let l = lm + lr * gauss_cp[i];

            // OCCT L1190-1196: D12d(l, Puv, Vuv) — the loaded (possibly
            // reversed) pcurve point and derivative.
            let (puv, vuv) = if reversed {
                let rp = a_c0.reversed_parameter(l);
                (a_c0.point_at(rp), -a_c0.derivative_at(rp))
            } else {
                (a_c0.point_at(l), a_c0.derivative_at(l))
            };

            let v = puv.y;
            let u2 = puv.x;

            // OCCT L1198-1201.
            let dul = vuv.y * gauss_cw[i];
            let um = 0.5 * (u2 + bu1);
            let ur = 0.5 * (u2 - bu1);

            // OCCT L1203-1214.
            let mut a_local_inertia = Inertia::default();
            for j in 0..nb_gaussgp_pnts {
                let u = add(um, mult(ur, gauss_spv[j]));
                let a_weight = dul * gauss_swv[j];

                let (a_point, a_normal) = face_normal(the_surf, u, v, s_reverse);

                compute_s_inertia_of_elementary_part(a_point, a_normal, the_location, a_weight, &mut a_local_inertia);
            }

            // OCCT L1216-1217.
            a_local_inertia.mult_restore(ur);
            a_c_inertia.add_restore(&a_local_inertia);
        }

        // OCCT L1219-1220.
        a_c_inertia.mult_restore(lr);
        an_inertia.add_restore(&a_c_inertia);
    }

    // OCCT L1213: convert (the Sinert form).
    Some(convert_sinert(&an_inertia))
}

// ---------------------------------------------------------------------------
// The infinite-aware arithmetic of BRepGProp_Gauss (BRepGProp_Gauss.cxx
// L36-155) is shared with the kernel surface translation
// (base::gprop::surface::add_inf / mult_inf, the faithful AddInf/MultInf
// forms switched in by checkBounds L418-431).
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Pool-free BRepGProp_Face / BRepTools::UVBounds reads.
// ---------------------------------------------------------------------------

/// OCCT TopExp_Explorer over the face wires — the domain edges of the face
/// with the composed orientation (BRepGProp_Domain::Next skips the
/// INTERNAL/EXTERNAL occurrences, BRepGProp_Domain.cxx L27-38).
fn face_domain_edges(the_face: &Shape) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    if let TShape::Face(fd) = the_face.data.as_ref() {
        for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
            if w.is_null() {
                continue;
            }
            if let TShape::Wire(wd) = w.data.as_ref() {
                for e in &wd.edges {
                    let ori = the_face.orientation.compose(w.orientation.compose(e.orientation));
                    // BRepGProp_Domain::Next: INTERNAL/EXTERNAL edges are not
                    // boundary edges.
                    if matches!(ori, Orientation::Internal | Orientation::External) {
                        continue;
                    }
                    let mut occ = e.clone();
                    occ.orientation = ori;
                    out.push(occ);
                }
            }
        }
    }
    out
}

/// OCCT BRepTools::UVBounds (BRepTools.cxx L64-80 + AddUVBounds L126-367) —
/// the pool-free form of base::gprop::surface::face_uv_bounds: the union of
/// every boundary edge pcurve 2D box, clamped to the natural domain for
/// non-periodic surfaces; a face without pcurves falls back to the surface
/// natural bounds (L139-153).  Returns [u_min, u_max, v_min, v_max].
pub fn face_uv_bounds_pool_free(the_face: &Shape, the_surf: &Surface3) -> [f64; 4] {
    // aS->Bounds (L191-192): the surface natural parameter domain.
    let [a_u_min, a_u_max, a_v_min, a_v_max] = the_surf.default_domain();
    // Bnd_Box2d aBoxS (L176).
    let mut b = [f64::INFINITY, f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY];
    let mut any = false;
    for an_edge in face_domain_edges(the_face) {
        // BRep_Tool::CurveOnSurface (L179) with the CurveOnPlane fallback
        // (L180-183 / BRep_Tool.cxx L367-372).
        let (c, a_t1, a_t2) = match curve_on_surface_pool_free(&an_edge, the_face) {
            Some(v) => v,
            None => match project_curve_on_plane_pool_free(&an_edge, the_surf) {
                Some(v) => v,
                None => continue,
            },
        };
        // BndLib_Add2dCurve::Add (L185): the 2D bounding box.
        let a_box_c = curve2d_bounding_box(&c, a_t1, a_t2, 0.0);
        let mut a_x_min = a_box_c[0];
        let mut a_x_max = a_box_c[1];
        let mut a_y_min = a_box_c[2];
        let mut a_y_max = a_box_c[3];
        // U periodicity (L202-281): clamp for a non-periodic surface
        // (L270-280).
        if !the_surf.is_u_periodic() {
            if (a_x_min < a_u_min) && (a_u_min < a_x_max) {
                a_x_min = a_u_min;
            }
            if (a_x_min < a_u_max) && (a_u_max < a_x_max) {
                a_x_max = a_u_max;
            }
        }
        // V periodicity (L283-362, L351-361).
        if !the_surf.is_v_periodic() {
            if (a_y_min < a_v_min) && (a_v_min < a_y_max) {
                a_y_min = a_v_min;
            }
            if (a_y_min < a_v_max) && (a_v_max < a_y_max) {
                a_y_max = a_v_max;
            }
        }
        // aBoxS.Update + aB.Add (L364-366).
        b[0] = b[0].min(a_x_min);
        b[1] = b[1].max(a_x_max);
        b[2] = b[2].min(a_y_min);
        b[3] = b[3].max(a_y_max);
        any = true;
    }
    if !any {
        // UVBounds: void box -> the surface natural bounds (L139-153).
        return [a_u_min, a_u_max, a_v_min, a_v_max];
    }
    b
}

/// OCCT BRep_Tool::CurveOnPlane (BRep_Tool.cxx L367-372 -> CurveOnPlane) —
/// the pool-free form of base::gprop::surface::project_curve_on_plane: the
/// edge 3D curve projected onto the PLANAR face along the plane normal
/// (ProjLib_Plane::Project keeps the curve type: a circle projects to the
/// same-parameter 2D circle, a line to the 2D line through the projected
/// frame coordinates; the 2D parameterization keeps the 3D parameter
/// range).  Architecture limitation: the identity-location form — the
/// located-instance transform is unreachable without the owning pool.
fn project_curve_on_plane_pool_free(
    the_edge: &Shape,
    the_surf: &Surface3,
) -> Option<(crate::geom::Curve2d, f64, f64)> {
    use crate::geom::CurveEval as _;
    // BRep_Tool::CurveOnPlane applies only to a planar surface (the OCCT
    // Plane() read of CurveOnPlane); other surface types return a null
    // pcurve — the Load fails.
    let crate::geom::Surface3::Plane(pl) = the_surf else {
        return None;
    };
    let (c3, r3_1, r3_2) = crate::topo::topods::curve_pool_free(the_edge)?;
    let r3 = [r3_1, r3_2];
    let o = pl.origin;
    if let crate::geom::Curve3::Circle(c) = &c3 {
        // ProjLib_Plane::Project(Circle): the center maps onto the plane
        // frame; the in-plane axes of the circle give the 2D image axes.
        let d = c.center - o;
        let c2 = glam::DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir));
        let vx = glam::DVec2::new(
            c.radius * c.x_dir.dot(pl.u_dir),
            c.radius * c.x_dir.dot(pl.v_dir),
        );
        let vy = glam::DVec2::new(
            c.radius * c.y_dir.dot(pl.u_dir),
            c.radius * c.y_dir.dot(pl.v_dir),
        );
        let lx = vx.length();
        let ly = vy.length();
        if lx > 1e-12
            && ly > 1e-12
            && vx.dot(vy).abs() < 1e-12 * lx * ly
            && (lx - ly).abs() < 1e-9 * lx.max(1.0)
        {
            let c2d = crate::geom::Curve2d::Circle(crate::geom::Circle2d {
                center: c2,
                x_dir: vx / lx,
                y_dir: vy / ly,
                radius: (lx + ly) * 0.5,
            });
            return Some((c2d, r3[0], r3[1]));
        }
        return None;
    }
    // ProjLib_Plane::Project(Line): the 2D line (u, v) = ((p - O).U, (p - O).V).
    let p0 = c3.point_at(r3[0]);
    let p1 = c3.point_at(r3[1]);
    let (u0, v0) = ((p0 - o).dot(pl.u_dir), (p0 - o).dot(pl.v_dir));
    let (u1, v1) = ((p1 - o).dot(pl.u_dir), (p1 - o).dot(pl.v_dir));
    let du = u1 - u0;
    let dv = v1 - v0;
    let len = (du * du + dv * dv).sqrt();
    if len < 1e-30 {
        return None;
    }
    let c = crate::geom::Curve2d::Line(crate::geom::Line2d::new(
        glam::DVec2::new(u0, v0),
        glam::DVec2::new(du / len, dv / len),
    ));
    Some((c, 0.0, len))
}

// ---------------------------------------------------------------------------
// BRepGProp::SurfaceProperties (BRepGProp.cxx L267-279 -> the static
// surfaceProperties L172-266).
// ---------------------------------------------------------------------------

/// OCCT roughBaryCenter(S) (BRepGProp.cxx L90-108): the average of the
/// vertex points; the triangulation fallback (L98-107) is not reachable in
/// the pool-free form (rcad carries no Poly_Triangulation here) and the
/// origin is returned when the shape has no vertices (OCCT xyz(0,0,0)).
fn rough_bary_center(the_s: &Shape) -> DVec3 {
    // OCCT L93-99: the average of TopExp_Explorer(S, TopAbs_VERTEX) points.
    let mut xyz = DVec3::ZERO;
    let mut i = 0usize;
    for a_v in explorer_shapes(the_s, ShapeType::Vertex) {
        if let TShape::Vertex(vd) = a_v.data.as_ref() {
            xyz += vd.point;
            i += 1;
        }
    }
    if i > 0 {
        xyz / i as f64
    } else {
        xyz
    }
}

/// OCCT BRepGProp::SurfaceProperties(S, Props, SkipShared = false,
/// UseTriangulation = false) (BRepGProp.cxx L267-279): the surface
/// properties of the faces of the shape accumulated into the GProp_GProps
/// form through the Eps = 1.0 FIXED-ORDER branch of surfaceProperties
/// (L231-253 — G.Perform(BF) for natural-restriction faces,
/// G.Perform(BF, BD) otherwise).  The Eps < 1.0 adaptive branch
/// (L231-242) and the triangulation branch (L216-222) are not translated
/// with moments (the kernel checkprops path, base::gprop::surface, is the
/// mass-only Sinert).
pub fn surface_properties(the_s: &Shape, the_skip_shared: bool) -> GProps {
    // OCCT L273-274: the origin — the identity-location form.
    let p = DVec3::ZERO;
    let mut props = GProps::with_location(p);
    surface_properties_eps(the_s, &mut props, 1.0, the_skip_shared);
    props
}

/// OCCT static surfaceProperties(S, Props, Eps, SkipShared,
/// UseTriangulation) (BRepGProp.cxx L172-266) — the Eps >= 1.0 fixed-order
/// branch.  Returns ErrorMax (0.0 on the fixed-order path: OCCT only fills
/// it from the adaptive G.Perform epsilon, L236-241).
fn surface_properties_eps(
    the_s: &Shape,
    the_props: &mut GProps,
    the_eps: f64,
    the_skip_shared: bool,
) -> f64 {
    // OCCT L177-180: the persistent integrator at the rough barycenter.
    let p = rough_bary_center(the_s);
    let _ = the_eps; // the adaptive branch is not translated (see above)
    // OCCT L188: NCollection_Map aFMap.
    let mut a_fmap: HashMap<(u64, u32), ()> = HashMap::new();
    // The reused BRepGProp_Sinert state G (L178-180): a face whose domain
    // Load fails leaves dim/g/inertia at the PREVIOUS face's converted
    // values, which Props.Add(G) re-adds (OCCT L1173-1176 + L254).
    let mut prev: (f64, DVec3, GpMat) = (0.0, DVec3::ZERO, GpMat::zero());

    for a_f in explorer_shapes(the_s, ShapeType::Face) {
        // OCCT L191-194.
        if the_skip_shared && a_fmap.insert((a_f.ptr_id(), a_f.location), ()).is_some() {
            continue;
        }
        // OCCT L197-213: NoSurf / NoTri — rcad carries no Poly_Triangulation
        // (UseTriangulation = false, NoTri = true), so NoSurf && NoTri skips
        // the face (L211-214) and the Gauss branch (L224+) needs a surface.
        let Some(a_surf) = (match a_f.data.as_ref() {
            TShape::Face(fd) => fd.surface.clone(),
            _ => None,
        }) else {
            continue;
        };
        // OCCT L226: IsNatRestr = (F.NbChildren() == 0) — no wires.
        let is_nat_restr = match a_f.data.as_ref() {
            TShape::Face(fd) => fd.outer_wire.is_null() && fd.inner_wires.is_empty(),
            _ => true,
        };
        // BRepGProp_Face::Load(F) L196: mySReverse = (F.Orientation() ==
        // TopAbs_REVERSED).
        let s_reverse = a_f.orientation.is_reversed();
        // OCCT L229-253: the Eps >= 1.0 fixed-order branch.
        let converted = if is_nat_restr {
            sinert_compute_natural(&a_surf, p, s_reverse)
        } else {
            match sinert_compute_domain(&a_f, &a_surf, p, s_reverse) {
                Some(v) => v,
                None => prev, // the Load-failed early return keeps G's state
            }
        };
        prev = converted;
        // OCCT L254: Props.Add(G) — the G state as a GProp_GProps at loc P.
        let mut g_item = GProps::with_location(p);
        g_item.dim = converted.0;
        g_item.g = converted.1;
        g_item.inertia = converted.2;
        the_props.add(&g_item, 1.0);
    }
    // OCCT L263: ErrorMax (the fixed-order path leaves it 0).
    0.0
}

// ---------------------------------------------------------------------------
// BRepGProp::VolumeProperties (BRepGProp.cxx L555-586) — the aggregation of
// the existing per-face Vinert components.
// ---------------------------------------------------------------------------

/// OCCT BRepGProp::VolumeProperties(S, VProps, OnlyClosed = false,
/// SkipShared = false, Eps = 1.0) — the pool (BRep) form: every face
/// occurrence is integrated with BRepGProp_Vinert about the shape origin
/// (aCoeff = {0, 0, 0}, isByPoint = true) and combined with
/// GProp_GProps::Add.
///
/// The accumulation formula mirrored line-by-line is the existing
/// `shape_vinert` (base::gprop::volume): BRepGProp.cxx L298-409
/// volumePropertiesFaces per-face loop + BRepGProp_Gauss.cxx L1215-1302
/// (Vinert domain Compute) / L1306-1393 (natural Compute) +
/// computeVInertiaOfElementaryPart isByPoint L306-339 — with the cumulative
/// TopExp_Explorer orientation handled by VinertFace::add/sub
/// (BRepGProp_Face::Load L196 mySReverse).  What this driver adds is the
/// GProp_GProps combination the OCCT flow applies per face
/// (GProp_GProps::Add L48-64 first branch — every face shares the location
/// P = the origin, so the item-wise adds collapse into the totals):
/// - dim = the accumulated Mass,
/// - g = (Ix, Iy, Iz) / dim — BRepGProp_Gauss::convert L498-511 with
///   aCoeff = {0, 0, 0} (isByPoint),
/// - inertia = the Vinert convert matrix L515-519 — the CONVERTED form
///   keeps the accumulated SIGN of the Vinert moments: the VinertFace
///   components are Ixx = +∫(y²+z²)·..., Ixy = −∫xy·... and the matrix
///   columns are (Ixx, Ixy, Ixz), (Ixy, Iyy, Iyz), (Ixz, Iyz, Izz).
///
/// The OnlyClosed / SkipShared / Eps != 1.0 variants (L590-725,
/// VolumePropertiesGK L727-940) are not translated — no consumer.
pub fn volume_properties(brep: &crate::topo::topods::BRep) -> GProps {
    // OCCT L560-561: the origin — gp_Pnt P(0, 0, 0); P.Transform(S.Location());
    // VProps = GProp_GProps(P).  The kernel Vinert integrates about the
    // world origin (the location form above is the identity-location one).
    let p = DVec3::ZERO;
    let mut v_props = GProps::with_location(p);
    // OCCT L566-586 -> volumePropertiesFaces L298-409: the per-face Vinert
    // accumulation (shape_vinert) + the GProp_GProps::Add combination.
    let v = crate::base::gprop::volume::shape_vinert(brep);
    // BRepGProp_Gauss::convert (L493-532, Vinert isByPoint aCoeff = {0,0,0}).
    v_props.dim = v.mass;
    if v.mass.abs() >= EPS_DIM {
        let an_inv_mass = 1.0 / v.mass;
        v_props.g = DVec3::new(
            v.ix * an_inv_mass,
            v.iy * an_inv_mass,
            v.iz * an_inv_mass,
        );
    } else {
        v_props.g = DVec3::ZERO;
        v_props.dim = 0.0;
    }
    v_props.inertia = GpMat::from_cols(
        DVec3::new(v.ixx, v.ixy, v.ixz),
        DVec3::new(v.ixy, v.iyy, v.iyz),
        DVec3::new(v.ixz, v.iyz, v.izz),
    );
    v_props
}
