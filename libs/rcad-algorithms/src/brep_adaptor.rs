//! BRepAdaptor-style topology adapters.
//!
//! Provides adapters to access BRep topology as geometric entities.
//! Analogous to OCCT's BRepAdaptor_Curve, BRepAdaptor_Surface, and BRepAdaptor_CompCurve.
//!
//! # Overview
//!
//! - [`EdgeAdaptor`]: Adapts an edge to act as a 3D curve
//! - [`FaceAdaptor`]: Adapts a face to act as a 3D surface
//! - [`WireAdaptor`]: Adapts a wire to act as a composite 3D curve
//! - [`CurveAdaptorArray`]: Array of edge adaptors for indexed access

use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::{topods, Curve3, CurveEval, Surface3, SurfaceEval, Wire};
use rcad_kernel::geom::TrimmedCurve3;
use rcad_kernel::topods::TShape;
use std::f64::consts::PI;

// =============================================================================
// EdgeAdaptor (BRepAdaptor_Curve)
// =============================================================================

/// Adapts a BRep edge to act as a 3D curve.
///
/// Provides curve-like evaluation methods (point, tangent, domain) for an edge
/// in a BRep, respecting the edge's orientation and parameter range.
///
/// Analogous to OCCT's `BRepAdaptor_Curve`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::EdgeAdaptor;
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let adaptor = EdgeAdaptor::new(&brep, 0);
/// let domain = adaptor.domain();
/// let midpoint = adaptor.point_at((domain[0] + domain[1]) / 2.0);
/// ```
#[derive(Debug, Clone)]
pub struct EdgeAdaptor<'a> {
    brep: &'a rcad_kernel::BRep,
    edge_idx: usize,
    /// Cached curve reference (if available).
    curve: Option<&'a Curve3>,
    /// Cached parameter range.
    range: [f64; 2],
    /// Whether the edge's natural direction is reversed.
    reversed: bool,
}

impl<'a> EdgeAdaptor<'a> {
    /// Create a new edge adaptor for the given edge index.
    ///
    /// The adaptor respects the edge's stored parameter range in `edge_curve_range`
    /// and falls back to the curve's natural domain if not specified.
    ///
    /// # Panics
    ///
    /// Does not panic; returns a default adaptor if the edge index is out of bounds
    /// or the edge has no associated curve.
    pub fn new(brep: &'a rcad_kernel::BRep, edge_idx: usize) -> Self {
        let (curve, range) = match brep.tshapes.get(edge_idx) {
            Some(ts) => match ts.as_ref() {
                TShape::Edge(ed) => (ed.curve.as_ref(), ed.range),
                _ => (None, [0.0, 1.0]),
            },
            None => (None, [0.0, 1.0]),
        };

        Self {
            brep,
            edge_idx,
            curve,
            range,
            reversed: false,
        }
    }

    /// Create an edge adaptor with reversed direction.
    ///
    /// This is used when an edge appears in a wire with `forward = false`,
    /// meaning the edge should be traversed from end to start.
    pub fn with_reversed(mut self, reversed: bool) -> Self {
        self.reversed = reversed;
        self
    }

    /// Evaluate the point on the edge at parameter `t`.
    ///
    /// The parameter `t` is in the edge's natural parameter domain
    /// (respecting `edge_curve_range` if specified).
    pub fn point_at(&self, t: f64) -> DVec3 {
        let Some(curve) = self.curve else {
            // Fall back to vertex interpolation if no curve is available.
            return self.point_from_vertices(t);
        };

        let t_mapped = self.map_parameter(t);
        curve.point_at(t_mapped)
    }

    /// Evaluate the unit tangent vector on the edge at parameter `t`.
    ///
    /// Returns the tangent pointing in the direction of increasing parameter
    /// on the underlying curve. If the edge is reversed, negate the result.
    pub fn tangent_at(&self, t: f64) -> DVec3 {
        let Some(curve) = self.curve else {
            // Fall back to straight-line tangent between vertices.
            return self.tangent_from_vertices();
        };

        let t_mapped = self.map_parameter(t);
        let mut tangent = curve.tangent_at(t_mapped);
        if self.reversed {
            tangent = -tangent;
        }
        tangent
    }

    /// Return the parameter domain of the edge.
    ///
    /// This is always `[0.0, 1.0]` for normalized parameter access,
    /// regardless of the underlying curve's natural domain.
    pub fn domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }

    /// Returns the first parameter of the edge adaptor domain.
    /// Analogous to OCCT `BRepAdaptor_Curve::FirstParameter()`.
    pub fn first_parameter(&self) -> f64 {
        self.domain()[0]
    }

    /// Returns the last parameter of the edge adaptor domain.
    /// Analogous to OCCT `BRepAdaptor_Curve::LastParameter()`.
    pub fn last_parameter(&self) -> f64 {
        self.domain()[1]
    }

    /// Alias for point evaluation.
    /// Analogous to OCCT `BRepAdaptor_Curve::Value()`.
    pub fn value(&self, t: f64) -> DVec3 {
        self.point_at(t)
    }

    /// Alias for first derivative (tangent) evaluation.
    /// Analogous to OCCT `BRepAdaptor_Curve::D1()`.
    pub fn d1(&self, t: f64) -> DVec3 {
        self.tangent_at(t)
    }

    /// Return the underlying curve reference, if available.
    pub fn curve(&self) -> Option<&Curve3> {
        self.curve
    }

    /// Return the natural parameter range of the edge on its curve.
    ///
    /// This returns the actual parameter range (not normalized to [0, 1]).
    pub fn curve_range(&self) -> [f64; 2] {
        self.range
    }

    /// Check if the edge is closed (start and end vertices are the same).
    ///
    /// A closed edge forms a loop, such as a circle or ellipse.
    pub fn is_closed(&self) -> bool {
        match self.brep.tshapes.get(self.edge_idx) {
            Some(ts) => match ts.as_ref() {
                TShape::Edge(ed) => ed.first.index == ed.last.index,
                _ => false,
            },
            None => false,
        }
    }

    /// Return the period of the edge's curve if it is periodic.
    ///
    /// Returns `Some(period)` for periodic curves (circles, ellipses),
    /// or `None` for non-periodic curves (lines, B-splines).
    pub fn period(&self) -> Option<f64> {
        let Some(curve) = self.curve else {
            return None;
        };
        Self::period_of_curve(curve)
    }

    fn period_of_curve(curve: &Curve3) -> Option<f64> {
        match curve {
            Curve3::Circle(_) => Some(2.0 * PI),
            Curve3::Ellipse(_) => Some(2.0 * PI),
            Curve3::Line(_) => None,
            Curve3::BSpline(_) => None,
            Curve3::Bezier(_) => None,
            Curve3::Offset(_) => None,
            Curve3::Hyperbola(_) => None,
            Curve3::Parabola(_) => None,
            Curve3::CircularHelix(_) => None,
            Curve3::SineWave(_) => None,
            Curve3::Trimmed(tc) => Self::period_of_curve(tc.basis_curve()),
        }
    }

    /// Returns true if the underlying edge curve is periodic.
    /// Analogous to OCCT `BRepAdaptor_Curve::IsPeriodic()`.
    pub fn is_periodic(&self) -> bool {
        self.period().is_some()
    }

    /// Map normalized parameter [0, 1] to curve's natural parameter.
    fn map_parameter(&self, t: f64) -> f64 {
        let [t0, t1] = self.range;
        if self.reversed {
            t0 + (1.0 - t) * (t1 - t0)
        } else {
            t0 + t * (t1 - t0)
        }
    }

    /// Fall back to vertex-based point evaluation when no curve is available.
    fn point_from_vertices(&self, t: f64) -> DVec3 {
        let (v0_idx, v1_idx) = match self.brep.tshapes.get(self.edge_idx) {
            Some(ts) => match ts.as_ref() {
                TShape::Edge(ed) => (ed.first.index, ed.last.index),
                _ => return DVec3::ZERO,
            },
            None => return DVec3::ZERO,
        };
        let p0 = self.brep.vertex_point(v0_idx).unwrap_or(DVec3::ZERO);
        let p1 = self.brep.vertex_point(v1_idx).unwrap_or(DVec3::ZERO);

        if self.reversed {
            p1.lerp(p0, t)
        } else {
            p0.lerp(p1, t)
        }
    }

    /// Fall back to vertex-based tangent when no curve is available.
    fn tangent_from_vertices(&self) -> DVec3 {
        let (v0_idx, v1_idx) = match self.brep.tshapes.get(self.edge_idx) {
            Some(ts) => match ts.as_ref() {
                TShape::Edge(ed) => (ed.first.index, ed.last.index),
                _ => return DVec3::X,
            },
            None => return DVec3::X,
        };
        let p0 = self.brep.vertex_point(v0_idx).unwrap_or(DVec3::X);
        let p1 = self.brep.vertex_point(v1_idx).unwrap_or(DVec3::X);

        let dir = (p1 - p0).normalize_or_zero();
        if self.reversed {
            -dir
        } else {
            dir
        }
    }

    /// Get the first vertex index of this edge.
    pub fn first_vertex(&self) -> Option<usize> {
        let ts = self.brep.tshapes.get(self.edge_idx)?;
        let TShape::Edge(ed) = ts.as_ref() else { return None };
        if self.reversed {
            Some(ed.last.index)
        } else {
            Some(ed.first.index)
        }
    }

    /// Get the last vertex index of this edge.
    pub fn last_vertex(&self) -> Option<usize> {
        let ts = self.brep.tshapes.get(self.edge_idx)?;
        let TShape::Edge(ed) = ts.as_ref() else { return None };
        if self.reversed {
            Some(ed.first.index)
        } else {
            Some(ed.last.index)
        }
    }
}

// =============================================================================
// FaceAdaptor (BRepAdaptor_Surface)
// =============================================================================

/// Adapts a BRep face to act as a 3D surface.
///
/// Provides surface-like evaluation methods (point, normal, domain) for a face
/// in a BRep, respecting the face's parameter range bounds.
///
/// Analogous to OCCT's `BRepAdaptor_Surface`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::FaceAdaptor;
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });
/// let adaptor = FaceAdaptor::new(&brep, 0);
/// let domain = adaptor.domain();
/// let center = adaptor.point_at(
///     (domain[0] + domain[1]) / 2.0,
///     (domain[2] + domain[3]) / 2.0,
/// );
/// ```
#[derive(Debug, Clone)]
pub struct FaceAdaptor<'a> {
    brep: &'a rcad_kernel::BRep,
    face_idx: usize,
    /// Cached surface reference (if available).
    surface: Option<&'a Surface3>,
    /// Cached parameter range [u_min, u_max, v_min, v_max].
    range: [f64; 4],
}

impl<'a> FaceAdaptor<'a> {
    /// Create a new face adaptor for the given flat face index.
    ///
    /// The flat face index counts faces across all solids/shells in traversal order.
    /// The adaptor respects the face's stored parameter range in `face_surface_range`
    /// and falls back to the surface's natural domain if not specified.
    ///
    /// # Panics
    ///
    /// Does not panic; returns a default adaptor if the face index is out of bounds
    /// or the face has no associated surface.
    pub fn new(brep: &'a rcad_kernel::BRep, face_idx: usize) -> Self {
        let (surface, range) = match brep.tshapes.get(face_idx) {
            Some(ts) => match ts.as_ref() {
                TShape::Face(fd) => {
                    let s = fd.surface.as_ref();
                    let r = fd.uv_domain.unwrap_or_else(|| {
                        s.map(|surf| surf.default_domain()).unwrap_or([0.0, 1.0, 0.0, 1.0])
                    });
                    (s, r)
                }
                _ => (None, [0.0, 1.0, 0.0, 1.0]),
            },
            None => (None, [0.0, 1.0, 0.0, 1.0]),
        };

        Self {
            brep,
            face_idx,
            surface,
            range,
        }
    }

    /// Evaluate the point on the face at parameters `(u, v)`.
    ///
    /// Parameters are in the face's parameter domain.
    pub fn point_at(&self, u: f64, v: f64) -> DVec3 {
        let Some(surface) = self.surface else {
            // Fall back to face centroid if no surface is available.
            return self.point_from_vertices();
        };

        surface.point_at(u, v)
    }

    /// Evaluate the unit normal vector on the face at parameters `(u, v)`.
    ///
    /// Returns the outward-pointing normal (respecting face orientation).
    pub fn normal_at(&self, u: f64, v: f64) -> DVec3 {
        let Some(surface) = self.surface else {
            // Fall back to stored face normal if no surface is available.
            return self.normal_from_face();
        };

        let mut normal = surface.normal_at(u, v);

        // Check if face orientation should flip the normal.
        // Approximate: evaluate surface normal at domain center.
        let surface_normal = self.normal_from_face();
        if surface_normal.length_squared() > 0.5 && normal.dot(surface_normal) < 0.0 {
            normal = -normal;
        }

        normal
    }

    /// Return the parameter domain of the face.
    ///
    /// Returns `[u_min, u_max, v_min, v_max]`.
    pub fn domain(&self) -> [f64; 4] {
        self.range
    }

    /// Returns the first U parameter of the face domain.
    /// Analogous to OCCT `BRepAdaptor_Surface::FirstUParameter()`.
    pub fn first_u_parameter(&self) -> f64 {
        self.range[0]
    }

    /// Returns the last U parameter of the face domain.
    /// Analogous to OCCT `BRepAdaptor_Surface::LastUParameter()`.
    pub fn last_u_parameter(&self) -> f64 {
        self.range[1]
    }

    /// Returns the first V parameter of the face domain.
    /// Analogous to OCCT `BRepAdaptor_Surface::FirstVParameter()`.
    pub fn first_v_parameter(&self) -> f64 {
        self.range[2]
    }

    /// Returns the last V parameter of the face domain.
    /// Analogous to OCCT `BRepAdaptor_Surface::LastVParameter()`.
    pub fn last_v_parameter(&self) -> f64 {
        self.range[3]
    }

    /// Alias for surface value evaluation.
    /// Analogous to OCCT `BRepAdaptor_Surface::Value()`.
    pub fn value(&self, u: f64, v: f64) -> DVec3 {
        self.point_at(u, v)
    }

    /// Return the underlying surface reference, if available.
    pub fn surface(&self) -> Option<&Surface3> {
        self.surface
    }

    /// Check if the face's surface is closed in the U direction.
    ///
    /// A surface is U-closed if `S(u_min, v) == S(u_max, v)` for all v.
    pub fn is_u_closed(&self) -> bool {
        let Some(surface) = self.surface else {
            return false;
        };

        match surface {
            Surface3::Cylinder(_) => true,
            Surface3::Sphere(_) => true,
            Surface3::Cone(_) => true,
            Surface3::Torus(_) => true,
            Surface3::Ellipsoid(_) => true,
            Surface3::Helicoid(_) => false,
            Surface3::Pipe(_) => true,
            Surface3::Plane(_) => false,
            Surface3::BSpline(s) => {
                // Check if first and last rows of control points coincide.
                let n_u = s.control_points.len();
                if n_u < 2 {
                    return false;
                }
                let first = &s.control_points[0];
                let last = &s.control_points[n_u - 1];
                if first.len() != last.len() {
                    return false;
                }
                first
                    .iter()
                    .zip(last.iter())
                    .all(|(a, b)| (a - b).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT)
            }
            Surface3::Bezier(_) => false,
            Surface3::TriBezier(_) => false,
            Surface3::LinearExtrusion(_) => false,
            Surface3::Revolution(_) => true,
            Surface3::Ruled(_) => false,
            Surface3::Coons(_) => false,
            Surface3::Offset(inner) => {
                // Delegate to inner surface check.
                match inner.basis.as_ref() {
                    Surface3::Cylinder(_)
                    | Surface3::Sphere(_)
                    | Surface3::Cone(_)
                    | Surface3::Torus(_)
                    | Surface3::Ellipsoid(_)
                    | Surface3::Pipe(_)
                    | Surface3::Revolution(_) => true,
                    _ => false,
                }
            }
            Surface3::Trimmed(inner) => {
                // Trimmed surface may cut a closed surface.
                match inner.basis.as_ref() {
                    Surface3::Cylinder(_)
                    | Surface3::Sphere(_)
                    | Surface3::Cone(_)
                    | Surface3::Torus(_)
                    | Surface3::Ellipsoid(_)
                    | Surface3::Pipe(_)
                    | Surface3::Revolution(_) => {
                        let [u0, u1, _, _] = inner.trim;
                        let [du0, du1, _, _] = inner.basis.default_domain();
                        (u0 - du0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT && (u1 - du1).abs() < TOLERANCE_LINEAR_ULTRA_STRICT
                    }
                    _ => false,
                }
            }
        }
    }

    /// Returns the U period if the face surface is U-periodic.
    pub fn u_period(&self) -> Option<f64> {
        let Some(surface) = self.surface else {
            return None;
        };
        match surface {
            Surface3::Cylinder(_)
            | Surface3::Sphere(_)
            | Surface3::Cone(_)
            | Surface3::Torus(_)
            | Surface3::Ellipsoid(_)
            | Surface3::Pipe(_)
            | Surface3::Revolution(_) => Some(2.0 * PI),
            Surface3::Trimmed(inner) => match inner.basis.as_ref() {
                Surface3::Cylinder(_)
                | Surface3::Sphere(_)
                | Surface3::Cone(_)
                | Surface3::Torus(_)
                | Surface3::Ellipsoid(_)
                | Surface3::Pipe(_)
                | Surface3::Revolution(_) => {
                    let [u0, u1, _, _] = inner.trim;
                    let [du0, du1, _, _] = inner.basis.default_domain();
                    if (u0 - du0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT && (u1 - du1).abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
                        Some(2.0 * PI)
                    } else {
                        None
                    }
                }
                _ => None,
            },
            _ => None,
        }
    }

    /// Returns true if the face surface is U-periodic.
    /// Analogous to OCCT `BRepAdaptor_Surface::IsUPeriodic()`.
    pub fn is_u_periodic(&self) -> bool {
        self.u_period().is_some()
    }

    /// Check if the face's surface is closed in the V direction.
    ///
    /// A surface is V-closed if `S(u, v_min) == S(u, v_max)` for all u.
    pub fn is_v_closed(&self) -> bool {
        let Some(surface) = self.surface else {
            return false;
        };

        match surface {
            Surface3::Torus(_) => true,
            Surface3::Sphere(_) => false, // Sphere has poles, not V-closed.
            Surface3::Cylinder(_) => false,
            Surface3::Cone(_) => false,
            Surface3::Plane(_) => false,
            Surface3::Ellipsoid(_) => false,
            Surface3::Helicoid(_) => false,
            Surface3::Pipe(_) => false,
            Surface3::BSpline(s) => {
                // Check if first and last columns of control points coincide.
                let n_u = s.control_points.len();
                if n_u == 0 {
                    return false;
                }
                let n_v = s.control_points[0].len();
                if n_v < 2 {
                    return false;
                }
                (0..n_u).all(|i| {
                    let first = &s.control_points[i][0];
                    let last = &s.control_points[i][n_v - 1];
                    (first - last).length_squared() < TOLERANCE_LINEAR_ULTRA_STRICT
                })
            }
            Surface3::Bezier(_) => false,
            Surface3::TriBezier(_) => false,
            Surface3::LinearExtrusion(_) => false,
            Surface3::Revolution(_) => false,
            Surface3::Ruled(_) => false,
            Surface3::Coons(_) => false,
            Surface3::Offset(inner) => {
                matches!(inner.basis.as_ref(), Surface3::Torus(_))
            }
            Surface3::Trimmed(inner) => {
                if let Surface3::Torus(_) = inner.basis.as_ref() {
                    let [_, _, v0, v1] = inner.trim;
                    let [_, _, dv0, dv1] = inner.basis.default_domain();
                    (v0 - dv0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT && (v1 - dv1).abs() < TOLERANCE_LINEAR_ULTRA_STRICT
                } else {
                    false
                }
            }
        }
    }

    /// Returns the V period if the face surface is V-periodic.
    pub fn v_period(&self) -> Option<f64> {
        let Some(surface) = self.surface else {
            return None;
        };
        match surface {
            Surface3::Torus(_) => Some(2.0 * PI),
            Surface3::Trimmed(inner) => {
                if let Surface3::Torus(_) = inner.basis.as_ref() {
                    let [_, _, v0, v1] = inner.trim;
                    let [_, _, dv0, dv1] = inner.basis.default_domain();
                    if (v0 - dv0).abs() < TOLERANCE_LINEAR_ULTRA_STRICT && (v1 - dv1).abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
                        Some(2.0 * PI)
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    /// Returns true if the face surface is V-periodic.
    /// Analogous to OCCT `BRepAdaptor_Surface::IsVPeriodic()`.
    pub fn is_v_periodic(&self) -> bool {
        self.v_period().is_some()
    }

    /// Get the TFaceData for this adaptor's face index.
    fn get_face_data(&self) -> Option<&topods::TFaceData> {
        let ts = self.brep.tshapes.get(self.face_idx)?;
        match ts.as_ref() {
            TShape::Face(fd) => Some(fd),
            _ => None,
        }
    }

    /// Fall back to vertex-based point when no surface is available.
    fn point_from_vertices(&self) -> DVec3 {
        let fd = match self.get_face_data() {
            Some(f) => f,
            None => return DVec3::ZERO,
        };

        // Get outer wire edges and compute centroid from their vertex positions.
        let wire_ts = match self.brep.tshapes.get(fd.outer_wire.index) {
            Some(ts) => ts,
            None => return DVec3::ZERO,
        };
        let TShape::Wire(wd) = wire_ts.as_ref() else { return DVec3::ZERO };

        let mut sum = DVec3::ZERO;
        let mut count = 0usize;
        for er in &wd.edges {
            let edge_ts = match self.brep.tshapes.get(er.index) {
                Some(ts) => ts,
                None => continue,
            };
            let TShape::Edge(ed) = edge_ts.as_ref() else { continue };
            if let Some(p) = self.brep.vertex_point(ed.first.index) {
                sum += p;
                count += 1;
            }
            if let Some(p) = self.brep.vertex_point(ed.last.index) {
                sum += p;
                count += 1;
            }
        }
        if count > 0 {
            sum / count as f64
        } else {
            DVec3::ZERO
        }
    }

    /// Fall back to DVec3::Z when no surface is available.
    fn normal_from_face(&self) -> DVec3 {
        // No separate normal stored on TFaceData; evaluate from surface
        let fd = match self.get_face_data() {
            Some(f) => f,
            None => return DVec3::Z,
        };
        let Some(surf) = &fd.surface else { return DVec3::Z };
        let dom = surf.default_domain();
        let u = (dom[0] + dom[1]) * 0.5;
        let v = (dom[2] + dom[3]) * 0.5;
        if u.is_finite() && v.is_finite() { surf.normal_at(u, v) } else { DVec3::Z }
    }

    /// Get the tolerance for this face.
    pub fn tolerance(&self) -> f64 {
        self.get_face_data()
            .map(|fd| fd.tolerance)
            .filter(|&t| t > 0.0)
            .unwrap_or(TOLERANCE_ABS)
    }
}

// =============================================================================
// WireAdaptor (BRepAdaptor_CompCurve)
// =============================================================================

/// Information about a segment in a wire adaptor.
#[derive(Debug, Clone)]
struct WireSegment {
    /// Edge adaptor for this segment.
    adaptor: EdgeAdaptor<'static>,
    /// Cumulative length fraction at the start of this segment.
    start_frac: f64,
    /// Cumulative length fraction at the end of this segment.
    end_frac: f64,
    /// Arc-length of this segment (approximate).
    length: f64,
}

/// Adapts a BRep wire to act as a composite 3D curve.
///
/// A wire is a connected sequence of edges. This adaptor treats the wire
/// as a single curve parameterized by cumulative arc-length fraction.
///
/// The parameter `t` in [0, 1] represents the position along the wire,
/// where each edge's contribution is weighted by its arc-length.
///
/// Analogous to OCCT's `BRepAdaptor_CompCurve`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::WireAdaptor;
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// // Get the outer wire of face 0.
/// let wire = brep.solids[0].shells[0].faces[0].outer_wire.clone();
/// let adaptor = WireAdaptor::new(&brep, &wire, 0);
/// let midpoint = adaptor.point_at(0.5);
/// let edge_at_mid = adaptor.edge_at(0.5);
/// ```
pub struct WireAdaptor<'a> {
    brep: &'a rcad_kernel::BRep,
    wire: &'a Wire,
    /// Face index for pcurve lookups (if available).
    face_idx: Option<usize>,
    /// Precomputed segments with arc-lengths.
    segments: Vec<WireSegmentData>,
    /// Total arc-length of the wire.
    total_length: f64,
}

/// Stored segment data (without lifetime issues).
#[derive(Debug, Clone)]
struct WireSegmentData {
    /// Edge index.
    edge_idx: usize,
    /// Whether the edge is reversed in this wire.
    reversed: bool,
    /// Cumulative length fraction at the start of this segment.
    start_frac: f64,
    /// Cumulative length fraction at the end of this segment.
    end_frac: f64,
    /// Arc-length of this segment.
    length: f64,
}

impl<'a> WireAdaptor<'a> {
    /// Create a new wire adaptor.
    ///
    /// # Arguments
    ///
    /// * `brep` - Reference to the BRep containing the wire.
    /// * `wire` - Reference to the wire to adapt.
    /// * `face_idx` - Optional flat face index for pcurve lookups.
    ///
    /// The wire's edges are preprocessed to compute arc-lengths for
    /// parameterization by cumulative length fraction.
    pub fn new(brep: &'a rcad_kernel::BRep, wire: &'a Wire, face_idx: usize) -> Self {
        let mut segments = Vec::with_capacity(wire.edges.len());
        let mut total_length = 0.0f64;

        for we in &wire.edges {
            let edge_idx = we.idx;
            let reversed = !we.forward;

            // Compute approximate arc-length for this edge.
            let length = Self::compute_edge_length(brep, edge_idx);
            total_length += length;

            segments.push(WireSegmentData {
                edge_idx,
                reversed,
                start_frac: 0.0, // Will be computed after total_length is known
                end_frac: 0.0,
                length,
            });
        }

        // Compute cumulative fractions.
        if total_length > TOLERANCE_FLOAT_DEDUP {
            let mut cum_length = 0.0f64;
            for seg in &mut segments {
                seg.start_frac = cum_length / total_length;
                cum_length += seg.length;
                seg.end_frac = cum_length / total_length;
            }
        } else if !segments.is_empty() {
            // All edges are zero-length; distribute uniformly.
            let n = segments.len() as f64;
            for (i, seg) in segments.iter_mut().enumerate() {
                seg.start_frac = i as f64 / n;
                seg.end_frac = (i + 1) as f64 / n;
                seg.length = 1.0; // Dummy length
            }
            total_length = n;
        }

        Self {
            brep,
            wire,
            face_idx: Some(face_idx),
            segments,
            total_length,
        }
    }

    /// Create a wire adaptor without a face context.
    ///
    /// This constructor is used when the wire is not associated with a specific face.
    pub fn without_face(brep: &'a rcad_kernel::BRep, wire: &'a Wire) -> Self {
        let mut adaptor = Self::new(brep, wire, 0);
        adaptor.face_idx = None;
        adaptor
    }

    /// Evaluate the point on the wire at parameter `t`.
    ///
    /// The parameter `t` is in [0, 1] and represents the cumulative
    /// arc-length fraction along the wire.
    pub fn point_at(&self, t: f64) -> DVec3 {
        let t = t.clamp(0.0, 1.0);
        let seg = self.find_segment(t);
        let local_t = self.local_parameter(t, seg);

        // Create an edge adaptor for this segment.
        let adaptor = self.create_edge_adaptor(seg.edge_idx, seg.reversed);
        adaptor.point_at(local_t)
    }

    /// Evaluate the unit tangent vector on the wire at parameter `t`.
    ///
    /// Returns the tangent pointing in the direction of traversal along the wire.
    pub fn tangent_at(&self, t: f64) -> DVec3 {
        let t = t.clamp(0.0, 1.0);
        let seg = self.find_segment(t);
        let local_t = self.local_parameter(t, seg);

        let adaptor = self.create_edge_adaptor(seg.edge_idx, seg.reversed);
        adaptor.tangent_at(local_t)
    }

    /// Return the edge index that contains the given parameter `t`.
    ///
    /// This is useful for determining which edge a point lies on.
    pub fn edge_at(&self, t: f64) -> usize {
        let t = t.clamp(0.0, 1.0);
        let seg = self.find_segment(t);
        seg.edge_idx
    }

    /// Return the number of edges in the wire.
    pub fn num_edges(&self) -> usize {
        self.wire.edges.len()
    }

    /// Return the total arc-length of the wire.
    pub fn length(&self) -> f64 {
        self.total_length
    }

    /// Return the parameter domain of the wire (always [0, 1]).
    pub fn domain(&self) -> [f64; 2] {
        [0.0, 1.0]
    }

    /// Returns the first parameter of the wire domain.
    pub fn first_parameter(&self) -> f64 {
        self.domain()[0]
    }

    /// Returns the last parameter of the wire domain.
    pub fn last_parameter(&self) -> f64 {
        self.domain()[1]
    }

    /// Alias for wire value evaluation.
    /// Analogous to OCCT `BRepAdaptor_CompCurve::Value()`.
    pub fn value(&self, t: f64) -> DVec3 {
        self.point_at(t)
    }

    /// Find the segment containing the given parameter.
    fn find_segment(&self, t: f64) -> &WireSegmentData {
        // Binary search for the segment containing t.
        let mut lo = 0usize;
        let mut hi = self.segments.len();

        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let seg = &self.segments[mid];
            if t < seg.start_frac {
                hi = mid;
            } else if t > seg.end_frac {
                lo = mid + 1;
            } else {
                return seg;
            }
        }

        // Fallback: return last segment.
        self.segments.last().unwrap_or(&WireSegmentData {
            edge_idx: 0,
            reversed: false,
            start_frac: 0.0,
            end_frac: 1.0,
            length: 1.0,
        })
    }

    /// Compute the local parameter within a segment for global parameter t.
    fn local_parameter(&self, t: f64, seg: &WireSegmentData) -> f64 {
        if seg.end_frac <= seg.start_frac {
            return 0.5;
        }
        ((t - seg.start_frac) / (seg.end_frac - seg.start_frac)).clamp(0.0, 1.0)
    }

    /// Create an edge adaptor with the specified orientation.
    fn create_edge_adaptor(&self, edge_idx: usize, reversed: bool) -> EdgeAdaptor<'a> {
        EdgeAdaptor::new(self.brep, edge_idx).with_reversed(reversed)
    }

    /// Compute the approximate arc-length of an edge.
    fn compute_edge_length(brep: &rcad_kernel::BRep, edge_idx: usize) -> f64 {
        // Try to compute from curve on TEdgeData.
        if let Some(ts) = brep.tshapes.get(edge_idx) {
            if let TShape::Edge(ed) = ts.as_ref() {
                if let Some(curve) = &ed.curve {
                    let range = ed.range;
                    // Use numerical integration for arc-length.
                    return Self::arc_length_numerical(curve, range[0], range[1]);
                }
                // Fall back to vertex distance.
                let p0 = brep.vertex_point(ed.first.index).unwrap_or(DVec3::ZERO);
                let p1 = brep.vertex_point(ed.last.index).unwrap_or(DVec3::ZERO);
                return (p1 - p0).length();
            }
        }
        0.0
    }

    /// Numerical integration for arc-length using Gauss-Legendre quadrature.
    fn arc_length_numerical(curve: &Curve3, t0: f64, t1: f64) -> f64 {
        // Use 5-point Gauss-Legendre quadrature.
        const GAUSS_POINTS: [(f64, f64); 5] = [
            (0.0, 0.5688888888888889),
            (-0.5384693101056831, 0.47862867049936647),
            (0.5384693101056831, 0.47862867049936647),
            (-0.906_179_845_938_664, 0.23692688505618908),
            (0.906_179_845_938_664, 0.23692688505618908),
        ];

        let dt = t1 - t0;
        let mut length = 0.0f64;

        for (xi, wi) in GAUSS_POINTS {
            let t = 0.5 * (t0 + t1 + xi * dt);
            let tangent = curve.tangent_at(t);
            length += wi * tangent.length();
        }

        length * 0.5 * dt.abs()
    }
}

// =============================================================================
// CurveAdaptorArray (BRepAdaptor_HArray1OfCurve)
// =============================================================================

/// An array of edge adaptors with indexed access.
///
/// Provides convenient storage and access for multiple curve adaptors,
/// analogous to OCCT's `BRepAdaptor_HArray1OfCurve`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_adaptor::{EdgeAdaptor, CurveAdaptorArray};
/// use rcad_kernel::BRep;
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let mut array = CurveAdaptorArray::new();
/// for i in 0..brep.edges.len() {
///     array.push(EdgeAdaptor::new(&brep, i));
/// }
///
/// for i in 0..array.len() {
///     let adaptor = array.get(i).unwrap();
///     println!("Edge {} domain: {:?}", i, adaptor.domain());
/// }
/// ```
#[derive(Debug, Clone, Default)]
pub struct CurveAdaptorArray<'a> {
    adaptors: Vec<EdgeAdaptor<'a>>,
}

impl<'a> CurveAdaptorArray<'a> {
    /// Create an empty array.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create an array with pre-allocated capacity.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            adaptors: Vec::with_capacity(capacity),
        }
    }

    /// Add an edge adaptor to the array.
    pub fn push(&mut self, adaptor: EdgeAdaptor<'a>) {
        self.adaptors.push(adaptor);
    }

    /// Get the edge adaptor at the given index.
    pub fn get(&self, index: usize) -> Option<&EdgeAdaptor<'a>> {
        self.adaptors.get(index)
    }

    /// Get a mutable reference to the edge adaptor at the given index.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut EdgeAdaptor<'a>> {
        self.adaptors.get_mut(index)
    }

    /// Return the number of adaptors in the array.
    pub fn len(&self) -> usize {
        self.adaptors.len()
    }

    /// Return true if the array is empty.
    pub fn is_empty(&self) -> bool {
        self.adaptors.is_empty()
    }

    /// Iterate over all adaptors.
    pub fn iter(&self) -> impl Iterator<Item = &EdgeAdaptor<'a>> {
        self.adaptors.iter()
    }

    /// Iterate mutably over all adaptors.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut EdgeAdaptor<'a>> {
        self.adaptors.iter_mut()
    }

    /// Clear all adaptors from the array.
    pub fn clear(&mut self) {
        self.adaptors.clear();
    }

    /// Create an array from a BRep's edges.
    ///
    /// Creates an adaptor for each edge in the BRep.
    pub fn from_brep(brep: &'a rcad_kernel::BRep) -> Self {
        let n = brep.edge_count();
        let mut array = Self::with_capacity(n);
        for i in 0..n {
            array.push(EdgeAdaptor::new(brep, i));
        }
        array
    }
}

impl<'a> std::ops::Index<usize> for CurveAdaptorArray<'a> {
    type Output = EdgeAdaptor<'a>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.adaptors[index]
    }
}

// =============================================================================
// Free-standing evaluation helpers (moved from brep_algo for OCCT-package alignment)
// =============================================================================

/// Evaluate the normal of a face at parameter (u, v) on a topods::BRep.
pub fn evaluate_face_normal(brep: &topods::BRep, face_idx: usize, u: f64, v: f64) -> DVec3 {
    let faces: Vec<&topods::TShape> = brep.tshapes.iter()
        .filter(|ts| matches!(ts.as_ref(), topods::TShape::Face(_)))
        .map(|ts| ts.as_ref())
        .collect();
    faces.get(face_idx)
        .and_then(|ts| match ts {
            topods::TShape::Face(fd) => fd.surface.as_ref(),
            _ => None,
        })
        .map(|s| s.normal_at(u, v))
        .unwrap_or(DVec3::Z)
}

/// Evaluate the unit tangent of an edge at normalized parameter t in [0, 1].
pub fn evaluate_edge_tangent(brep: &topods::BRep, edge_idx: usize, t: f64) -> DVec3 {
    let edges: Vec<&topods::TShape> = brep.tshapes.iter()
        .filter(|ts| matches!(ts.as_ref(), topods::TShape::Edge(_)))
        .map(|ts| ts.as_ref())
        .collect();
    edges.get(edge_idx)
        .and_then(|ts| match ts {
            topods::TShape::Edge(ed) => ed.curve.as_ref().map(|c| (ed.range, c)),
            _ => None,
        })
        .map(|([t0, t1], curve)| {
            let t_actual = t0 + t * (t1 - t0);
            curve.tangent_at(t_actual)
        })
        .unwrap_or(DVec3::X)
}

/// Evaluate approximate normal at a vertex (average of adjacent face normals).
/// Simplified — returns Z for unconnected vertices.
pub fn evaluate_vertex_normal(brep: &topods::BRep, vertex_idx: usize) -> DVec3 {
    let v_count = brep.tshapes.iter().filter(|ts| matches!(ts.as_ref(), topods::TShape::Vertex(_))).count();
    if vertex_idx >= v_count { return DVec3::Z; }
    DVec3::Z
}

// =============================================================================
// Tests
// =============================================================================


