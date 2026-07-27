//! TopLoc-style location utilities for coordinate system management.
//!
//! This module provides utilities analogous to OCCT's `TopLoc` package:
//!
//! - **Location**: Coordinate transformations with composition and inversion
//! - **Datum**: Reference coordinate systems (local frames)
//! - **LocationManager**: Efficient storage and composition of multiple locations
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::top_loc::{Location, Datum, LocationManager};
//! use glam::{DVec3, DAffine3};
//!
//! // Create a translation
//! let loc = Location::from_translation(DVec3::new(5.0, 0.0, 0.0));
//! assert!(!loc.is_identity());
//!
//! // Compose locations
//! let rotation = Location::from_rotation(DVec3::Z, std::f64::consts::FRAC_PI_2);
//! let combined = loc.compose(&rotation);
//!
//! // Create a datum (local coordinate system)
//! let datum = Datum::from_origin_and_normal(DVec3::new(1.0, 2.0, 3.0), DVec3::Z);
//! let local_point = DVec3::new(0.5, 0.5, 0.0);
//! let world_point = datum.to_world(local_point);
//! ```

use crate::tolerance::*;
use glam::{DAffine3, DMat3, DMat4, DQuat, DVec3};
use rcad_kernel::topods;
use std::collections::HashMap;

// =============================================================================
// Location - Coordinate Transformation
// =============================================================================

/// A coordinate transformation with efficient composition and inversion.
///
/// `Location` wraps an affine transformation and provides methods for
/// creating common transformations, composition, and inversion.
///
/// This is analogous to OCCT's `TopLoc_Location` class.
#[derive(Debug, Clone, PartialEq)]
pub struct Location {
    /// The underlying affine transformation.
    transform: DAffine3,
    /// Cached flag for identity check.
    is_identity_cache: bool,
}

impl Location {
    /// Create the identity location (no transformation).
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    ///
    /// let identity = Location::identity();
    /// assert!(identity.is_identity());
    /// ```
    pub fn identity() -> Self {
        Location {
            transform: DAffine3::IDENTITY,
            is_identity_cache: true,
        }
    }

    /// Create a location from an affine transformation.
    ///
    /// # Arguments
    ///
    /// * `transform` - The affine transformation matrix.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::{DAffine3, DVec3};
    ///
    /// let affine = DAffine3::from_translation(DVec3::new(1.0, 2.0, 3.0));
    /// let loc = Location::from_transform(affine);
    /// assert_eq!(*loc.transform(), affine);
    /// ```
    pub fn from_transform(transform: DAffine3) -> Self {
        let is_identity = transform == DAffine3::IDENTITY;
        Location {
            transform,
            is_identity_cache: is_identity,
        }
    }

    /// Create a location from a translation vector.
    ///
    /// # Arguments
    ///
    /// * `translation` - The translation vector.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::DVec3;
    ///
    /// let loc = Location::from_translation(DVec3::new(5.0, 0.0, 0.0));
    /// assert!(!loc.is_identity());
    /// ```
    pub fn from_translation(translation: DVec3) -> Self {
        Location {
            transform: DAffine3::from_translation(translation),
            is_identity_cache: translation == DVec3::ZERO,
        }
    }

    /// Create a location from a rotation around an axis.
    ///
    /// # Arguments
    ///
    /// * `axis` - The axis of rotation (will be normalized).
    /// * `angle` - The rotation angle in radians.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::DVec3;
    /// use std::f64::consts::FRAC_PI_2;
    ///
    /// // Rotate 90 degrees around Z axis
    /// let loc = Location::from_rotation(DVec3::Z, FRAC_PI_2);
    /// assert!(!loc.is_identity());
    /// ```
    pub fn from_rotation(axis: DVec3, angle: f64) -> Self {
        if angle.abs() < TOLERANCE_LEN_MIN {
            return Self::identity();
        }
        let normalized_axis = axis.normalize_or(DVec3::Z);
        Location {
            transform: DAffine3::from_axis_angle(normalized_axis, angle),
            is_identity_cache: false,
        }
    }

    /// Create a location from a uniform scale factor.
    ///
    /// # Arguments
    ///
    /// * `scale` - The uniform scale factor.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    ///
    /// let loc = Location::from_scale(2.0);
    /// assert!(!loc.is_identity());
    /// ```
    pub fn from_scale(scale: f64) -> Self {
        if (scale - 1.0).abs() < TOLERANCE_LEN_MIN {
            return Self::identity();
        }
        Location {
            transform: DAffine3::from_scale(DVec3::splat(scale)),
            is_identity_cache: false,
        }
    }

    /// Create a location from a quaternion rotation.
    ///
    /// # Arguments
    ///
    /// * `quat` - The quaternion representing the rotation.
    pub fn from_quaternion(quat: DQuat) -> Self {
        Location {
            transform: DAffine3::from_quat(quat),
            is_identity_cache: quat == DQuat::IDENTITY,
        }
    }

    /// Create a location from translation, rotation (as quaternion), and scale.
    ///
    /// # Arguments
    ///
    /// * `translation` - The translation vector.
    /// * `rotation` - The rotation quaternion.
    /// * `scale` - The uniform scale factor.
    pub fn from_trs(translation: DVec3, rotation: DQuat, scale: f64) -> Self {
        let is_identity = translation == DVec3::ZERO
            && rotation == DQuat::IDENTITY
            && (scale - 1.0).abs() < TOLERANCE_LEN_MIN;

        let scale_vec = DVec3::splat(scale);
        Location {
            transform: DAffine3::from_scale_rotation_translation(scale_vec, rotation, translation),
            is_identity_cache: is_identity,
        }
    }

    /// Create a location from a look-at transformation.
    ///
    /// Creates a coordinate frame where:
    /// - The origin is at `eye`
    /// - The -Z axis points toward `target`
    /// - The Y axis is as close to `up` as possible
    ///
    /// # Arguments
    ///
    /// * `eye` - The position of the camera/observer.
    /// * `target` - The point to look at.
    /// * `up` - The approximate up direction.
    pub fn from_look_at(eye: DVec3, target: DVec3, up: DVec3) -> Self {
        Location {
            transform: DAffine3::look_at_rh(eye, target, up),
            is_identity_cache: false,
        }
    }

    /// Compose this location with another location.
    ///
    /// The resulting transformation applies `self` first, then `other`.
    /// This is equivalent to matrix multiplication: `other * self`.
    ///
    /// # Arguments
    ///
    /// * `other` - The location to compose with.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::DVec3;
    ///
    /// let translate = Location::from_translation(DVec3::new(1.0, 0.0, 0.0));
    /// let scale = Location::from_scale(2.0);
    ///
    /// let combined = translate.compose(&scale);
    /// // First translates by 1, then scales by 2
    /// ```
    pub fn compose(&self, other: &Location) -> Location {
        if self.is_identity_cache {
            return other.clone();
        }
        if other.is_identity_cache {
            return self.clone();
        }
        // The resulting transformation applies `self` first, then `other`.
        // This is equivalent to matrix multiplication: `other * self`.
        let transform = other.transform * self.transform;
        // Check if the result is identity (within tolerance)
        let is_identity = transform
            .to_cols_array()
            .iter()
            .zip(DAffine3::IDENTITY.to_cols_array().iter())
            .all(|(a, b)| (a - b).abs() < TOLERANCE_LEN_MIN);
        Location {
            transform,
            is_identity_cache: is_identity,
        }
    }

    /// Compute the inverse of this location.
    ///
    /// The inverse location undoes the transformation of this location.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::DVec3;
    ///
    /// let loc = Location::from_translation(DVec3::new(5.0, 0.0, 0.0));
    /// let inv = loc.inverse();
    ///
    /// let combined = loc.compose(&inv);
    /// assert!(combined.is_identity());
    /// ```
    pub fn inverse(&self) -> Location {
        if self.is_identity_cache {
            return Self::identity();
        }
        Location {
            transform: self.transform.inverse(),
            is_identity_cache: false,
        }
    }

    /// Get a reference to the underlying affine transformation.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::{DAffine3, DVec3};
    ///
    /// let loc = Location::from_translation(DVec3::new(1.0, 2.0, 3.0));
    /// let transform = loc.transform();
    /// ```
    pub fn transform(&self) -> &DAffine3 {
        &self.transform
    }

    /// Check if this is the identity transformation.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    ///
    /// let identity = Location::identity();
    /// assert!(identity.is_identity());
    ///
    /// let translate = Location::from_translation(glam::DVec3::new(1.0, 0.0, 0.0));
    /// assert!(!translate.is_identity());
    /// ```
    pub fn is_identity(&self) -> bool {
        self.is_identity_cache
    }

    /// Transform a point by this location.
    ///
    /// # Arguments
    ///
    /// * `point` - The point to transform.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::DVec3;
    ///
    /// let loc = Location::from_translation(DVec3::new(5.0, 0.0, 0.0));
    /// let point = DVec3::new(1.0, 2.0, 3.0);
    /// let transformed = loc.transform_point(point);
    /// assert_eq!(transformed, DVec3::new(6.0, 2.0, 3.0));
    /// ```
    pub fn transform_point(&self, point: DVec3) -> DVec3 {
        self.transform.transform_point3(point)
    }

    /// Transform a vector by this location (ignoring translation).
    ///
    /// # Arguments
    ///
    /// * `vector` - The vector to transform.
    pub fn transform_vector(&self, vector: DVec3) -> DVec3 {
        self.transform.transform_vector3(vector)
    }

    /// Extract the translation component.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Location;
    /// use glam::DVec3;
    ///
    /// let loc = Location::from_translation(DVec3::new(1.0, 2.0, 3.0));
    /// assert_eq!(loc.translation(), DVec3::new(1.0, 2.0, 3.0));
    /// ```
    pub fn translation(&self) -> DVec3 {
        self.transform.translation
    }

    /// Extract the rotation as a quaternion.
    ///
    /// Note: This assumes no scale component.
    pub fn rotation_quat(&self) -> DQuat {
        DQuat::from_mat3(&self.transform.matrix3)
    }

    /// Extract the scale component (assuming uniform scale).
    ///
    /// Returns the scale factor, or 1.0 if the matrix is not uniformly scaled.
    pub fn scale(&self) -> f64 {
        // Extract scale from the matrix columns
        let sx = self.transform.matrix3.x_axis.length();
        let sy = self.transform.matrix3.y_axis.length();
        let sz = self.transform.matrix3.z_axis.length();
        // Return average scale if approximately uniform, otherwise 1.0
        let avg = (sx + sy + sz) / 3.0;
        if (sx - avg).abs() < TOLERANCE_COORD_SUB
            && (sy - avg).abs() < TOLERANCE_COORD_SUB
            && (sz - avg).abs() < TOLERANCE_COORD_SUB
        {
            avg
        } else {
            1.0
        }
    }

    /// Create a location from a rotation around the X axis.
    pub fn from_rotation_x(angle: f64) -> Self {
        Self::from_rotation(DVec3::X, angle)
    }

    /// Create a location from a rotation around the Y axis.
    pub fn from_rotation_y(angle: f64) -> Self {
        Self::from_rotation(DVec3::Y, angle)
    }

    /// Create a location from a rotation around the Z axis.
    pub fn from_rotation_z(angle: f64) -> Self {
        Self::from_rotation(DVec3::Z, angle)
    }

    /// Create a location from Euler angles (ZYX order).
    ///
    /// # Arguments
    ///
    /// * `roll` - Rotation around X axis.
    /// * `pitch` - Rotation around Y axis.
    /// * `yaw` - Rotation around Z axis.
    pub fn from_euler_angles(roll: f64, pitch: f64, yaw: f64) -> Self {
        let quat = DQuat::from_euler(glam::EulerRot::ZYX, yaw, pitch, roll);
        Self::from_quaternion(quat)
    }

    /// Create a location from a 4x4 matrix.
    ///
    /// # Arguments
    ///
    /// * `matrix` - The 4x4 transformation matrix.
    ///
    /// # Note
    ///
    /// This extracts the 3x4 affine part of the matrix. The last row is ignored.
    pub fn from_matrix(matrix: DMat4) -> Self {
        let affine = DAffine3::from_cols(
            matrix.x_axis.truncate(),
            matrix.y_axis.truncate(),
            matrix.z_axis.truncate(),
            matrix.w_axis.truncate(),
        );
        Self::from_transform(affine)
    }

    /// Convert to a 4x4 matrix.
    ///
    /// The last row is set to [0, 0, 0, 1].
    pub fn to_matrix(&self) -> DMat4 {
        DMat4::from_cols(
            self.transform.matrix3.x_axis.extend(0.0),
            self.transform.matrix3.y_axis.extend(0.0),
            self.transform.matrix3.z_axis.extend(0.0),
            self.transform.translation.extend(1.0),
        )
    }
}

impl Default for Location {
    fn default() -> Self {
        Self::identity()
    }
}

impl std::ops::Mul for &Location {
    type Output = Location;

    fn mul(self, rhs: &Location) -> Location {
        self.compose(rhs)
    }
}

impl std::ops::Mul for Location {
    type Output = Location;

    fn mul(self, rhs: Location) -> Location {
        self.compose(&rhs)
    }
}

// =============================================================================
// Datum - Reference Coordinate System
// =============================================================================

/// A reference coordinate system (local frame).
///
/// A `Datum` represents a local coordinate system defined by an origin
/// and three orthogonal axes (X, Y, Z).
///
/// This is analogous to OCCT's `TopLoc_Datum` concept.
#[derive(Debug, Clone, PartialEq)]
pub struct Datum {
    /// The origin of the coordinate system.
    origin: DVec3,
    /// The X direction (normalized).
    x_dir: DVec3,
    /// The Y direction (normalized).
    y_dir: DVec3,
    /// The Z direction (normalized).
    z_dir: DVec3,
}

impl Datum {
    /// Create a new datum from origin and axis directions.
    ///
    /// The axes will be normalized and made orthogonal if necessary.
    ///
    /// # Arguments
    ///
    /// * `origin` - The origin of the coordinate system.
    /// * `x_dir` - The X direction (will be normalized).
    /// * `y_dir` - The Y direction (will be normalized).
    /// * `z_dir` - The Z direction (will be normalized).
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Datum;
    /// use glam::DVec3;
    ///
    /// let datum = Datum::new(
    /// DVec3::new(1.0, 2.0, 3.0),
    /// DVec3::X,
    /// DVec3::Y,
    /// DVec3::Z,
    /// );
    /// assert_eq!(datum.origin(), DVec3::new(1.0, 2.0, 3.0));
    /// ```
    pub fn new(origin: DVec3, x_dir: DVec3, y_dir: DVec3, z_dir: DVec3) -> Self {
        let x_dir = x_dir.normalize_or(DVec3::X);
        let y_dir = y_dir.normalize_or(DVec3::Y);
        let z_dir = z_dir.normalize_or(DVec3::Z);

        Datum {
            origin,
            x_dir,
            y_dir,
            z_dir,
        }
    }

    /// Create a datum from an origin and a normal direction.
    ///
    /// The Z axis will be aligned with the normal. X and Y will be
    /// automatically chosen to form an orthogonal right-handed frame.
    ///
    /// # Arguments
    ///
    /// * `origin` - The origin of the coordinate system.
    /// * `normal` - The Z direction (normal to the XY plane).
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Datum;
    /// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
    /// use glam::DVec3;
    ///
    /// let datum = Datum::from_origin_and_normal(
    /// DVec3::new(1.0, 2.0, 3.0),
    /// DVec3::Z,
    /// );
    /// // Z axis points along world Z
    /// # use rcad_algorithms::tolerance::*;
    /// assert!((datum.z_direction().dot(DVec3::Z) - 1.0).abs() < TOLERANCE_COORD_SUB);
    /// ```
    pub fn from_origin_and_normal(origin: DVec3, normal: DVec3) -> Self {
        let z_dir = normal.normalize_or(DVec3::Z);

        // Find a vector not parallel to z_dir
        let temp = if z_dir.dot(DVec3::X).abs() < 0.9 {
            DVec3::X
        } else {
            DVec3::Y
        };

        // For a right-handed coordinate system with z as normal:
        // y = z x, so x = temp projected onto the plane perpendicular to z
        // First get y as perpendicular to both z and temp
        // Then x = y z (completing the right-handed system)
        let y_dir = z_dir.cross(temp).normalize_or(DVec3::Y);
        let x_dir = y_dir.cross(z_dir).normalize_or(DVec3::X);

        Datum {
            origin,
            x_dir,
            y_dir,
            z_dir,
        }
    }

    /// Create a datum from an origin, normal, and X direction hint.
    ///
    /// The X axis will be as close as possible to the hint while remaining
    /// orthogonal to the normal.
    ///
    /// # Arguments
    ///
    /// * `origin` - The origin of the coordinate system.
    /// * `normal` - The Z direction (normal to the XY plane).
    /// * `x_hint` - Hint for the X direction (projected onto the plane).
    pub fn from_origin_normal_and_x(origin: DVec3, normal: DVec3, x_hint: DVec3) -> Self {
        let z_dir = normal.normalize_or(DVec3::Z);
        let x_dir = x_hint.normalize_or(DVec3::X);

        // Project x_hint onto the plane perpendicular to normal
        let x_proj = (x_dir - z_dir * z_dir.dot(x_dir)).normalize_or(DVec3::X);
        let y_dir = z_dir.cross(x_proj).normalize_or(DVec3::Y);

        Datum {
            origin,
            x_dir: x_proj,
            y_dir,
            z_dir,
        }
    }

    /// Create a datum from three points.
    ///
    /// The origin is at `origin`, X points toward `x_point`,
    /// and the XY plane contains `y_point`.
    ///
    /// # Arguments
    ///
    /// * `origin` - The origin of the coordinate system.
    /// * `x_point` - A point on the X axis (not at origin).
    /// * `y_point` - A point in the XY plane (not on X axis).
    pub fn from_three_points(origin: DVec3, x_point: DVec3, y_point: DVec3) -> Self {
        let x_dir = (x_point - origin).normalize_or(DVec3::X);

        // Vector from origin to y_point
        let oy = y_point - origin;
        // Project oy onto the plane perpendicular to x_dir
        let y_proj = oy - x_dir * x_dir.dot(oy);
        let y_dir = y_proj.normalize_or(DVec3::Y);

        let z_dir = x_dir.cross(y_dir).normalize_or(DVec3::Z);

        Datum {
            origin,
            x_dir,
            y_dir,
            z_dir,
        }
    }

    /// Create a datum representing the world coordinate system.
    pub fn world() -> Self {
        Datum {
            origin: DVec3::ZERO,
            x_dir: DVec3::X,
            y_dir: DVec3::Y,
            z_dir: DVec3::Z,
        }
    }

    /// Get the origin of this datum.
    pub fn origin(&self) -> DVec3 {
        self.origin
    }

    /// Get the X direction of this datum.
    pub fn x_direction(&self) -> DVec3 {
        self.x_dir
    }

    /// Get the Y direction of this datum.
    pub fn y_direction(&self) -> DVec3 {
        self.y_dir
    }

    /// Get the Z direction of this datum.
    pub fn z_direction(&self) -> DVec3 {
        self.z_dir
    }

    /// Transform a local point to world coordinates.
    ///
    /// # Arguments
    ///
    /// * `local` - A point in local coordinates.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Datum;
    /// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
    /// use glam::DVec3;
    ///
    /// let datum = Datum::from_origin_and_normal(DVec3::new(1.0, 0.0, 0.0), DVec3::Z);
    /// let local = DVec3::new(0.5, 0.0, 0.0);
    /// let world = datum.to_world(local);
    /// # use rcad_algorithms::tolerance::*;
    /// assert!((world.x - 1.5).abs() < TOLERANCE_COORD_SUB); // 1.0 (origin) + 0.5 (local x)
    /// ```
    pub fn to_world(&self, local: DVec3) -> DVec3 {
        self.origin + self.x_dir * local.x + self.y_dir * local.y + self.z_dir * local.z
    }

    /// Transform a world point to local coordinates.
    ///
    /// # Arguments
    ///
    /// * `world` - A point in world coordinates.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::Datum;
    /// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
    /// use glam::DVec3;
    ///
    /// let datum = Datum::from_origin_and_normal(DVec3::new(1.0, 0.0, 0.0), DVec3::Z);
    /// let world = DVec3::new(1.5, 0.0, 0.0);
    /// let local = datum.to_local(world);
    /// # use rcad_algorithms::tolerance::*;
    /// assert!((local.x - 0.5).abs() < TOLERANCE_COORD_SUB); // world x - origin x
    /// ```
    pub fn to_local(&self, world: DVec3) -> DVec3 {
        let diff = world - self.origin;
        DVec3::new(
            diff.dot(self.x_dir),
            diff.dot(self.y_dir),
            diff.dot(self.z_dir),
        )
    }

    /// Transform a local vector to world coordinates (no translation).
    ///
    /// # Arguments
    ///
    /// * `local` - A vector in local coordinates.
    pub fn vector_to_world(&self, local: DVec3) -> DVec3 {
        self.x_dir * local.x + self.y_dir * local.y + self.z_dir * local.z
    }

    /// Transform a world vector to local coordinates (no translation).
    ///
    /// # Arguments
    ///
    /// * `world` - A vector in world coordinates.
    pub fn vector_to_local(&self, world: DVec3) -> DVec3 {
        DVec3::new(
            world.dot(self.x_dir),
            world.dot(self.y_dir),
            world.dot(self.z_dir),
        )
    }

    /// Convert this datum to a location.
    ///
    /// The resulting location transforms points from the local frame
    /// to the world frame.
    pub fn to_location(&self) -> Location {
        let matrix = DMat3::from_cols(self.x_dir, self.y_dir, self.z_dir);
        Location::from_transform(DAffine3::from_mat3_translation(matrix, self.origin))
    }

    /// Create a datum from a location.
    ///
    /// The origin and axes are extracted from the transformation.
    pub fn from_location(loc: &Location) -> Self {
        let transform = loc.transform();
        let matrix = transform.matrix3;

        Datum {
            origin: transform.translation,
            x_dir: matrix.x_axis.normalize_or(DVec3::X),
            y_dir: matrix.y_axis.normalize_or(DVec3::Y),
            z_dir: matrix.z_axis.normalize_or(DVec3::Z),
        }
    }

    /// Get the axes as a 3x3 rotation matrix.
    pub fn axes_matrix(&self) -> DMat3 {
        DMat3::from_cols(self.x_dir, self.y_dir, self.z_dir)
    }

    /// Check if the datum represents the world coordinate system.
    pub fn is_world(&self) -> bool {
        self.origin == DVec3::ZERO
            && (self.x_dir - DVec3::X).length() < TOLERANCE_COORD_SUB
            && (self.y_dir - DVec3::Y).length() < TOLERANCE_COORD_SUB
            && (self.z_dir - DVec3::Z).length() < TOLERANCE_COORD_SUB
    }

    /// Create a new datum with a different origin.
    pub fn with_origin(&self, origin: DVec3) -> Self {
        Datum {
            origin,
            x_dir: self.x_dir,
            y_dir: self.y_dir,
            z_dir: self.z_dir,
        }
    }
}

impl Default for Datum {
    fn default() -> Self {
        Self::world()
    }
}

// =============================================================================
// LocationManager
// =============================================================================

/// Manages a collection of locations with efficient composition.
///
/// `LocationManager` provides storage and retrieval of multiple locations,
/// with efficient composition by index.
///
/// This is analogous to OCCT's location management in `TopLoc`.
#[derive(Debug, Clone, Default)]
pub struct LocationManager {
    /// Stored locations.
    locations: Vec<Location>,
    /// Cache for composed locations.
    compose_cache: HashMap<Vec<usize>, Location>,
}

impl LocationManager {
    /// Create a new empty location manager.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::LocationManager;
    ///
    /// let manager = LocationManager::new();
    /// assert_eq!(manager.len(), 0);
    /// ```
    pub fn new() -> Self {
        LocationManager {
            locations: Vec::new(),
            compose_cache: HashMap::new(),
        }
    }

    /// Create a location manager with pre-allocated capacity.
    ///
    /// # Arguments
    ///
    /// * `capacity` - Initial capacity for the location storage.
    pub fn with_capacity(capacity: usize) -> Self {
        LocationManager {
            locations: Vec::with_capacity(capacity),
            compose_cache: HashMap::new(),
        }
    }

    /// Add a location and return its index.
    ///
    /// # Arguments
    ///
    /// * `loc` - The location to add.
    ///
    /// # Returns
    ///
    /// The index of the added location.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::{LocationManager, Location};
    /// use glam::DVec3;
    ///
    /// let mut manager = LocationManager::new();
    /// let idx = manager.add_location(Location::from_translation(DVec3::new(1.0, 0.0, 0.0)));
    /// assert_eq!(idx, 0);
    /// ```
    pub fn add_location(&mut self, loc: Location) -> usize {
        let idx = self.locations.len();
        self.locations.push(loc);
        idx
    }

    /// Get a reference to a location by index.
    ///
    /// # Arguments
    ///
    /// * `idx` - The index of the location.
    ///
    /// # Returns
    ///
    /// `Some(&Location)` if the index is valid, `None` otherwise.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::{LocationManager, Location};
    /// use glam::DVec3;
    ///
    /// let mut manager = LocationManager::new();
    /// let idx = manager.add_location(Location::from_translation(DVec3::new(1.0, 0.0, 0.0)));
    /// let loc = manager.get_location(idx);
    /// assert!(loc.is_some());
    /// assert!(manager.get_location(999).is_none());
    /// ```
    pub fn get_location(&self, idx: usize) -> Option<&Location> {
        self.locations.get(idx)
    }

    /// Get a mutable reference to a location by index.
    ///
    /// # Arguments
    ///
    /// * `idx` - The index of the location.
    pub fn get_location_mut(&mut self, idx: usize) -> Option<&mut Location> {
        self.locations.get_mut(idx)
    }

    /// Compose multiple locations by their indices.
    ///
    /// The locations are composed in order: indices[0] * indices[1] * ...
    ///
    /// # Arguments
    ///
    /// * `indices` - The indices of the locations to compose.
    ///
    /// # Returns
    ///
    /// The composed location, or the identity if `indices` is empty.
    ///
    /// # Example
    ///
    /// ```
    /// use rcad_algorithms::top_loc::{LocationManager, Location};
    /// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
    /// use glam::DVec3;
    ///
    /// let mut manager = LocationManager::new();
    /// let t1 = manager.add_location(Location::from_translation(DVec3::new(1.0, 0.0, 0.0)));
    /// let t2 = manager.add_location(Location::from_translation(DVec3::new(0.0, 2.0, 0.0)));
    ///
    /// let composed = manager.compose_locations(&[t1, t2]);
    /// let result = composed.transform_point(DVec3::ZERO);
    /// assert!((result.x - 1.0).abs() < TOLERANCE_COORD_SUB);
    /// # use rcad_algorithms::tolerance::*;
    /// assert!((result.y - 2.0).abs() < TOLERANCE_COORD_SUB);
    /// ```
    pub fn compose_locations(&self, indices: &[usize]) -> Location {
        if indices.is_empty() {
            return Location::identity();
        }

        // Check cache
        if let Some(cached) = self.compose_cache.get(indices) {
            return cached.clone();
        }

        // Compose all locations
        let mut result = Location::identity();
        for &idx in indices {
            if let Some(loc) = self.locations.get(idx) {
                result = result.compose(loc);
            }
        }

        result
    }

    /// Compose locations and cache the result.
    ///
    /// This is useful when the same composition is used repeatedly.
    pub fn compose_locations_cached(&mut self, indices: &[usize]) -> Location {
        if indices.is_empty() {
            return Location::identity();
        }

        // Check cache first
        if let Some(cached) = self.compose_cache.get(indices) {
            return cached.clone();
        }

        // Compute and cache
        let result = self.compose_locations(indices);
        self.compose_cache.insert(indices.to_vec(), result.clone());
        result
    }

    /// Clear the composition cache.
    pub fn clear_cache(&mut self) {
        self.compose_cache.clear();
    }

    /// Get the number of stored locations.
    pub fn len(&self) -> usize {
        self.locations.len()
    }

    /// Check if the manager is empty.
    pub fn is_empty(&self) -> bool {
        self.locations.is_empty()
    }

    /// Get an iterator over all stored locations.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Location)> {
        self.locations.iter().enumerate()
    }

    /// Remove all locations and clear the cache.
    pub fn clear(&mut self) {
        self.locations.clear();
        self.compose_cache.clear();
    }

    /// Reserve capacity for additional locations.
    pub fn reserve(&mut self, additional: usize) {
        self.locations.reserve(additional);
    }

    /// Find locations that match the given predicate.
    ///
    /// # Arguments
    ///
    /// * `predicate` - A function that returns `true` for matching locations.
    ///
    /// # Returns
    ///
    /// A vector of indices for matching locations.
    pub fn find<F>(&self, predicate: F) -> Vec<usize>
    where
        F: Fn(&Location) -> bool,
    {
        self.locations
            .iter()
            .enumerate()
            .filter_map(|(i, loc)| if predicate(loc) { Some(i) } else { None })
            .collect()
    }

    /// Find the index of a specific location.
    ///
    /// # Arguments
    ///
    /// * `loc` - The location to find.
    ///
    /// # Returns
    ///
    /// `Some(index)` if found, `None` otherwise.
    pub fn find_location(&self, loc: &Location) -> Option<usize> {
        self.locations.iter().position(|l| l == loc)
    }

    /// Add a location only if it doesn't already exist.
    ///
    /// # Arguments
    ///
    /// * `loc` - The location to add.
    ///
    /// # Returns
    ///
    /// The index of the location (existing or newly added).
    pub fn add_location_dedup(&mut self, loc: Location) -> usize {
        if let Some(idx) = self.find_location(&loc) {
            return idx;
        }
        self.add_location(loc)
    }
}

// =============================================================================
// Shape Location Application
// =============================================================================

/// Apply a location transformation to a BRep shape.
///
/// This function transforms all geometry in the BRep according to the
/// given location.
///
/// # Arguments
///
/// * `brep` - The BRep to transform (modified in place).
/// * `loc` - The location (transformation) to apply.
///
/// # Example
///
/// ```
/// # use rcad_algorithms::tolerance::*;
/// use rcad_algorithms::top_loc::{Location, apply_location_to_shape};
/// use rcad_algorithms::tolerance::TOLERANCE_COORD_SUB;
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use glam::DVec3;
///
/// let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
/// width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let loc = Location::from_translation(DVec3::new(5.0, 0.0, 0.0));
/// apply_location_to_shape(&mut brep, &loc);
///
/// // The box is now centered at (5.5, 0.5, 0.5)
/// assert!((brep.vertex(0).point.x - 5.0).abs() < TOLERANCE_COORD_SUB);
/// ```
pub fn apply_location_to_shape(brep: &mut topods::BRep, loc: &Location) {
    if loc.is_identity() {
        return;
    }
    brep.apply_transform(*loc.transform());
}

/// Apply a location to a BRep and return a new transformed BRep.
///
/// # Arguments
///
/// * `brep` - The BRep to transform.
/// * `loc` - The location (transformation) to apply.
///
/// # Returns
///
/// A new BRep with the transformation applied.
pub fn apply_location_to_shape_owned(brep: &topods::BRep, loc: &Location) -> topods::BRep {
    if loc.is_identity() {
        return brep.clone();
    }
    let mut result = brep.clone();
    result.apply_transform(*loc.transform());
    result
}

// =============================================================================
// Tests
// =============================================================================
