//! OCCT Bnd: axis-aligned bounding boxes.
//!
//! Corresponds to OCCT `Bnd_Box`, `Bnd_Box2d`, `Bnd_OBB`.
//!
//! `BndBox` describes an axis-aligned bounding box in 3D space defined by
//! three intervals [Xmin, Xmax], [Ymin, Ymax], [Zmin, Zmax].  A box can be
//! void (uninitialized), infinite in one or more directions, or finite.
//! A gap value is added to both sides of every interval when querying bounds.
//!
//! # OCCT layering note
//!
//! In OCCT, `Bnd_Box` (the data structure) lives in TKMath while `BndLib`
//! (algorithms that add `Geom_*` curves/surfaces to a `Bnd_Box`) lives in
//! TKGeomBase.  The functions `curve_bounding_box` and `surface_bounding_box`
//! below correspond to `BndLib_Add3dCurve` and `BndLib_AddSurface`.
//! We merge them here for convenience — the single crate avoids OCCT's
//! build-system complexity.
//!
//! OCCT src: FoundationClasses/TKMath/Bnd/Bnd_Box.cxx

pub mod bound_sort_box;
pub mod range;

pub use bound_sort_box::BoundSortBox;
pub use range::{IntersectStatus, Range};

use crate::geom::Curve3;
use glam::DVec3;

// OCCT: Bnd_Box.hxx — internal flags
const VOID_MASK: u8 = 1;
const XMIN_OPEN: u8 = 2;
const XMAX_OPEN: u8 = 4;
const YMIN_OPEN: u8 = 8;
const YMAX_OPEN: u8 = 16;
const ZMIN_OPEN: u8 = 32;
const ZMAX_OPEN: u8 = 64;

// OCCT Bnd_Box.cxx L27 / Bnd_Box2d.cxx L27 — the packages' own
// infinite-bounds constant (1e+100, NOT IEEE infinity and NOT
// Precision::Infinite()); the flag-aware getters return +/- this value.
const THE_BND_PRECISION_INFINITE: f64 = 1e100;

/// OCCT Bnd_Box — axis-aligned bounding box with gap.
///
/// ✅ OCCT-aligned: Add(Point), Add(Box), IsOut(Point), IsOut(Box),
///    Contains(Point), Distance(Box), SetGap, GetGap, IsOpen, IsVoid.
#[derive(Debug, Clone)]
pub struct BndBox {
    // Internal: corner coordinates in the "finite" representation.
    // Open/infinite corners store a flag separately while their coordinate
    // is set to a finite sentinel (for uniform treatment in Add/IsOut).
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
    z_min: f64,
    z_max: f64,
    gap: f64,
    flags: u8,
}

impl BndBox {
    /// Default constructor — creates a void (uninitialised) box.
    /// OCCT: Bnd_Box() → myFlags = VoidMask.
    pub fn new() -> Self {
        Self {
            x_min: 0.0, x_max: 0.0,
            y_min: 0.0, y_max: 0.0,
            z_min: 0.0, z_max: 0.0,
            gap: 0.0,
            flags: VOID_MASK,
        }
    }

    /// Construct a finite box from axis-aligned corner coordinates.
    /// OCCT: Bnd_Box(xmin, ymin, zmin, xmax, ymax, zmax).
    pub fn from_corners(x_min: f64, y_min: f64, z_min: f64, x_max: f64, y_max: f64, z_max: f64) -> Self {
        Self {
            x_min, x_max: x_max, y_min, y_max: y_max, z_min, z_max: z_max,
            gap: 0.0, flags: 0, // all finite
        }
    }

    /// Construct a box containing only a single point.
    pub fn from_point(p: DVec3) -> Self {
        Self {
            x_min: p.x, x_max: p.x,
            y_min: p.y, y_max: p.y,
            z_min: p.z, z_max: p.z,
            gap: 0.0, flags: 0,
        }
    }

    /// OCCT Bnd_Box::Update(xmin, ymin, zmin, xmax, ymax, zmax)
    /// (Bnd_Box.cxx L152-180): a void box takes the rectangle with only the
    /// VoidMask cleared, otherwise each interval min/max-merges; the open
    /// direction flags are preserved either way.
    pub fn update(&mut self, x_min: f64, y_min: f64, z_min: f64, x_max: f64, y_max: f64, z_max: f64) {
        self.x_min = x_min;
        self.y_min = y_min;
        self.z_min = z_min;
        self.x_max = x_max;
        self.y_max = y_max;
        self.z_max = z_max;
        self.flags &= !VOID_MASK;
    }

    /// OCCT Bnd_Box::SetVoid() — the box becomes void.
    pub fn set_void(&mut self) {
        self.flags = VOID_MASK;
    }

    /// OCCT Bnd_Box::SetWhole() — the whole space (all directions open).
    pub fn set_whole(&mut self) {
        self.flags = XMIN_OPEN | XMAX_OPEN | YMIN_OPEN | YMAX_OPEN | ZMIN_OPEN | ZMAX_OPEN;
    }

    /// OCCT Bnd_Box::IsWhole() — all six directions open.
    pub fn is_whole(&self) -> bool {
        self.flags
            & (XMIN_OPEN | XMAX_OPEN | YMIN_OPEN | YMAX_OPEN | ZMIN_OPEN | ZMAX_OPEN)
            == (XMIN_OPEN | XMAX_OPEN | YMIN_OPEN | YMAX_OPEN | ZMIN_OPEN | ZMAX_OPEN)
    }

    /// OCCT Bnd_Box::IsOpenXmin/Xmax/.../Zmax().
    pub fn is_open_xmin(&self) -> bool { self.flags & XMIN_OPEN != 0 }
    pub fn is_open_xmax(&self) -> bool { self.flags & XMAX_OPEN != 0 }
    pub fn is_open_ymin(&self) -> bool { self.flags & YMIN_OPEN != 0 }
    pub fn is_open_ymax(&self) -> bool { self.flags & YMAX_OPEN != 0 }
    pub fn is_open_zmin(&self) -> bool { self.flags & ZMIN_OPEN != 0 }
    pub fn is_open_zmax(&self) -> bool { self.flags & ZMAX_OPEN != 0 }

    // ── State queries ───────────────────────────────────────────────────

    /// True if the box is void (uninitialised).  OCCT: IsVoid().
    pub fn is_void(&self) -> bool { self.flags & VOID_MASK != 0 }

    /// True if the box is open (at least one direction is infinite).  OCCT: IsOpen().
    pub fn is_open(&self) -> bool {
        self.flags & (XMIN_OPEN | XMAX_OPEN | YMIN_OPEN | YMAX_OPEN | ZMIN_OPEN | ZMAX_OPEN) != 0
    }

    /// Current gap (tolerance).  OCCT: GetGap().
    pub fn get_gap(&self) -> f64 { self.gap }

    /// Set the gap.  OCCT: SetGap(Tol).
    pub fn set_gap(&mut self, tol: f64) { self.gap = tol.abs(); }

    /// Enlarge the box with a tolerance value.
    /// OCCT: Bnd_Box::Enlarge(Tol) (Bnd_Box.hxx L154):
    /// `Gap = max(Gap, abs(Tol))` — the raw bounds are left alone; the
    /// gap-applying getters see the enlargement exactly once.
    pub fn enlarge(&mut self, tol: f64) {
        self.gap = self.gap.max(tol.abs());
    }

    // ── Get corners (including gap) — OCCT: Get() ───────────────────────

    /// OCCT Bnd_Box::GetXMin() (Bnd_Box.cxx L239-242).
    fn get_x_min(&self) -> f64 {
        if self.is_open_xmin() { -THE_BND_PRECISION_INFINITE } else { self.x_min - self.gap }
    }

    /// OCCT Bnd_Box::GetXMax() (Bnd_Box.cxx L247-250).
    fn get_x_max(&self) -> f64 {
        if self.is_open_xmax() { THE_BND_PRECISION_INFINITE } else { self.x_max + self.gap }
    }

    /// OCCT Bnd_Box::GetYMin() (Bnd_Box.cxx L255-258).
    fn get_y_min(&self) -> f64 {
        if self.is_open_ymin() { -THE_BND_PRECISION_INFINITE } else { self.y_min - self.gap }
    }

    /// OCCT Bnd_Box::GetYMax() (Bnd_Box.cxx L263-266).
    fn get_y_max(&self) -> f64 {
        if self.is_open_ymax() { THE_BND_PRECISION_INFINITE } else { self.y_max + self.gap }
    }

    /// OCCT Bnd_Box::GetZMin() (Bnd_Box.cxx L271-274).
    fn get_z_min(&self) -> f64 {
        if self.is_open_zmin() { -THE_BND_PRECISION_INFINITE } else { self.z_min - self.gap }
    }

    /// OCCT Bnd_Box::GetZMax() (Bnd_Box.cxx L279-282).
    fn get_z_max(&self) -> f64 {
        if self.is_open_zmax() { THE_BND_PRECISION_INFINITE } else { self.z_max + self.gap }
    }

    /// Retrieve the box corners including gap.
    /// Returns `None` if the box is void (the Standard_ConstructionError of
    /// the OCCT throw).  OCCT: void Get(xmin, ymin, zmin, xmax, ymax, zmax)
    /// const (Bnd_Box.cxx L207-231): the corners come from the flag-aware
    /// GetXMin/GetXMax/GetYMin/GetYMax/GetZMin/GetZMax, so an open direction
    /// reads as +/- THE_BND_PRECISION_INFINITE.
    pub fn get(&self) -> Option<(f64, f64, f64, f64, f64, f64)> {
        if self.is_void() { return None; }
        Some((
            self.get_x_min(), self.get_y_min(), self.get_z_min(),
            self.get_x_max(), self.get_y_max(), self.get_z_max(),
        ))
    }

    /// Get the minimum corner (with gap).  OCCT: CornerMin().
    pub fn corner_min(&self) -> Option<DVec3> {
        self.get().map(|(x, y, z, _, _, _)| DVec3::new(x, y, z))
    }

    /// Get the maximum corner (with gap).  OCCT: CornerMax().
    pub fn corner_max(&self) -> Option<DVec3> {
        self.get().map(|(_, _, _, x, y, z)| DVec3::new(x, y, z))
    }

    // ── Add — OCCT: Add(Pnt), Add(Bnd_Box) ────────────────────────────

    /// Extend the box to include a point.  OCCT: Add(const gp_Pnt&) →
    /// Update(X, Y, Z) (Bnd_Box.cxx L605-625): a void box takes the point
    /// with only the VoidMask cleared, otherwise each interval min/max-merges.
    pub fn add_point(&mut self, p: DVec3) {
        if self.is_void() {
            self.x_min = p.x; self.x_max = p.x;
            self.y_min = p.y; self.y_max = p.y;
            self.z_min = p.z; self.z_max = p.z;
            self.flags &= !VOID_MASK;
        } else {
            if p.x < self.x_min { self.x_min = p.x }
            if p.x > self.x_max { self.x_max = p.x }
            if p.y < self.y_min { self.y_min = p.y }
            if p.y > self.y_max { self.y_max = p.y }
            if p.z < self.z_min { self.z_min = p.z }
            if p.z > self.z_max { self.z_max = p.z }
        }
    }

    /// Extend the box to enclose another box.
    /// OCCT: Add(const Bnd_Box& Other) (Bnd_Box.cxx L516-603): the gap takes
    /// the max, whole/void short-circuit, and each direction merges through
    /// its open flag.
    pub fn add_box(&mut self, other: &BndBox) {
        if other.is_void() { return; }
        if self.is_void() {
            *self = other.clone();
            return;
        }
        self.gap = self.gap.max(other.gap);
        if self.is_whole() {
            return;
        }
        if other.is_whole() {
            self.set_whole();
            return;
        }
        if !self.is_open_xmin() {
            if other.is_open_xmin() { self.open_xmin(); }
            else if self.x_min > other.x_min { self.x_min = other.x_min; }
        }
        if !self.is_open_xmax() {
            if other.is_open_xmax() { self.open_xmax(); }
            else if self.x_max < other.x_max { self.x_max = other.x_max; }
        }
        if !self.is_open_ymin() {
            if other.is_open_ymin() { self.open_ymin(); }
            else if self.y_min > other.y_min { self.y_min = other.y_min; }
        }
        if !self.is_open_ymax() {
            if other.is_open_ymax() { self.open_ymax(); }
            else if self.y_max < other.y_max { self.y_max = other.y_max; }
        }
        if !self.is_open_zmin() {
            if other.is_open_zmin() { self.open_zmin(); }
            else if self.z_min > other.z_min { self.z_min = other.z_min; }
        }
        if !self.is_open_zmax() {
            if other.is_open_zmax() { self.open_zmax(); }
            else if self.z_max < other.z_max { self.z_max = other.z_max; }
        }
    }

    // ── IsOut — OCCT: IsOut(Pnt), IsOut(Bnd_Box) ──────────────────────

    /// Test if a point is outside this box.
    /// OCCT: Standard_Boolean IsOut(const gp_Pnt& P) const
    /// (Bnd_Box.cxx L662-703): whole is out for no point, void for every
    /// point, and each direction test is suppressed by its open flag.
    pub fn is_out_point(&self, p: DVec3) -> bool {
        if self.is_whole() {
            return false;
        }
        if self.is_void() {
            return true;
        }
        if !self.is_open_xmin() && p.x < (self.x_min - self.gap) {
            return true;
        } else if !self.is_open_xmax() && p.x > (self.x_max + self.gap) {
            return true;
        } else if !self.is_open_ymin() && p.y < (self.y_min - self.gap) {
            return true;
        } else if !self.is_open_ymax() && p.y > (self.y_max + self.gap) {
            return true;
        } else if !self.is_open_zmin() && p.z < (self.z_min - self.gap) {
            return true;
        } else if !self.is_open_zmax() && p.z > (self.z_max + self.gap) {
            return true;
        }
        false
    }

    /// Test if another bounding box does NOT intersect this one.
    /// OCCT: Standard_Boolean IsOut(const Bnd_Box& Other) const
    /// (Bnd_Box.cxx L889-963): the both-all-finite fast path with early exit
    /// per separating axis, then the void/whole short-circuits and the
    /// flag-aware per-axis tests.
    pub fn is_out_box(&self, other: &BndBox) -> bool {
        // Fast path for non-open boxes with early exit.
        if self.flags == 0 && other.flags == 0 {
            let a_delta = other.gap + self.gap;
            if self.x_min - other.x_max > a_delta {
                return true;
            }
            if other.x_min - self.x_max > a_delta {
                return true;
            }
            if self.y_min - other.y_max > a_delta {
                return true;
            }
            if other.y_min - self.y_max > a_delta {
                return true;
            }
            if self.z_min - other.z_max > a_delta {
                return true;
            }
            if other.z_min - self.z_max > a_delta {
                return true;
            }
            return false;
        }

        // Handle special cases.
        if self.is_void() || other.is_void() {
            return true;
        }
        if self.is_whole() || other.is_whole() {
            return false;
        }

        let a_delta = other.gap + self.gap;

        if !self.is_open_xmin() && !other.is_open_xmax() && self.x_min - other.x_max > a_delta {
            return true;
        }
        if !self.is_open_xmax() && !other.is_open_xmin() && other.x_min - self.x_max > a_delta {
            return true;
        }
        if !self.is_open_ymin() && !other.is_open_ymax() && self.y_min - other.y_max > a_delta {
            return true;
        }
        if !self.is_open_ymax() && !other.is_open_ymin() && other.y_min - self.y_max > a_delta {
            return true;
        }
        if !self.is_open_zmin() && !other.is_open_zmax() && self.z_min - other.z_max > a_delta {
            return true;
        }
        if !self.is_open_zmax() && !other.is_open_zmin() && other.z_min - self.z_max > a_delta {
            return true;
        }
        false
    }

    // ── Contains — OCCT: Contains(Pnt) ─────────────────────────────────

    /// Test if a point is inside or on the boundary of this box.
    /// OCCT: Contains(const gp_Pnt& P) const is the inline `!IsOut(P)`
    /// (Bnd_Box.hxx).
    pub fn contains(&self, p: DVec3) -> bool {
        !self.is_out_point(p)
    }

    // ── Distance — OCCT: Distance(Bnd_Box) ─────────────────────────────

    /// Minimum Euclidean distance between this box and another.
    /// Returns 0 if they intersect or either is void.
    /// OCCT: Standard_Real Distance(const Bnd_Box& Other) const
    /// (Bnd_Box.cxx L1249-1268): the gap- and flag-aware Get() corners of
    /// both boxes feed the per-dimension squared distances; the sum is
    /// square-rooted once (OCCT's exact formula).
    pub fn distance(&self, other: &BndBox) -> f64 {
        if self.is_void() || other.is_void() {
            return 0.0;
        }
        let (axmin1, aymin1, azmin1, axmax1, aymax1, azmax1) = self.get().unwrap();
        let (axmin2, aymin2, azmin2, axmax2, aymax2, azmax2) = other.get().unwrap();
        let a_dist_x = distance_in_dimension(axmin1, axmax1, axmin2, axmax2);
        let a_dist_y = distance_in_dimension(aymin1, aymax1, aymin2, aymax2);
        let a_dist_z = distance_in_dimension(azmin1, azmax1, azmin2, azmax2);
        (a_dist_x + a_dist_y + a_dist_z).sqrt()
    }

    // ── Transform — OCCT: Transformed(Trsf) ────────────────────────────

    /// Return a transformed copy of this box (axis-aligned result, not OBB).
    /// OCCT: Bnd_Box Transformed(const gp_Trsf& T) const
    /// (Bnd_Box.cxx L411-508): identity returns the box as-is; otherwise the
    /// finite part (when any) contributes its 8 transformed corners, the gap
    /// is copied, and each open direction re-opens through the transformed
    /// direction (Add(const gp_Dir&), Bnd_Box.cxx L578-603).
    ///
    /// Architecture note: OCCT dispatches on gp_Trsf::Form() with a
    /// translation fast path that returns an open box untranslated; rcad
    /// carries a DAffine3 without a Form tag, so the general path is used
    /// (same construction, OCCT's fast-path shortcut for open boxes is
    /// unreachable by design here).
    pub fn transformed(&self, transform: &glam::DAffine3) -> Self {
        if self.is_void() { return Self::new(); }
        if *transform == glam::DAffine3::IDENTITY {
            return self.clone();
        }
        let mut out = Self::new();
        if self.has_finite_part() {
            let corners = [
                DVec3::new(self.x_min, self.y_min, self.z_min),
                DVec3::new(self.x_max, self.y_min, self.z_min),
                DVec3::new(self.x_min, self.y_max, self.z_min),
                DVec3::new(self.x_max, self.y_max, self.z_min),
                DVec3::new(self.x_min, self.y_min, self.z_max),
                DVec3::new(self.x_max, self.y_min, self.z_max),
                DVec3::new(self.x_min, self.y_max, self.z_max),
                DVec3::new(self.x_max, self.y_max, self.z_max),
            ];
            for &p in &corners {
                out.add_point(transform.transform_point3(p));
            }
        }
        out.gap = self.gap;
        if !self.is_open() {
            return out;
        }
        let dirs = [
            (self.is_open_xmin(), DVec3::new(-1.0, 0.0, 0.0)),
            (self.is_open_xmax(), DVec3::new(1.0, 0.0, 0.0)),
            (self.is_open_ymin(), DVec3::new(0.0, -1.0, 0.0)),
            (self.is_open_ymax(), DVec3::new(0.0, 1.0, 0.0)),
            (self.is_open_zmin(), DVec3::new(0.0, 0.0, -1.0)),
            (self.is_open_zmax(), DVec3::new(0.0, 0.0, 1.0)),
        ];
        for (open, dir) in dirs {
            if open {
                out.add_dir(transform.transform_vector3(dir).normalize_or_zero());
            }
        }
        out
    }

    /// OCCT Bnd_Box::HasFinitePart() const (Bnd_Box.hxx L374):
    /// `!IsVoid() && Xmax >= Xmin`.
    fn has_finite_part(&self) -> bool {
        !self.is_void() && self.x_max >= self.x_min
    }

    /// OCCT Bnd_Box::Add(const gp_Dir& D) (Bnd_Box.cxx L578-603): opens the
    /// directions a (unit) direction vector points toward, thresholded at
    /// the real epsilon (gp::RealEpsilon == DBL_EPSILON).
    fn add_dir(&mut self, d: DVec3) {
        if d.x < -f64::EPSILON {
            self.open_xmin();
        } else if d.x > f64::EPSILON {
            self.open_xmax();
        }
        if d.y < -f64::EPSILON {
            self.open_ymin();
        } else if d.y > f64::EPSILON {
            self.open_ymax();
        }
        if d.z < -f64::EPSILON {
            self.open_zmin();
        } else if d.z > f64::EPSILON {
            self.open_zmax();
        }
    }

    /// Clear the box (set to void).  OCCT: void Clear().
    pub fn clear(&mut self) {
        self.x_min = 0.0; self.x_max = 0.0;
        self.y_min = 0.0; self.y_max = 0.0;
        self.z_min = 0.0; self.z_max = 0.0;
        self.gap = 0.0;
        self.flags = VOID_MASK;
    }

    /// OCCT: void OpenXmin() / OpenXmax() etc.
    /// Mark a direction as open (infinite).
    pub fn open_xmin(&mut self) { self.flags |= XMIN_OPEN; }
    pub fn open_xmax(&mut self) { self.flags |= XMAX_OPEN; }
    pub fn open_ymin(&mut self) { self.flags |= YMIN_OPEN; }
    pub fn open_ymax(&mut self) { self.flags |= YMAX_OPEN; }
    pub fn open_zmin(&mut self) { self.flags |= ZMIN_OPEN; }
    pub fn open_zmax(&mut self) { self.flags |= ZMAX_OPEN; }

    // ── BVH query helpers (not in OCCT Bnd_Box, but needed for BVH) ─────

    /// Center of the box. Returns `DVec3::ZERO` if void.
    pub fn center(&self) -> DVec3 {
        if self.is_void() { return DVec3::ZERO; }
        DVec3::new(
            0.5 * (self.x_min + self.x_max),
            0.5 * (self.y_min + self.y_max),
            0.5 * (self.z_min + self.z_max),
        )
    }

    /// Surface area of the box (for SAH). Returns 0 if void.
    pub fn surface_area(&self) -> f64 {
        if self.is_void() { return 0.0; }
        let dx = self.x_max - self.x_min;
        let dy = self.y_max - self.y_min;
        let dz = self.z_max - self.z_min;
        2.0 * (dx * dy + dx * dz + dy * dz)
    }

    /// True if this box intersects another (inverse of is_out_box).
    pub fn intersects(&self, other: &BndBox) -> bool {
        !self.is_out_box(other)
    }

    /// Ray-box intersection using the slab method.
    /// `inv_dir` = 1.0 / ray_direction (component-wise; INF for zero components).
    /// Returns `Some(t)` where `t ≥ 0` is the entry distance, or `None`.
    pub fn ray_intersect(&self, origin: DVec3, inv_dir: DVec3) -> Option<f64> {
        if self.is_void() { return None; }
        let (xmin, ymin, zmin, xmax, ymax, zmax) = self.get()?;
        let t1 = (xmin - origin.x) * inv_dir.x;
        let t2 = (xmax - origin.x) * inv_dir.x;
        let t3 = (ymin - origin.y) * inv_dir.y;
        let t4 = (ymax - origin.y) * inv_dir.y;
        let t5 = (zmin - origin.z) * inv_dir.z;
        let t6 = (zmax - origin.z) * inv_dir.z;
        let t_enter = t1.min(t2).max(t3.min(t4)).max(t5.min(t6));
        let t_exit = t1.max(t2).min(t3.max(t4)).min(t5.max(t6));
        if t_enter <= t_exit && t_exit >= 0.0 { Some(t_enter.max(0.0)) } else { None }
    }

    /// Minimum squared distance from a point to this box.  Returns 0 if inside.
    pub fn point_dist_sq(&self, point: DVec3) -> f64 {
        if self.is_void() { return f64::INFINITY; }
        let g = self.gap;
        let dx = if point.x < self.x_min - g { self.x_min - g - point.x }
                 else if point.x > self.x_max + g { point.x - self.x_max - g }
                 else { 0.0 };
        let dy = if point.y < self.y_min - g { self.y_min - g - point.y }
                 else if point.y > self.y_max + g { point.y - self.y_max - g }
                 else { 0.0 };
        let dz = if point.z < self.z_min - g { self.z_min - g - point.z }
                 else if point.z > self.z_max + g { point.z - self.z_max - g }
                 else { 0.0 };
        dx * dx + dy * dy + dz * dz
    }

    /// Alias for `add_point` — extend to include a point.
    pub fn expand_point(&mut self, p: DVec3) { self.add_point(p); }

    /// Alias for `add_box` — extend to enclose another box.
    pub fn expand_aabb(&mut self, other: &BndBox) { self.add_box(other); }

    /// Raw minimum corner (without gap). For BVH SAH computations.
    pub fn raw_min(&self) -> DVec3 { DVec3::new(self.x_min, self.y_min, self.z_min) }

    /// Raw maximum corner (without gap). For BVH SAH computations.
    pub fn raw_max(&self) -> DVec3 { DVec3::new(self.x_max, self.y_max, self.z_max) }

    /// Raw x-axis span.
    pub fn dx(&self) -> f64 { self.x_max - self.x_min }

    /// Raw y-axis span.
    pub fn dy(&self) -> f64 { self.y_max - self.y_min }

    /// Raw z-axis span.
    pub fn dz(&self) -> f64 { self.z_max - self.z_min }
}

impl Default for BndBox {
    fn default() -> Self { Self::new() }
}

// OCCT Bnd_Box.cxx L45-52 — DistMini2Box: the minimum squared distance
// between two 1D intervals.
fn dist_mini_2_box(r1_min: f64, r1_max: f64, r2_min: f64, r2_max: f64) -> f64 {
    let a_r1 = (r1_min - r2_max) * (r1_min - r2_max);
    let a_r2 = (r1_max - r2_min) * (r1_max - r2_min);
    a_r1.min(a_r2)
}

// OCCT Bnd_Box.cxx L53-64 — DistanceInDimension: the squared distance in one
// dimension, 0 when the intervals overlap.
fn distance_in_dimension(min1: f64, max1: f64, min2: f64, max2: f64) -> f64 {
    if (min1 <= min2 && min2 <= max1) || (min2 <= min1 && min1 <= max2) {
        0.0
    } else {
        dist_mini_2_box(min1, max1, min2, max2)
    }
}

// ── OCCT Bnd_Box2d (Bnd_Box2d.cxx) — the 2D axis-aligned box ──────────────

/// OCCT Bnd_Box2d — axis-aligned bounding box in 2D with a gap (tolerance).
///
/// State flags mirror Bnd_Box2d.hxx: VoidMask + the four open directions.
#[derive(Debug, Clone)]
pub struct BndBox2d {
    x_min: f64,
    y_min: f64,
    x_max: f64,
    y_max: f64,
    gap: f64,
    flags: u8,
}

// OCCT Bnd_Box2d.hxx flags
const VOID2D_MASK: u8 = 1;
const XMIN2D_OPEN: u8 = 2;
const XMAX2D_OPEN: u8 = 4;
const YMIN2D_OPEN: u8 = 8;
const YMAX2D_OPEN: u8 = 16;

impl BndBox2d {
    /// Default constructor — a void (uninitialised) box.
    /// OCCT: Bnd_Box2d() → myFlags = VoidMask.
    pub fn new() -> Self {
        BndBox2d {
            x_min: 0.0,
            y_min: 0.0,
            x_max: 0.0,
            y_max: 0.0,
            gap: 0.0,
            flags: VOID2D_MASK,
        }
    }

    /// OCCT Bnd_Box2d::Update(aXmin, aYmin, aXmax, aYmax) — Bnd_Box2d.cxx
    /// L39-70: a void box takes the rectangle (VoidMask cleared), otherwise
    /// each still-closed direction takes the min/max.
    pub fn update(&mut self, x: f64, y: f64, x_max: f64, y_max: f64) {
        if (self.flags & VOID2D_MASK) != 0 {
            self.x_min = x;
            self.y_min = y;
            self.x_max = x_max;
            self.y_max = y_max;
            self.flags &= !VOID2D_MASK;
        } else {
            if (self.flags & XMIN2D_OPEN) == 0 {
                self.x_min = self.x_min.min(x);
            }
            if (self.flags & XMAX2D_OPEN) == 0 {
                self.x_max = self.x_max.max(x_max);
            }
            if (self.flags & YMIN2D_OPEN) == 0 {
                self.y_min = self.y_min.min(y);
            }
            if (self.flags & YMAX2D_OPEN) == 0 {
                self.y_max = self.y_max.max(y_max);
            }
        }
    }

    /// OCCT Bnd_Box2d::Update(X, Y) — Bnd_Box2d.cxx L72-100.
    pub fn update_xy(&mut self, x: f64, y: f64) {
        if (self.flags & VOID2D_MASK) != 0 {
            self.x_min = x;
            self.y_min = y;
            self.x_max = x;
            self.y_max = y;
            self.flags &= !VOID2D_MASK;
        } else {
            if (self.flags & XMIN2D_OPEN) == 0 {
                self.x_min = self.x_min.min(x);
            }
            if (self.flags & XMAX2D_OPEN) == 0 {
                self.x_max = self.x_max.max(x);
            }
            if (self.flags & YMIN2D_OPEN) == 0 {
                self.y_min = self.y_min.min(y);
            }
            if (self.flags & YMAX2D_OPEN) == 0 {
                self.y_max = self.y_max.max(y);
            }
        }
    }

    /// OCCT Bnd_Box2d::SetVoid() — the box becomes void.
    pub fn set_void(&mut self) {
        self.flags = VOID2D_MASK;
    }

    /// OCCT Bnd_Box2d::SetWhole() — the whole plane (all directions open).
    pub fn set_whole(&mut self) {
        self.flags = XMIN2D_OPEN | XMAX2D_OPEN | YMIN2D_OPEN | YMAX2D_OPEN;
    }

    /// OCCT Bnd_Box2d::IsVoid().
    pub fn is_void(&self) -> bool {
        self.flags & VOID2D_MASK != 0
    }

    /// OCCT Bnd_Box2d::IsWhole() — all four directions open.
    pub fn is_whole(&self) -> bool {
        self.flags & (XMIN2D_OPEN | XMAX2D_OPEN | YMIN2D_OPEN | YMAX2D_OPEN)
            == (XMIN2D_OPEN | XMAX2D_OPEN | YMIN2D_OPEN | YMAX2D_OPEN)
    }

    /// OCCT Bnd_Box2d::OpenXmin() (Bnd_Box2d.hxx L164).
    pub fn open_xmin(&mut self) {
        self.flags |= XMIN2D_OPEN;
    }
    /// OCCT Bnd_Box2d::OpenXmax() (Bnd_Box2d.hxx L167).
    pub fn open_xmax(&mut self) {
        self.flags |= XMAX2D_OPEN;
    }
    /// OCCT Bnd_Box2d::OpenYmin() (Bnd_Box2d.hxx L170).
    pub fn open_ymin(&mut self) {
        self.flags |= YMIN2D_OPEN;
    }
    /// OCCT Bnd_Box2d::OpenYmax() (Bnd_Box2d.hxx L173).
    pub fn open_ymax(&mut self) {
        self.flags |= YMAX2D_OPEN;
    }

    /// OCCT Bnd_Box2d::IsOpenXmin/Xmax/Ymin/Ymax().
    pub fn is_open_xmin(&self) -> bool {
        self.flags & XMIN2D_OPEN != 0
    }
    pub fn is_open_xmax(&self) -> bool {
        self.flags & XMAX2D_OPEN != 0
    }
    pub fn is_open_ymin(&self) -> bool {
        self.flags & YMIN2D_OPEN != 0
    }
    pub fn is_open_ymax(&self) -> bool {
        self.flags & YMAX2D_OPEN != 0
    }

    /// Current gap (tolerance).  OCCT: GetGap().
    pub fn get_gap(&self) -> f64 {
        self.gap
    }

    /// Set the gap.  OCCT: SetGap(Tol).
    pub fn set_gap(&mut self, tol: f64) {
        self.gap = tol.abs();
    }

    /// OCCT Bnd_Box2d::Enlarge(Tol) (Bnd_Box2d.hxx L127):
    /// `Gap = max(Gap, abs(Tol))` — the raw bounds are left alone; the
    /// gap-applying getters see the enlargement exactly once.
    pub fn enlarge(&mut self, tol: f64) {
        self.gap = self.gap.max(tol.abs());
    }

    /// OCCT Bnd_Box2d::GetXMin() — Bnd_Box2d.cxx L127-130.
    fn get_x_min(&self) -> f64 {
        if self.flags & XMIN2D_OPEN != 0 {
            -THE_BND_PRECISION_INFINITE
        } else {
            self.x_min - self.gap
        }
    }

    /// OCCT Bnd_Box2d::GetXMax() — Bnd_Box2d.cxx L134-137.
    fn get_x_max(&self) -> f64 {
        if self.flags & XMAX2D_OPEN != 0 {
            THE_BND_PRECISION_INFINITE
        } else {
            self.x_max + self.gap
        }
    }

    /// OCCT Bnd_Box2d::GetYMin() — Bnd_Box2d.cxx L141-144.
    fn get_y_min(&self) -> f64 {
        if self.flags & YMIN2D_OPEN != 0 {
            -THE_BND_PRECISION_INFINITE
        } else {
            self.y_min - self.gap
        }
    }

    /// OCCT Bnd_Box2d::GetYMax() — Bnd_Box2d.cxx L148-151.
    fn get_y_max(&self) -> f64 {
        if self.flags & YMAX2D_OPEN != 0 {
            THE_BND_PRECISION_INFINITE
        } else {
            self.y_max + self.gap
        }
    }

    /// OCCT Bnd_Box2d::Get(xmin, ymin, xmax, ymax) — Bnd_Box2d.cxx
    /// L105-116: throws for a void box; otherwise the finite corners through
    /// GetXMin/GetXMax/GetYMin/GetYMax (gap applied, open directions
    /// infinite).  `None` models the Standard_ConstructionError.
    pub fn get(&self) -> Option<(f64, f64, f64, f64)> {
        if self.is_void() {
            return None;
        }
        Some((
            self.get_x_min(),
            self.get_y_min(),
            self.get_x_max(),
            self.get_y_max(),
        ))
    }

    /// OCCT Bnd_Box2d::Add(const gp_Pnt2d& thePnt) (Bnd_Box2d.hxx L206) —
    /// the header inlines Add to Update(thePnt.X(), thePnt.Y()).
    pub fn add_point(&mut self, p: glam::DVec2) {
        self.update_xy(p.x, p.y);
    }

    /// OCCT Bnd_Box2d::Add(const Bnd_Box2d& Other) — Bnd_Box2d.cxx L258-330.
    pub fn add_box(&mut self, other: &BndBox2d) {
        if self.is_whole() {
            return;
        } else if other.is_void() {
            return;
        } else if other.is_whole() {
            self.set_whole();
        } else if self.is_void() {
            *self = other.clone();
        } else {
            if !self.is_open_xmin() {
                if other.is_open_xmin() {
                    self.open_xmin();
                } else if self.x_min > other.x_min {
                    self.x_min = other.x_min;
                }
            }
            if !self.is_open_xmax() {
                if other.is_open_xmax() {
                    self.open_xmax();
                } else if self.x_max < other.x_max {
                    self.x_max = other.x_max;
                }
            }
            if !self.is_open_ymin() {
                if other.is_open_ymin() {
                    self.open_ymin();
                } else if self.y_min > other.y_min {
                    self.y_min = other.y_min;
                }
            }
            if !self.is_open_ymax() {
                if other.is_open_ymax() {
                    self.open_ymax();
                } else if self.y_max < other.y_max {
                    self.y_max = other.y_max;
                }
            }
            self.gap = self.gap.max(other.gap);
        }
    }

    /// OCCT Bnd_Box2d::IsOut(const gp_Pnt2d& P) — Bnd_Box2d.cxx L354-390.
    /// A whole box is out for no point, a void box for every point, and the
    /// per-direction open flags suppress the corresponding test.
    pub fn is_out_point(&self, the_p: glam::DVec2) -> bool {
        if self.is_whole() {
            return false;
        }
        if self.is_void() {
            return true;
        }
        if (self.flags & XMIN2D_OPEN) == 0 && the_p.x < (self.x_min - self.gap) {
            return true;
        }
        if (self.flags & XMAX2D_OPEN) == 0 && the_p.x > (self.x_max + self.gap) {
            return true;
        }
        if (self.flags & YMIN2D_OPEN) == 0 && the_p.y < (self.y_min - self.gap) {
            return true;
        }
        if (self.flags & YMAX2D_OPEN) == 0 && the_p.y > (self.y_max + self.gap) {
            return true;
        }
        false
    }

    /// OCCT Bnd_Box2d::IsOut(const gp_Lin2d& theL) — Bnd_Box2d.cxx L393-416.
    /// The signed area of the parallelogram (direction, box-center offset)
    /// against the box half-extents projected on the line direction.
    pub fn is_out_line(&self, the_l: &crate::geom::Line2d) -> bool {
        if self.is_whole() {
            return false;
        }
        if self.is_void() {
            return true;
        }
        let Some((a_x_min, a_y_min, a_x_max, a_y_max)) = self.get() else {
            return true;
        };

        let a_center = glam::DVec2::new((a_x_min + a_x_max) / 2.0, (a_y_min + a_y_max) / 2.0);
        let a_heigh = glam::DVec2::new(
            (a_x_max - a_center.x).abs(),
            (a_y_max - a_center.y).abs(),
        );

        let a_dir = the_l.direction;
        let a_loc = the_l.origin;
        // gp_XY::operator^ is the 2D cross product X1*Y2 - Y1*X2.
        let a_prod = [
            a_dir.x * (a_center.y - a_loc.y) - a_dir.y * (a_center.x - a_loc.x),
            a_dir.x * a_heigh.y,
            a_dir.y * a_heigh.x,
        ];
        a_prod[0].abs() > (a_prod[1].abs() + a_prod[2].abs())
    }

    /// OCCT Bnd_Box2d::IsOut(const gp_Pnt2d& theP0, const gp_Pnt2d& theP1) —
    /// Bnd_Box2d.cxx L418-454: the segment-vs-box rejection used by
    /// BRepClass_Intersector::IsInter (BRepClass_Intersector.cxx L122-136).
    pub fn is_out_segment(&self, the_p0: glam::DVec2, the_p1: glam::DVec2) -> bool {
        if self.is_whole() {
            return false;
        }
        if self.is_void() {
            return true;
        }
        let Some((a_loc_x_min, a_loc_y_min, a_loc_x_max, a_loc_y_max)) = self.get() else {
            return true;
        };

        // Intersect the line containing the segment.
        let a_seg_delta = the_p1 - the_p0;

        let a_center = glam::DVec2::new(
            (a_loc_x_min + a_loc_x_max) / 2.0,
            (a_loc_y_min + a_loc_y_max) / 2.0,
        );
        let a_heigh = glam::DVec2::new(
            (a_loc_x_max - a_center.x).abs(),
            (a_loc_y_max - a_center.y).abs(),
        );

        let a_prod = [
            a_seg_delta.x * (a_center.y - the_p0.y) - a_seg_delta.y * (a_center.x - the_p0.x),
            a_seg_delta.x * a_heigh.y,
            a_seg_delta.y * a_heigh.x,
        ];

        if a_prod[0].abs() <= (a_prod[1].abs() + a_prod[2].abs()) {
            // Intersection with the line detected; check the segment as a
            // bounding box around its own center.
            let a_h_seg = glam::DVec2::new(0.5 * a_seg_delta.x, 0.5 * a_seg_delta.y);
            let a_h_seg_abs = glam::DVec2::new(a_h_seg.x.abs(), a_h_seg.y.abs());
            let a_mid = the_p0 + a_h_seg - a_center;
            return a_mid.x.abs() > (a_heigh.x + a_h_seg_abs.x)
                || a_mid.y.abs() > (a_heigh.y + a_h_seg_abs.y);
        }
        true
    }

    /// OCCT Bnd_Box2d::IsOut(const Bnd_Box2d& Other) — Bnd_Box2d.cxx
    /// L456-511: fast path for non-open/non-void/non-whole boxes, then the
    /// per-flag general path.
    pub fn is_out_box(&self, other: &BndBox2d) -> bool {
        // Fast path for non-open, non-void, non-whole boxes.
        if self.flags == 0 && other.flags == 0 {
            let a_delta = other.gap + self.gap;
            if self.x_min - other.x_max > a_delta {
                return true;
            }
            if other.x_min - self.x_max > a_delta {
                return true;
            }
            if self.y_min - other.y_max > a_delta {
                return true;
            }
            if other.y_min - self.y_max > a_delta {
                return true;
            }
            return false;
        }

        // Handle special cases.
        if self.is_void() || other.is_void() {
            return true;
        }
        if self.is_whole() || other.is_whole() {
            return false;
        }

        let Some((oxmin, oymin, oxmax, oymax)) = other.get() else {
            return true;
        };
        if (self.flags & XMIN2D_OPEN) == 0 && oxmax < (self.x_min - self.gap) {
            return true;
        }
        if (self.flags & XMAX2D_OPEN) == 0 && oxmin > (self.x_max + self.gap) {
            return true;
        }
        if (self.flags & YMIN2D_OPEN) == 0 && oymax < (self.y_min - self.gap) {
            return true;
        }
        if (self.flags & YMAX2D_OPEN) == 0 && oymin > (self.y_max + self.gap) {
            return true;
        }
        false
    }
}

impl Default for BndBox2d {
    fn default() -> Self {
        Self::new()
    }
}



// ══════════════════════════════════════════════════════════════════════════
// Tests (OCCT-aligned: Bnd_Box_Test.cxx)
// ══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_box_is_void() {
        let b = BndBox::new();
        assert!(b.is_void());
        assert!(b.is_out_point(DVec3::ZERO));
    }

    #[test]
    fn from_point_contains_that_point() {
        let p = DVec3::new(1.0, 2.0, 3.0);
        let b = BndBox::from_point(p);
        assert!(!b.is_void());
        assert!(b.contains(p));
        assert!(!b.is_out_point(p));
    }

    #[test]
    fn add_point_expands_box() {
        let mut b = BndBox::from_point(DVec3::ZERO);
        b.add_point(DVec3::new(1.0, 1.0, 1.0));
        assert!(b.contains(DVec3::new(0.5, 0.5, 0.5)));
        assert!(!b.is_out_point(DVec3::new(0.5, 0.5, 0.5)));
        assert!(b.is_out_point(DVec3::new(2.0, 2.0, 2.0)));
    }

    #[test]
    fn is_out_box_detects_separation() {
        let a = BndBox::from_corners(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = BndBox::from_corners(2.0, 0.0, 0.0, 3.0, 1.0, 1.0);
        assert!(a.is_out_box(&b));
        assert!(b.is_out_box(&a));
    }

    #[test]
    fn is_out_box_false_for_overlapping() {
        let a = BndBox::from_corners(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let b = BndBox::from_corners(1.0, 1.0, 1.0, 3.0, 3.0, 3.0);
        assert!(!a.is_out_box(&b));
    }

    #[test]
    fn distance_between_separated_boxes() {
        let a = BndBox::from_corners(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = BndBox::from_corners(4.0, 0.0, 0.0, 5.0, 1.0, 1.0);
        let d = a.distance(&b);
        assert!((d - 3.0).abs() < 1e-12, "distance={}", d);
    }

    #[test]
    fn gap_expands_box() {
        let mut b = BndBox::from_point(DVec3::ZERO);
        b.set_gap(1.0);
        assert!(b.contains(DVec3::new(0.9, 0.9, 0.9)));
        assert!(b.is_out_point(DVec3::new(1.1, 0.0, 0.0)));
    }

    /// OCCT Bnd_Box::Enlarge raises the gap only (hxx L154): the raw bounds
    /// stay untouched and the getter applies the gap exactly once.
    #[test]
    fn enlarge_raises_gap_only() {
        let mut b = BndBox::from_corners(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        b.set_gap(0.5);
        b.enlarge(2.0);
        assert_eq!(b.get_gap(), 2.0);
        let (xmin, ymin, zmin, xmax, ymax, zmax) = b.get().unwrap();
        assert_eq!(xmin, -2.0);
        assert_eq!(ymin, -2.0);
        assert_eq!(zmin, -2.0);
        assert_eq!(xmax, 3.0);
        assert_eq!(ymax, 3.0);
        assert_eq!(zmax, 3.0);
        // A smaller tolerance does not shrink the gap.
        b.enlarge(1.0);
        assert_eq!(b.get_gap(), 2.0);
    }

    /// OCCT Bnd_Box2d::Enlarge mirrors the 3D gap-only semantics (hxx L127).
    #[test]
    fn enlarge2d_raises_gap_only() {
        let mut b = BndBox2d::new();
        b.update(0.0, 0.0, 1.0, 1.0);
        b.set_gap(0.5);
        b.enlarge(2.0);
        assert_eq!(b.get_gap(), 2.0);
        assert_eq!(b.get(), Some((-2.0, -2.0, 3.0, 3.0)));
    }

    #[test]
    fn transformed_box_still_contains_corners() {
        let b = BndBox::from_corners(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let t = glam::DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0));
        let bt = b.transformed(&t);
        assert!(bt.contains(DVec3::new(10.5, 0.5, 0.5)));
    }

    #[test]
    fn add_box_union() {
        let mut a = BndBox::from_corners(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = BndBox::from_corners(2.0, 2.0, 2.0, 3.0, 3.0, 3.0);
        a.add_box(&b);
        assert!(a.contains(DVec3::new(0.5, 0.5, 0.5)));
        assert!(a.contains(DVec3::new(2.5, 2.5, 2.5)));
    }

    #[test]
    fn get_returns_corners_with_gap() {
        let mut b = BndBox::from_corners(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        b.set_gap(0.5);
        let (xmin, ymin, zmin, xmax, ymax, zmax) = b.get().unwrap();
        assert!((xmin + 0.5).abs() < 1e-12);
        assert!((xmax - 1.5).abs() < 1e-12);
    }

    #[test]
    fn clear_resets_to_void() {
        let mut b = BndBox::from_point(DVec3::new(1.0, 2.0, 3.0));
        assert!(!b.is_void());
        b.clear();
        assert!(b.is_void());
    }

    #[test]
    fn open_flags_work() {
        let mut b = BndBox::from_point(DVec3::ZERO);
        assert!(!b.is_open());
        b.open_xmax();
        assert!(b.is_open());
    }

    #[test]
    fn curve_circle_bounding_box() {
        use crate::geom::{Circle3, Curve3};
        let c = Curve3::Circle(Circle3::new(DVec3::ZERO, DVec3::Z, 2.0));
        let bb = crate::curve_bounding_box(&c).unwrap();
        assert!((bb[0].x + 2.0).abs() < 1e-12);
        assert!((bb[1].x - 2.0).abs() < 1e-12);
        assert!((bb[0].z).abs() < 1e-12); // circle in XY plane
        assert!((bb[1].z).abs() < 1e-12);
    }
}

// Bnd_Range tests (OCCT Bnd_Range.hxx/.cxx anchors) live in range.rs.
#[cfg(test)]
mod range_tests {
    use super::range::{IntersectStatus, Range};

    /// OCCT Bnd_Range void semantics: default ctor is VOID; Delta() negative;
    /// Get* return empty.
    #[test]
    fn void_range_semantics() {
        let r = Range::new();
        assert!(r.is_void());
        assert!(r.delta() < 0.0);
        assert_eq!(r.get_min(), None);
        assert_eq!(r.get_bounds(), None);
        assert_eq!(r.center(), None);
        // Add revives a void range to a single point.
        let mut r2 = Range::new();
        r2.add_parameter(3.5);
        assert_eq!(r2.get_bounds(), Some((3.5, 3.5)));
    }

    /// OCCT Bnd_Range::Common / Union / Add (cxx L21-61, hxx L114-127).
    #[test]
    fn common_union_add() {
        let mut a = Range::from_bounds(0.0, 10.0);
        a.common(&Range::from_bounds(5.0, 20.0));
        assert_eq!(a.get_bounds(), Some((5.0, 10.0)));

        let mut b = Range::from_bounds(8.0, 12.0);
        assert!(b.union_with(&Range::from_bounds(11.0, 15.0)));
        assert_eq!(b.get_bounds(), Some((8.0, 15.0)));
        // Separated ranges cannot be united.
        assert!(!b.union_with(&Range::from_bounds(20.0, 25.0)));
        // Add merges unconditionally.
        b.add_range(&Range::from_bounds(20.0, 25.0));
        assert_eq!(b.get_bounds(), Some((8.0, 25.0)));
    }

    /// OCCT Bnd_Range::IsIntersected non-periodic (cxx L65-90): Boundary at
    /// the ends, In strictly inside, Out outside.
    #[test]
    fn is_intersected_statuses() {
        let r = Range::from_bounds(0.0, 10.0);
        assert_eq!(r.is_intersected(5.0, 0.0), IntersectStatus::In);
        assert_eq!(r.is_intersected(0.0, 0.0), IntersectStatus::Boundary);
        assert_eq!(r.is_intersected(10.0, 0.0), IntersectStatus::Boundary);
        assert_eq!(r.is_intersected(11.0, 0.0), IntersectStatus::Out);
        // Periodic: 12 == 2 mod 10 lies inside the shifted lattice.
        assert_eq!(r.is_intersected(12.0, 10.0), IntersectStatus::In);
        assert_eq!(r.is_intersected(20.0, 10.0), IntersectStatus::Boundary);
    }

    /// OCCT Bnd_Range::Split (cxx L132-171): [3,15] by 5 -> [3,5],[5,15];
    /// periodic split of [3,15] by value 5 period 4 -> [3,5],[5,9],[9,13],[13,15].
    #[test]
    fn split_simple_and_periodic() {
        let r = Range::from_bounds(3.0, 15.0);
        let mut parts = Vec::new();
        r.split(5.0, &mut parts, 0.0);
        assert_eq!(parts.len(), 2);
        assert_eq!(parts[0].get_bounds(), Some((3.0, 5.0)));
        assert_eq!(parts[1].get_bounds(), Some((5.0, 15.0)));

        let mut pparts = Vec::new();
        r.split(5.0, &mut pparts, 4.0);
        assert_eq!(pparts.len(), 4);
        assert_eq!(pparts[0].get_bounds(), Some((3.0, 5.0)));
        assert_eq!(pparts[1].get_bounds(), Some((5.0, 9.0)));
        assert_eq!(pparts[2].get_bounds(), Some((9.0, 13.0)));
        assert_eq!(pparts[3].get_bounds(), Some((13.0, 15.0)));
    }
}
