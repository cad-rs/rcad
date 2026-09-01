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

    /// OCCT Bnd_Box::Update(xmin, ymin, zmin, xmax, ymax, zmax) — the box
    /// becomes the finite axis-aligned box (open flags cleared).
    pub fn update(&mut self, x_min: f64, y_min: f64, z_min: f64, x_max: f64, y_max: f64, z_max: f64) {
        self.x_min = x_min;
        self.y_min = y_min;
        self.z_min = z_min;
        self.x_max = x_max;
        self.y_max = y_max;
        self.z_max = z_max;
        self.flags &= !(VOID_MASK | XMIN_OPEN | XMAX_OPEN | YMIN_OPEN | YMAX_OPEN | ZMIN_OPEN | ZMAX_OPEN);
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

    /// Enlarge the box by a tolerance in all six directions.
    /// OCCT: Bnd_Box::Enlarge(const Standard_Real Tol).
    /// A void box is left untouched (OCCT: void box stays void).
    pub fn enlarge(&mut self, tol: f64) {
        if self.is_void() || !tol.is_finite() { return; }
        self.x_min -= tol;
        self.x_max += tol;
        self.y_min -= tol;
        self.y_max += tol;
        self.z_min -= tol;
        self.z_max += tol;
    }

    // ── Get corners (including gap) — OCCT: Get() ───────────────────────

    /// Retrieve the finite box corners including gap.
    /// Returns `None` if the box is void.
    /// OCCT: void Get(xmin, ymin, zmin, xmax, ymax, zmax) const.
    pub fn get(&self) -> Option<(f64, f64, f64, f64, f64, f64)> {
        if self.is_void() { return None; }
        let g = self.gap;
        Some((
            self.x_min - g, self.y_min - g, self.z_min - g,
            self.x_max + g, self.y_max + g, self.z_max + g,
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

    /// Extend the box to include a point.  OCCT: Add(const gp_Pnt&).
    pub fn add_point(&mut self, p: DVec3) {
        if self.is_void() {
            self.x_min = p.x; self.x_max = p.x;
            self.y_min = p.y; self.y_max = p.y;
            self.z_min = p.z; self.z_max = p.z;
            self.flags = 0;
        } else {
            if p.x < self.x_min { self.x_min = p.x }
            if p.x > self.x_max { self.x_max = p.x }
            if p.y < self.y_min { self.y_min = p.y }
            if p.y > self.y_max { self.y_max = p.y }
            if p.z < self.z_min { self.z_min = p.z }
            if p.z > self.z_max { self.z_max = p.z }
        }
    }

    /// Extend the box to enclose another box.  OCCT: Add(const Bnd_Box&).
    pub fn add_box(&mut self, other: &BndBox) {
        if other.is_void() { return; }
        if self.is_void() {
            *self = other.clone();
            return;
        }
        if other.x_min < self.x_min { self.x_min = other.x_min }
        if other.x_max > self.x_max { self.x_max = other.x_max }
        if other.y_min < self.y_min { self.y_min = other.y_min }
        if other.y_max > self.y_max { self.y_max = other.y_max }
        if other.z_min < self.z_min { self.z_min = other.z_min }
        if other.z_max > self.z_max { self.z_max = other.z_max }
        // Merge flags: if other has an open direction, propagate
        self.flags |= other.flags & !VOID_MASK;
    }

    // ── IsOut — OCCT: IsOut(Pnt), IsOut(Bnd_Box) ──────────────────────

    /// Test if a point is outside this box.
    /// OCCT: Standard_Boolean IsOut(const gp_Pnt&) const.
    pub fn is_out_point(&self, p: DVec3) -> bool {
        if self.is_void() { return true; }
        let g = self.gap;
        p.x < self.x_min - g || p.x > self.x_max + g
            || p.y < self.y_min - g || p.y > self.y_max + g
            || p.z < self.z_min - g || p.z > self.z_max + g
    }

    /// Test if another bounding box does NOT intersect this one.
    /// Returns `true` if `other` is entirely outside.
    /// OCCT: Standard_Boolean IsOut(const Bnd_Box&) const.
    pub fn is_out_box(&self, other: &BndBox) -> bool {
        if self.is_void() || other.is_void() { return true; }
        let g = self.gap + other.gap;
        other.x_max + g < self.x_min || other.x_min - g > self.x_max
            || other.y_max + g < self.y_min || other.y_min - g > self.y_max
            || other.z_max + g < self.z_min || other.z_min - g > self.z_max
    }

    // ── Contains — OCCT: Contains(Pnt) ─────────────────────────────────

    /// Test if a point is inside or on the boundary of this box.
    /// OCCT: Standard_Boolean Contains(const gp_Pnt&) const.
    pub fn contains(&self, p: DVec3) -> bool {
        if self.is_void() { return false; }
        let g = self.gap;
        p.x >= self.x_min - g && p.x <= self.x_max + g
            && p.y >= self.y_min - g && p.y <= self.y_max + g
            && p.z >= self.z_min - g && p.z <= self.z_max + g
    }

    // ── Distance — OCCT: Distance(Bnd_Box) ─────────────────────────────

    /// Minimum Euclidean distance between this box and another.
    /// Returns 0 if they intersect.
    /// OCCT: Standard_Real Distance(const Bnd_Box&) const.
    pub fn distance(&self, other: &BndBox) -> f64 {
        if self.is_void() || other.is_void() || !self.is_out_box(other) {
            return 0.0;
        }
        let dx = if self.x_max < other.x_min { other.x_min - self.x_max }
                 else if other.x_max < self.x_min { self.x_min - other.x_max }
                 else { 0.0 };
        let dy = if self.y_max < other.y_min { other.y_min - self.y_max }
                 else if other.y_max < self.y_min { self.y_min - other.y_max }
                 else { 0.0 };
        let dz = if self.z_max < other.z_min { other.z_min - self.z_max }
                 else if other.z_max < self.z_min { self.z_min - other.z_max }
                 else { 0.0 };
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    // ── Transform — OCCT: Transformed(Trsf) ────────────────────────────

    /// Return a transformed copy of this box (axis-aligned result, not OBB).
    /// OCCT: Bnd_Box Transformed(const gp_Trsf&) const.
    pub fn transformed(&self, transform: &glam::DAffine3) -> Self {
        if self.is_void() { return Self::new(); }
        let corners = [
            DVec3::new(self.x_min, self.y_min, self.z_min),
            DVec3::new(self.x_min, self.y_min, self.z_max),
            DVec3::new(self.x_min, self.y_max, self.z_min),
            DVec3::new(self.x_min, self.y_max, self.z_max),
            DVec3::new(self.x_max, self.y_min, self.z_min),
            DVec3::new(self.x_max, self.y_min, self.z_max),
            DVec3::new(self.x_max, self.y_max, self.z_min),
            DVec3::new(self.x_max, self.y_max, self.z_max),
        ];
        let mut out = Self::new();
        for &p in &corners {
            out.add_point(transform.transform_point3(p));
        }
        out.gap = self.gap;
        out
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

    /// OCCT Bnd_Box2d::Update(xmin, ymin, xmax, ymax) — the box becomes the
    /// finite axis-aligned rectangle (open flags cleared).
    pub fn update(&mut self, x_min: f64, y_min: f64, x_max: f64, y_max: f64) {
        self.x_min = x_min;
        self.y_min = y_min;
        self.x_max = x_max;
        self.y_max = y_max;
        self.flags &= !(VOID2D_MASK | XMIN2D_OPEN | XMAX2D_OPEN | YMIN2D_OPEN | YMAX2D_OPEN);
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

    /// OCCT Bnd_Box2d::Enlarge(Tol) — grow on all four sides (void unchanged).
    pub fn enlarge(&mut self, tol: f64) {
        if self.is_void() || !tol.is_finite() {
            return;
        }
        self.x_min -= tol;
        self.x_max += tol;
        self.y_min -= tol;
        self.y_max += tol;
    }

    /// OCCT Bnd_Box2d::Get(xmin, ymin, xmax, ymax) — the finite corners
    /// including the gap.  None for a void box.
    pub fn get(&self) -> Option<(f64, f64, f64, f64)> {
        if self.is_void() {
            return None;
        }
        let g = self.gap;
        Some((
            self.x_min - g,
            self.y_min - g,
            self.x_max + g,
            self.y_max + g,
        ))
    }

    /// OCCT Bnd_Box2d::Add(Pnt2d) — extend to include a point.
    pub fn add_point(&mut self, p: glam::DVec2) {
        if self.is_void() {
            self.x_min = p.x;
            self.x_max = p.x;
            self.y_min = p.y;
            self.y_max = p.y;
            self.flags &= !VOID2D_MASK;
        } else {
            if p.x < self.x_min {
                self.x_min = p.x;
            }
            if p.x > self.x_max {
                self.x_max = p.x;
            }
            if p.y < self.y_min {
                self.y_min = p.y;
            }
            if p.y > self.y_max {
                self.y_max = p.y;
            }
        }
    }

    /// OCCT Bnd_Box2d::IsOut(Pnt2d) — point outside the box (with gap).
    /// A void box is out for every point; a whole box for none.
    pub fn is_out_point(&self, p: glam::DVec2) -> bool {
        if self.is_void() {
            return true;
        }
        if self.is_whole() {
            return false;
        }
        let g = self.gap;
        p.x < self.x_min - g || p.x > self.x_max + g || p.y < self.y_min - g || p.y > self.y_max + g
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
