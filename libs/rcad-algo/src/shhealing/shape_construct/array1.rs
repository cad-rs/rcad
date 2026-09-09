//! Architecture bridge: OCCT `NCollection_Array1<T>` re-host.
//!
//! OCCT TKShHealing code stores per-point working arrays in
//! `NCollection_Array1` with explicit bounds (`ArrayOfPnt`, `ArrayOfPnt2d`,
//! `ArrayOfReal`, `CacheArray` in
//! `ShapeConstruct_ProjectCurveOnSurface.hxx` L55-61).  The rcad kernel has no
//! equivalent bounded-array foundation type, so this module provides the
//! minimal Vec-backed, explicitly-bounded array the 1:1 translation indexes
//! through.  It is a foundation (TKernel-level) bridge, not a TKShHealing
//! concept; indexing stays OCCT-form (`Value(i)` / `SetValue(i, v)` relative
//! to the array's own lower bound).

/// Minimal re-host of `NCollection_Array1<T>`: a contiguous run
/// `[lower, upper]` where `upper = lower + length - 1`.
#[derive(Debug, Clone)]
pub struct Array1<T> {
    lower: usize,
    data: Vec<T>,
}

impl<T: Clone> Array1<T> {
    /// OCCT `NCollection_Array1(theLower, theUpper)` — every element is
    /// default-constructed in OCCT; `the_fill` supplies the default here.
    pub fn new(lower: usize, upper: usize, the_fill: T) -> Self {
        let len = if upper + 1 >= lower {
            upper + 1 - lower
        } else {
            0
        };
        Array1 {
            lower,
            data: vec![the_fill; len],
        }
    }

    /// OCCT default constructor — empty array with lower bound 1.
    pub fn empty() -> Self {
        Array1 {
            lower: 1,
            data: Vec::new(),
        }
    }

    /// OCCT `Lower()`.
    pub fn lower(&self) -> usize {
        self.lower
    }

    /// OCCT `Upper()`.
    pub fn upper(&self) -> usize {
        self.lower + self.data.len() - 1
    }

    /// OCCT `Length()`.
    pub fn length(&self) -> usize {
        self.data.len()
    }

    /// OCCT `IsEmpty()`.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// OCCT `Value(theIndex)` — indexed relative to the array's own bounds.
    pub fn value(&self, the_index: usize) -> &T {
        &self.data[the_index - self.lower]
    }

    /// OCCT `SetValue(theIndex, theValue)`.
    pub fn set_value(&mut self, the_index: usize, the_value: T) {
        self.data[the_index - self.lower] = the_value;
    }

    /// OCCT `ChangeValue(theIndex)` — mutable element access.
    pub fn change_value(&mut self, the_index: usize) -> &mut T {
        &mut self.data[the_index - self.lower]
    }

    /// OCCT `First()`.
    pub fn first(&self) -> &T {
        self.data.first().unwrap()
    }

    /// OCCT `Last()`.
    pub fn last(&self) -> &T {
        self.data.last().unwrap()
    }

    /// OCCT `Resize(theLower, theUpper, theToCopy = false)` — old values are
    /// dropped and every element default-constructed (OCCT); `the_fill`
    /// supplies the default.  Every rcad call site either refills the whole
    /// array afterwards or reads none of the new slots before writing them.
    pub fn resize(&mut self, the_lower: usize, the_upper: usize, the_fill: T) {
        let new_len = if the_upper + 1 >= the_lower {
            the_upper + 1 - the_lower
        } else {
            0
        };
        self.lower = the_lower;
        self.data.clear();
        self.data.resize(new_len, the_fill);
    }
}
