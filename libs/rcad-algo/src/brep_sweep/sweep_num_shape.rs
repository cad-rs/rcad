//! OCCT Sweep_NumShape (TKPrim/Sweep) — the simple indexed representation of
//! a Directing Edge topology.
//!
//! Sources:
//! - Sweep_NumShape.hxx L30-95
//! - Sweep_NumShape.cxx L21-80
//! - Sweep_NumShape.lxx L19-50 (the inline accessors)
//!
//! Architecture differences:
//! - `TopAbs_ShapeEnum` -> `rcad_kernel::topods::ShapeType`.
//! - `TopAbs_Orientation` -> `rcad_kernel::topods::Orientation`.

use rcad_kernel::topods::{Orientation, ShapeType};

/// OCCT Sweep_NumShape (Sweep_NumShape.hxx L30-91).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SweepNumShape {
    my_type: ShapeType,
    my_index: i32,
    my_closed: bool,
    my_beg_inf: bool,
    my_end_inf: bool,
}

impl Default for SweepNumShape {
    /// OCCT Sweep_NumShape::Sweep_NumShape() (cxx L21-29).
    fn default() -> Self {
        Self::new()
    }
}

impl SweepNumShape {
    /// OCCT Sweep_NumShape::Sweep_NumShape() (cxx L21-29) — creates a dummy
    /// indexed edge.
    pub fn new() -> Self {
        SweepNumShape {
            my_type: ShapeType::Shape,
            my_index: 0,
            my_closed: false,
            my_beg_inf: false,
            my_end_inf: false,
        }
    }

    /// OCCT Sweep_NumShape::Sweep_NumShape(Index, Type, Closed, BegInf, EndInf)
    /// (cxx L33-44) — creates a new simple indexed edge.  OCCT has default
    /// arguments (Closed = BegInf = EndInf = false); Rust has no overloading:
    /// use `with_defaults` for the OCCT default-argument form.
    pub fn with_all(
        index: i32,
        the_type: ShapeType,
        closed: bool,
        beg_inf: bool,
        end_inf: bool,
    ) -> Self {
        SweepNumShape {
            my_type: the_type,
            my_index: index,
            my_closed: closed,
            my_beg_inf: beg_inf,
            my_end_inf: end_inf,
        }
    }

    /// OCCT Sweep_NumShape::Sweep_NumShape(Index, Type) with the OCCT default
    /// arguments Closed = BegInf = EndInf = false (hxx L49-53).
    pub fn new_indexed(index: i32, the_type: ShapeType) -> Self {
        Self::with_all(index, the_type, false, false, false)
    }

    /// OCCT Sweep_NumShape::Init(Index, Type, Closed, BegInf, EndInf)
    /// (cxx L48-59) — reinitializes a simple indexed edge.
    pub fn init(&mut self, index: i32, the_type: ShapeType, closed: bool, beg_inf: bool, end_inf: bool) {
        self.my_index = index;
        self.my_type = the_type;
        self.my_closed = closed;
        self.my_beg_inf = beg_inf;
        self.my_end_inf = end_inf;
    }

    /// OCCT Sweep_NumShape::Init(Index, Type) with the OCCT default arguments
    /// (hxx L67-71).
    pub fn init_indexed(&mut self, index: i32, the_type: ShapeType) {
        self.init(index, the_type, false, false, false);
    }

    /// OCCT Sweep_NumShape::Index() (lxx L19-22).
    pub fn index(&self) -> i32 {
        self.my_index
    }

    /// OCCT Sweep_NumShape::Type() (lxx L26-29).
    pub fn type_(&self) -> ShapeType {
        self.my_type
    }

    /// OCCT Sweep_NumShape::Closed() (lxx L33-36).
    pub fn closed(&self) -> bool {
        self.my_closed
    }

    /// OCCT Sweep_NumShape::BegInfinite() (lxx L40-43).
    pub fn beg_infinite(&self) -> bool {
        self.my_beg_inf
    }

    /// OCCT Sweep_NumShape::EndInfinite() (lxx L47-50).
    pub fn end_infinite(&self) -> bool {
        self.my_end_inf
    }

    /// OCCT Sweep_NumShape::Orientation() (cxx L63-80).
    pub fn orientation(&self) -> Orientation {
        if self.my_type == ShapeType::Edge {
            Orientation::Forward
        } else if self.my_index == 2 {
            Orientation::Forward
        } else {
            Orientation::Reversed
        }
    }
}
