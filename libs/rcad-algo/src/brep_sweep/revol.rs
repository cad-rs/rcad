//! OCCT BRepSweep_Revol (TKPrim/BRepSweep) — natural constructors to build
//! BRepSweep rotated swept Primitives.
//!
//! Sources:
//! - BRepSweep_Revol.hxx L33-95
//! - BRepSweep_Revol.cxx L29-163

use rcad_kernel::precision::ANGULAR;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};

use super::rotation::GpAx1;
use super::rotation::BRepSweepRotation;
use super::num_linear_regular_sweep::NumLinearRegularSweepSlots;
use super::sweep_num_shape::SweepNumShape;

/// OCCT BRepSweep_Revol (BRepSweep_Revol.hxx L33-95).
pub struct BRepSweepRevol {
    /// OCCT: myRotation (BRepSweep_Rotation).
    pub my_rotation: BRepSweepRotation,
}

impl BRepSweepRevol {
    /// OCCT BRepSweep_Revol::BRepSweep_Revol(S, A, D, C) (cxx L29-37) —
    /// builds the Revol of meridian S axis A and angle D.  If C is true S is
    /// copied.
    pub fn with_angle(s: &Shape, a: GpAx1, d: f64, c: bool) -> Self {
        // OCCT L33: myRotation(S.Oriented(TopAbs_FORWARD), NumShape(D),
        // Location(Ax, D), Axe(Ax, D), Angle(D), C).
        let mut oriented = s.clone();
        oriented.orientation = Orientation::Forward;
        let my_rotation = BRepSweepRotation::new(
            &oriented,
            &Self::num_shape(d),
            Self::location(a, d),
            Self::axe_ax1(a, d),
            Self::angle_of(d),
            c,
        );
        // OCCT L35-36: Standard_ConstructionError_Raise_if(Angle(D) <=
        // Precision::Angular(), "BRepSweep_Revol::Constructor").
        if Self::angle_of(d) <= ANGULAR {
            panic!("Standard_ConstructionError: BRepSweep_Revol::Constructor");
        }
        BRepSweepRevol { my_rotation }
    }

    /// OCCT BRepSweep_Revol::BRepSweep_Revol(S, A, C) (cxx L41-50) — builds
    /// the Revol of meridian S axis A and angle 2*Pi.
    pub fn with_full_circle(s: &Shape, a: GpAx1, c: bool) -> Self {
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut oriented = s.clone();
        oriented.orientation = Orientation::Forward;
        let my_rotation = BRepSweepRotation::new(
            &oriented,
            &Self::num_shape(two_pi),
            Self::location(a, two_pi),
            Self::axe_ax1(a, two_pi),
            Self::angle_of(two_pi),
            c,
        );
        BRepSweepRevol { my_rotation }
    }

    /// OCCT BRepSweep_Revol::Shape() (cxx L54-57) — the TopoDS Shape attached
    /// to the Revol.
    pub fn shape(&mut self) -> Shape {
        self.my_rotation.shape_full()
    }

    /// OCCT BRepSweep_Revol::Shape(aGenS) (cxx L61-64).
    pub fn shape_of(&mut self, a_gen_s: &Shape) -> Shape {
        self.my_rotation.shape_of_gen(a_gen_s)
    }

    /// OCCT BRepSweep_Revol::FirstShape() (cxx L68-71) — the first shape of
    /// the revol (coinciding with the generating shape).
    pub fn first_shape(&mut self) -> Shape {
        self.my_rotation.first_shape_full()
    }

    /// OCCT BRepSweep_Revol::FirstShape(aGenS) (cxx L75-78).
    pub fn first_shape_of(&mut self, a_gen_s: &Shape) -> Shape {
        self.my_rotation.first_shape(a_gen_s)
    }

    /// OCCT BRepSweep_Revol::LastShape() (cxx L82-85).
    pub fn last_shape(&mut self) -> Shape {
        self.my_rotation.last_shape_full()
    }

    /// OCCT BRepSweep_Revol::LastShape(aGenS) (cxx L89-92).
    pub fn last_shape_of(&mut self, a_gen_s: &Shape) -> Shape {
        self.my_rotation.last_shape(a_gen_s)
    }

    /// OCCT BRepSweep_Revol::Axe() (cxx L153-156) — returns the axis.
    pub fn axe(&self) -> GpAx1 {
        self.my_rotation.axe()
    }

    /// OCCT BRepSweep_Revol::Angle() (cxx L146-149) — returns the angle.
    pub fn angle(&self) -> f64 {
        self.my_rotation.angle()
    }

    /// OCCT BRepSweep_Revol::IsUsed(aGenS) (cxx L160-163).
    pub fn is_used(&self, a_gen_s: &Shape) -> bool {
        self.my_rotation.is_used(a_gen_s)
    }

    /// OCCT BRepSweep_Revol::NumShape(D) (cxx L96-108) — builds the NumShape
    /// (private).
    fn num_shape(d: f64) -> SweepNumShape {
        let mut n = SweepNumShape::new();
        if (Self::angle_of(d) - 2.0 * std::f64::consts::PI).abs() <= ANGULAR {
            n.init(2, ShapeType::Edge, true, false, false);
        } else {
            n.init_indexed(2, ShapeType::Edge);
        }
        n
    }

    /// OCCT BRepSweep_Revol::Location(Ax, D) (cxx L112-118) — builds the
    /// Location (private).  Architecture difference: the TopLoc_Location
    /// travels as a DAffine3 (the rotation transform).
    fn location(ax: GpAx1, d: f64) -> glam::DAffine3 {
        let a = Self::axe_ax1(ax, d);
        // OCCT gp_Trsf::SetRotation(Ax1, Ang): the rotation ABOUT the axis
        // line (the location point included).
        glam::DAffine3::from_translation(a.0)
            * glam::DAffine3::from_axis_angle(a.1.normalize_or_zero(), Self::angle_of(d))
            * glam::DAffine3::from_translation(-a.0)
    }

    /// OCCT BRepSweep_Revol::Axe(Ax, D) (cxx L122-130) — builds the axis
    /// (private): reversed for a negative angle.
    fn axe_ax1(ax: GpAx1, d: f64) -> GpAx1 {
        let mut a = ax;
        if d < 0.0 {
            // OCCT: A.Reverse() — the gp_Ax1 reversal (location unchanged,
            // direction negated).
            a.1 = -a.1;
        }
        a
    }

    /// OCCT BRepSweep_Revol::Angle(D) (cxx L134-142) — computes the angle
    /// (private).
    fn angle_of(d: f64) -> f64 {
        let mut d = d.abs();
        while d > 2.0 * std::f64::consts::PI + ANGULAR {
            d -= 2.0 * std::f64::consts::PI;
        }
        d
    }
}
