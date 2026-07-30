// OCCT BRepPrimAPI_MakeBox 1:1 translation.
//
// OCCT: ModelingAlgorithms/TKPrim/BRepPrimAPI/BRepPrimAPI_MakeBox.hxx/.cxx
//
// Constructs a parallelepiped box. Internally delegates to a box builder
// that creates the TopoDS_Shape-compatible BRep topology.

use rcad_kernel::{topods, BRep};
use glam::DVec3;

/// OCCT BRepPrimAPI_MakeBox — box primitive builder.
///
/// OCCT L50: class BRepPrimAPI_MakeBox : public BRepBuilderAPI_MakeShape
pub struct MakeBox {
    // OCCT: BRepPrim_Wedge myWedge;
    pmin: DVec3,
    dx: f64,
    dy: f64,
    dz: f64,
}

impl MakeBox {
    /// OCCT L56: default constructor
    pub fn new() -> Self {
        MakeBox { pmin: DVec3::ZERO, dx: 0.0, dy: 0.0, dz: 0.0 }
    }

    /// OCCT L59: Make a box with a corner at 0,0,0 and the other dx,dy,dz
    /// OCCT: BRepPrimAPI_MakeBox(double dx, double dy, double dz)
    pub fn new_at_origin(dx: f64, dy: f64, dz: f64) -> Self {
        // OCCT L46-51: myWedge(gp_Ax2(pmin(P,dx,dy,dz), Z, X), abs(dx), abs(dy), abs(dz))
        // pmin adjusts for negative sizes
        MakeBox {
            pmin: DVec3::new(if dx < 0.0 { dx } else { 0.0 }, if dy < 0.0 { dy } else { 0.0 }, if dz < 0.0 { dz } else { 0.0 }),
            dx: dx.abs(), dy: dy.abs(), dz: dz.abs(),
        }
    }

    /// OCCT L62-65: Make a box with a corner at P and size dx, dy, dz
    /// OCCT: BRepPrimAPI_MakeBox(const gp_Pnt& P, double dx, double dy, double dz)
    pub fn new_at(p: DVec3, dx: f64, dy: f64, dz: f64) -> Self {
        // OCCT L55-64: myWedge(gp_Ax2(pmin(P,dx,dy,dz), Z, X), abs(dx), abs(dy), abs(dz))
        MakeBox {
            pmin: DVec3::new(
                p.x + if dx < 0.0 { dx } else { 0.0 },
                p.y + if dy < 0.0 { dy } else { 0.0 },
                p.z + if dz < 0.0 { dz } else { 0.0 },
            ),
            dx: dx.abs(), dy: dy.abs(), dz: dz.abs(),
        }
    }

    /// OCCT L68: Make a box with corners P1, P2
    /// OCCT: BRepPrimAPI_MakeBox(const gp_Pnt& P1, const gp_Pnt& P2)
    pub fn new_between(p1: DVec3, p2: DVec3) -> Self {
        // OCCT L73-79: myWedge(gp_Ax2(pmin(P1,P2), Z, X), |dx|, |dy|, |dz|)
        MakeBox {
            pmin: DVec3::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z)),
            dx: (p2.x - p1.x).abs(),
            dy: (p2.y - p1.y).abs(),
            dz: (p2.z - p1.z).abs(),
        }
    }

    /// Build the box — OCCT L98-99: Build(const Message_ProgressRange&)
    pub fn build(&self) -> Result<BRep, crate::BuildError> {
        // Delegates to the existing box builder.
        // OCCT: myWedge.Shell() builds the actual topology via BRepPrim_Builder.
        // rcad: crate::builder::make_box_brep builds the BRep directly.
        crate::builder::make_box_brep(
            self.pmin, DVec3::X, DVec3::Y,
            self.dx, self.dy, self.dz,
        )
    }

    /// Returns the constructed box as a BRep — OCCT L101-103: Shell()
    pub fn shell(&self) -> Result<BRep, crate::BuildError> {
        self.build()
    }

    /// Returns the constructed box as a BRep — OCCT L106-107: Solid()
    pub fn solid(&self) -> Result<BRep, crate::BuildError> {
        self.build()
    }
}

/// Free function: build a box with corner at origin, size dx×dy×dz.
/// Equivalent to OCCT: BRepPrimAPI_MakeBox(dx, dy, dz).Solid()
pub fn box_brep(dx: f64, dy: f64, dz: f64) -> Result<BRep, crate::BuildError> {
    MakeBox::new_at_origin(dx, dy, dz).build()
}

/// Legacy: box_brep with origin, axes, dimensions (matches old builder API).
pub fn make_box_brep(
    origin: DVec3, x_dir: DVec3, y_dir: DVec3,
    width: f64, height: f64, depth: f64,
) -> Result<BRep, crate::BuildError> {
    crate::builder::make_box_brep(origin, x_dir, y_dir, width, height, depth)
}
