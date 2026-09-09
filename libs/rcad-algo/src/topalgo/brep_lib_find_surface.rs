//! OCCT BRepLib_FindSurface (TKTopAlgo/BRepLib — BRepLib_FindSurface.cxx,
//! 621 lines) — the planar-surface finder over a shape.
//!
//! GAP carrier (plan §0.6): the full body (BRepLib_FindSurface::Init —
//! the face/wire exploration, the BRepAdaptor_Surface least-squares plane
//! fit and the GeomLib::FindPlane checks) is out of this dispatch; the
//! OCCT failure path is preserved (Init raises / Found() stays false from
//! construction).  Consumers: BRepFill_Sweep::BuildFace (part B static,
//! cxx L584) and BRepFill_Axe.

use rcad_kernel::geom::Surface3;
use rcad_kernel::topo::topods::Shape;

/// OCCT BRepLib_FindSurface (BRepLib_FindSurface.hxx L45-100).
pub struct BRepLibFindSurface {
    /// OCCT myFound.
    my_found: bool,
    /// OCCT mySurface.
    my_surface: Surface3,
    /// OCCT myTol.
    my_tol: f64,
}

impl BRepLibFindSurface {
    /// OCCT BRepLib_FindSurface(S, Tol = -1, OnlyPlane = false)
    /// (BRepLib_FindSurface.cxx L44-52) — the deferred Init form.
    pub fn new(_brep: &mut rcad_kernel::topo::topods::BRep, s: &Shape, tol: f64, only_plane: bool) -> Self {
        let mut this = BRepLibFindSurface {
            my_found: false,
            my_surface: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: glam::DVec3::ZERO,
                normal: glam::DVec3::Z,
                u_dir: glam::DVec3::X,
                v_dir: glam::DVec3::Y,
            }),
            my_tol: if tol < 0.0 { 0.0 } else { tol },
        };
        this.init(_brep, s, tol, only_plane);
        this
    }

    /// OCCT Init(S, Tol, OnlyPlane) (BRepLib_FindSurface.cxx L58-560).
    pub fn init(
        &mut self,
        _brep: &mut rcad_kernel::topo::topods::BRep,
        _s: &Shape,
        _tol: f64,
        _only_plane: bool,
    ) {
        panic!(
            "GAP: BRepLib_FindSurface::Init (TKTopAlgo/BRepLib, \
             BRepLib_FindSurface.cxx L58-560) is not translated — \
             see file header (plan section 0.6)"
        );
    }

    /// OCCT Found() (cxx L565-568).
    pub fn found(&self) -> bool {
        self.my_found
    }

    /// OCCT Surface() (cxx L572-575).
    pub fn surface(&self) -> Surface3 {
        self.my_surface.clone()
    }

    /// OCCT Tolerance() (cxx L581-584).
    pub fn tolerance(&self) -> f64 {
        self.my_tol
    }
}
