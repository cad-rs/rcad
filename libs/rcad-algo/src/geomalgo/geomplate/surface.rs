//! OCCT GeomPlate_Surface (TKGeomAlgo/GeomPlate) — 1:1 port of the
//! anchor-relevant subset of GeomPlate_Surface.cxx: ctor (L41-50),
//! Bounds/SetBounds (L97-112, L278-293), RealBounds (L296-301),
//! Constraints (L303-309), EvalD0 (L197-203).
//!
//! OCCT wraps `handle<Geom_Surface>` for the basis; rcad stores the
//! double-precision `Surface3` instead (architecture mapping).  The
//! Geom_Surface-inherited API not exercised by the plate pipeline
//! (UReverse/VReverse/TransformParameters/D1../ParametricTransformation)
//! is left unported.  "IsNull()" at the caller side maps to Option<...>.

use rcad_kernel::geom::{Surface3, SurfaceEval};

use super::super::plate::Plate;
use glam::{DVec2, DVec3};

/// OCCT GeomPlate_Surface — a surface defined by an initial basis surface
/// plus a Plate deformation computed in its parameter space.
#[derive(Clone)]
pub struct GeomPlateSurface {
    my_surfinter: Plate,
    my_surfinit: Surface3,
    my_umin: f64,
    my_umax: f64,
    my_vmin: f64,
    my_vmax: f64,
}

impl GeomPlateSurface {
    /// OCCT ctor (GeomPlate_Surface.cxx L41-50).
    pub fn new(surfinit: Surface3, surfinter: Plate) -> Self {
        GeomPlateSurface {
            my_surfinter: surfinter,
            my_surfinit: surfinit,
            my_umin: 0.0,
            my_umax: 0.0,
            my_vmin: 0.0,
            my_vmax: 0.0,
        }
    }

    /// OCCT Bounds (GeomPlate_Surface.cxx L97-111).
    pub fn bounds(&self) -> (f64, f64, f64, f64) {
        (
            self.my_umin,
            self.my_umax,
            self.my_vmin,
            self.my_vmax,
        )
    }

    /// OCCT SetBounds (GeomPlate_Surface.cxx L278-293).
    pub fn set_bounds(&mut self, u_min: f64, u_max: f64, v_min: f64, v_max: f64) {
        if u_min > u_max || v_min > v_max {
            panic!("Bounds haven't the good sense");
        }
        if u_min == u_max || v_min == v_max {
            panic!("Bounds are equal");
        }
        self.my_umin = u_min;
        self.my_umax = u_max;
        self.my_vmin = v_min;
        self.my_vmax = v_max;
    }

    /// OCCT RealBounds (GeomPlate_Surface.cxx L296-301).
    pub fn real_bounds(&self) -> (f64, f64, f64, f64) {
        let (mut u1, mut u2, mut v1, mut v2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        self.my_surfinter.uv_box(&mut u1, &mut u2, &mut v1, &mut v2);
        (u1, u2, v1, v2)
    }

    /// OCCT Constraints (GeomPlate_Surface.cxx L303-309).
    pub fn constraints(&self, seq: &mut Vec<DVec2>) {
        self.my_surfinter.uv_constraints(seq);
    }

    /// OCCT EvalD0 (GeomPlate_Surface.cxx L197-203): basis surface point plus
    /// the plate deformation at (U, V).
    pub fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        let a_surf_p = self.my_surfinit.point_at(u, v);
        let p3 = self.my_surfinter.evaluate(DVec2::new(u, v));
        p3 + a_surf_p
    }

    /// OCCT basis surface access (hxx `Surface()`).
    pub fn basis_surface(&self) -> &Surface3 {
        &self.my_surfinit
    }

    /// OCCT mySurfinter access (hxx `Plate()`).
    pub fn plate(&self) -> &Plate {
        &self.my_surfinter
    }
}
