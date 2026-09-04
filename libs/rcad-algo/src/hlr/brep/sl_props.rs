// OCCT HLRBRep_SLProps + HLRBRep_SLPropsATool (TKHLR) — the surface local
// properties engine instantiated over HLRBRep_Surface.
//
// HLRBRep_SLProps.hxx L20-24:
//   using HLRBRep_SLProps = GeomLProp_SLPropsBase<HLRBRep_SurfacePtr,
//       LProp_SurfaceUtils::ToolAccess<HLRBRep_SLPropsATool>>;
//
// HLRBRep_SLPropsATool.hxx (the tool statics delegate to the
// HLRBRep_Surface adaptor methods; Bounds reads the First/Last U/V
// parameters).  The rcad tool is the `SLPropsSurface` impl for [`Surface`].

use rcad_kernel::base::geom_lprop::{SLPropsSurface, SlPropsBase};

use super::surface::Surface;

/// OCCT HLRBRep_SLProps — the surface local-properties engine over the
/// projected face adaptor.
pub type SLProps<'a> = SlPropsBase<'a, Surface<'a>>;

/// OCCT HLRBRep_SLPropsATool (the ToolAccess policy for HLRBRep_Surface).
impl SLPropsSurface for Surface<'_> {
    /// OCCT Tool::Value(A, U, V, P) — the surface point.
    fn eval_d0(&self, u: f64, v: f64) -> DVec3 {
        self.value(u, v)
    }

    /// OCCT Tool::D1(A, U, V, P, D1u, D1v).
    fn eval_d1(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3) {
        self.d1(u, v)
    }

    /// OCCT Tool::D2(A, U, V, P, D1u, D1v, D2u, D2v, Duv).
    fn eval_d2(&self, u: f64, v: f64) -> (DVec3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        self.d2(u, v)
    }

    /// OCCT Tool::Bounds(A, U1, V1, U2, V2) — the First/Last U/V parameters.
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (
            self.first_u_parameter(),
            self.first_v_parameter(),
            self.last_u_parameter(),
            self.last_v_parameter(),
        )
    }
}

use rcad_kernel::geom::Vec3;
use glam::DVec3;
#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{CylindricalSurface, Surface3};

    use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;
    use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

    /// OCCT anchor: SLProps over an HLRBRep_Surface wrapping a unit
    /// cylinder — the normal is radial, the principal curvatures are
    /// -1 (circular) and 0 (axial) with the outward radial normal.
    #[test]
    fn hlr_slprops_cylinder() {
        let mut brep = rcad_kernel::BRep::new();
        let mut b = rcad_kernel::topo::topods::BRepBuilder::new();
        let v1 = b.add_vertex(&mut brep, DVec3::ZERO, 1e-7);
        let v2 = b.add_vertex(&mut brep, DVec3::new(0.0, 0.0, 2.0), 1e-7);
        let e = b.add_edge(
            &mut brep,
            Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::new(0.0, 0.0, 1.0),
            })),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        let face = brep.add_tface(
            Some(Surface3::Cylinder(CylindricalSurface {
                origin: DVec3::ZERO,
                axis: DVec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
                ref_dir: DVec3::new(1.0, 0.0, 0.0),
                y_dir: None,
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, std::f64::consts::TAU, 0.0, 2.0]),
            Vec::new(),
            true,
        );

        // The HLRBRep_Surface face adaptor (the SLProps curve type is the
        // OCCT HLRBRep_SurfacePtr).
        let mut hsurf = Surface::new();
        hsurf.load(&brep, &face);

        let mut slp = SLProps::new(2, 1e-12);
        slp.set_surface(&hsurf);
        slp.set_parameters(0.0, 1.0);

        // The point (1, 0, 1); the radial normal.
        assert!((slp.value().x - 1.0).abs() < 1e-12);
        let n = slp.normal().expect("normal defined");
        assert!((n.x - 1.0).abs() < 1e-12 && n.y.abs() < 1e-12 && n.z.abs() < 1e-12);

        // Principal curvatures: 0 (axial) and -1 (circular, inward D2u);
        // the kernel trait-default D2 is a finite difference (h = 1e-5),
        // so the tolerance is 1e-6.
        assert!(slp.max_curvature().abs() < 1e-6, "max={}", slp.max_curvature());
        assert!((slp.min_curvature() + 1.0).abs() < 1e-6);
        assert!((slp.mean_curvature() + 0.5).abs() < 1e-6);
        assert!(slp.gaussian_curvature().abs() < 1e-6);
    }
}
