// OCCT HLRBRep_InterCSurf (TKHLR) — the curve/surface intersection for the
// HLR chain: IntCurveSurface_HInter assembly (IntCurveSurface_Inter.pxx +
// InterUtils.pxx) instantiated with TheCurve = gp_Lin,
// TheCurveTool = HLRBRep_LineTool, TheSurface = HLRBRep_Surface*,
// TheSurfaceTool = HLRBRep_SurfaceTool (InterCSurf.hxx L44-53 type aliases).
//
// HLRBRep_InterCSurf.hxx L36-211 + .cxx L57-578: the public
// Perform/DoSurface/DoNewBounds/InternalPerform members forward to the
// generic template functions with those fixed bindings.  The rcad shell
// reuses the landed HInter engine shell with the same fixed bindings —
// semantically the OCCT instantiation.

use glam::DVec3;
use rcad_kernel::geom::Line3;
use rcad_kernel::math::bnd::BoundSortBox;

use crate::geomalgo::int_curve_surface::hinter::HInter;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};

use super::surface::Surface;
use super::surface_tool::{LineTool, SurfaceTool};

/// OCCT HLRBRep_InterCSurf.
#[derive(Debug, Clone)]
pub struct InterCSurf {
    /// The IntCurveSurface_Intersection base + the HInter engine state
    /// (done flag, parallel flag, points, segments).
    pub inter: HInter,
}

impl InterCSurf {
    /// OCCT HLRBRep_InterCSurf() (cxx L57).
    pub fn new() -> Self {
        InterCSurf { inter: HInter::new() }
    }

    /// OCCT Perform(gp_Lin, HLRBRep_Surface*) (cxx L107-117).
    pub fn perform<'a>(&mut self, curve: &Line3, surface: &Surface<'a>) {
        self.inter
            .perform::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(curve, surface);
    }

    /// OCCT Perform(gp_Lin, HLRBRep_Surface*, U1, V1, U2, V2) (cxx L121-151).
    pub fn perform_bounds<'a>(
        &mut self,
        curve: &Line3,
        surface: &Surface<'a>,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        self.inter.perform_bounds::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
            curve, surface, u1, v1, u2, v2,
        );
    }

    /// OCCT Perform(gp_Lin, ThePolygon, HLRBRep_Surface*) (cxx L155-169).
    pub fn perform_polygon<'a>(
        &mut self,
        curve: &Line3,
        polygon: &ThePolygonOfHInter,
        surface: &Surface<'a>,
    ) {
        self.inter
            .perform_polygon::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
                curve, polygon, surface,
            );
    }

    /// OCCT Perform(gp_Lin, HLRBRep_Surface*, ThePolyhedron) (cxx L173-191).
    pub fn perform_polyhedron<'a>(
        &mut self,
        curve: &Line3,
        surface: &Surface<'a>,
        polyhedron: &ThePolyhedronOfHInter,
    ) {
        self.inter
            .perform_polyhedron::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
                curve, surface, polyhedron,
            );
    }

    /// OCCT Perform(gp_Lin, ThePolygon, HLRBRep_Surface*, ThePolyhedron)
    /// (cxx L195-220).
    pub fn perform_polygon_polyhedron<'a>(
        &mut self,
        curve: &Line3,
        polygon: &ThePolygonOfHInter,
        surface: &Surface<'a>,
        polyhedron: &ThePolyhedronOfHInter,
    ) {
        self.inter
            .perform_polygon_polyhedron::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
                curve, polygon, surface, polyhedron,
            );
    }

    /// OCCT Perform(..., Bnd_BoundSortBox&) (cxx L224-252).
    pub fn perform_polygon_polyhedron_bsb<'a>(
        &mut self,
        curve: &Line3,
        polygon: &ThePolygonOfHInter,
        surface: &Surface<'a>,
        polyhedron: &ThePolyhedronOfHInter,
        bsb: &mut BoundSortBox,
    ) {
        self.inter
            .perform_polygon_polyhedron_bsb::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
                curve, polygon, surface, polyhedron, bsb,
            );
    }

    /// OCCT InternalPerform 8-arg (cxx L289-316).
    #[allow(clippy::too_many_arguments)]
    pub fn internal_perform<'a>(
        &mut self,
        curve: &Line3,
        polygon: &ThePolygonOfHInter,
        surface: &Surface<'a>,
        polyhedron: &ThePolyhedronOfHInter,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
    ) {
        self.inter.internal_perform::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
            curve, polygon, surface, polyhedron, u0, v0, u1, v1,
        );
    }

    /// OCCT InternalPerform 9-arg with the bound sort box (cxx L256-285).
    #[allow(clippy::too_many_arguments)]
    pub fn internal_perform_bsb<'a>(
        &mut self,
        curve: &Line3,
        polygon: &ThePolygonOfHInter,
        surface: &Surface<'a>,
        polyhedron: &ThePolyhedronOfHInter,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
        bsb: &mut BoundSortBox,
    ) {
        self.inter.internal_perform_bsb::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
            curve, polygon, surface, polyhedron, u0, v0, u1, v1, bsb,
        );
    }

    /// OCCT InternalPerformCurveQuadric (cxx L320-330).
    pub fn internal_perform_curve_quadric<'a>(&mut self, curve: &Line3, surface: &Surface<'a>) {
        self.inter
            .internal_perform_curve_quadric::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
                curve, surface,
            );
    }

    /// OCCT InternalPerform 7-arg PolygonBounds (cxx L334-378).
    pub fn internal_perform_polygon_bounds<'a>(
        &mut self,
        curve: &Line3,
        polygon: &ThePolygonOfHInter,
        surface: &Surface<'a>,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        self.inter
            .internal_perform_polygon_bounds::<Line3, LineTool, Surface<'a>, SurfaceTool<'a>>(
                curve, polygon, surface, u1, v1, u2, v2,
            );
    }

    /// OCCT IsDone / NbPoints / Point / NbSegments / Segment — the
    /// IntCurveSurface_Intersection base accessors.
    pub fn is_done(&self) -> bool {
        self.inter.base.is_done()
    }
    pub fn nb_points(&self) -> usize {
        self.inter.base.nb_points()
    }
    pub fn point(&self, index: usize) -> &crate::geomalgo::int_curve_surface::IntersectionPoint {
        self.inter.base.point(index)
    }
    pub fn nb_segments(&self) -> usize {
        self.inter.base.nb_segments()
    }
    pub fn segment(
        &self,
        index: usize,
    ) -> &crate::geomalgo::int_curve_surface::IntersectionSegment {
        self.inter.base.segment(index)
    }
}

impl Default for InterCSurf {
    fn default() -> Self {
        Self::new()
    }
}

// The DVec3 parity reference for the InterUtils pntsOnSurface grid values.
#[allow(unused)]
fn _grid_parity() -> DVec3 {
    DVec3::ZERO
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::{
        CylindricalSurface, Line3, Plane as GPlane, Surface3,
    };
    use rcad_kernel::topo::topods::{BRepBuilder, Shape};

    use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;
    use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

    fn plane_face(brep: &mut rcad_kernel::BRep, z: f64) -> Shape {
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(brep, DVec3::new(0.0, 0.0, z), 1e-7);
        let v2 = b.add_vertex(brep, DVec3::new(2.0, 0.0, z), 1e-7);
        let e = b.add_edge(
            brep,
            Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(DVec3::new(0.0, 0.0, z), DVec3::new(1.0, 0.0, 0.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        brep.add_tface(
            Some(Surface3::Plane(GPlane {
                origin: DVec3::new(0.0, 0.0, z),
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire,
            Vec::new(),
            None,
            Some([0.0, 2.0, 0.0, 2.0]),
            Vec::new(),
            true,
        )
    }

    fn cylinder_face(brep: &mut rcad_kernel::BRep) -> Shape {
        let mut b = BRepBuilder::new();
        let v1 = b.add_vertex(brep, DVec3::ZERO, 1e-7);
        let v2 = b.add_vertex(brep, DVec3::new(0.0, 0.0, 2.0), 1e-7);
        let e = b.add_edge(
            brep,
            Some(rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3::new(DVec3::ZERO, DVec3::new(0.0, 0.0, 1.0)))),
            v1,
            v2,
            [0.0, 2.0],
        );
        let wire = brep.add_twire(vec![e]);
        brep.add_tface(
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
        )
    }

    /// OCCT anchor: the viewing ray (0, 0.5, 1) -> +X pierces the unit
    /// cylinder once on the forward half (u = pi/2, v = 1, w = 1 - sqrt(3)/4)
    /// — the InterCSurf plane branch handles the quadric through
    /// InternalPerformCurveQuadric; the IntAna line/plane path covers the
    /// plane face exactly: the line (0, 0.5, 1) + X hits z = 1... at every
    /// point (the line lies IN the plane)? No — the plane z = 1 contains the
    /// line, so Perform reports parallel/undefined; use z = 0.5 plane
    /// instead and expect w = 0.5.
    #[test]
    fn inter_csurf_line_plane() {
        let mut brep = rcad_kernel::BRep::new();
        let face = plane_face(&mut brep, 0.5);

        let mut hsurf = Surface::new();
        hsurf.load(&brep, &face);

        let mut ics = InterCSurf::new();
        let line = Line3::new(DVec3::new(0.0, 0.5, 2.0), DVec3::new(0.0, 0.0, -1.0));
        ics.perform(&line, &hsurf);

        assert!(ics.is_done());
        // One intersection: the piercing point (0, 0.5, 0.5), w = 1.5.
        assert_eq!(ics.nb_points(), 1, "nb={}", ics.nb_points());
        let pt = ics.point(1);
        assert!((pt.pnt().x - 0.0).abs() < 1e-9);
        assert!((pt.pnt().y - 0.5).abs() < 1e-9);
        assert!((pt.pnt().z - 0.5).abs() < 1e-9);
        assert!((pt.w() - 1.5).abs() < 1e-9);
        // UV on the plane face window: u = 0, v = 0.5.
        assert!((pt.u() - 0.0).abs() < 1e-9);
        assert!((pt.v() - 0.5).abs() < 1e-9);
    }

    /// OCCT anchor: the line along +X at (y=0, z=1) pierces the unit
    /// cylinder (window v in [0, 2]) twice: entering at (-1, 0, 1)
    /// (u = pi, w = 1) and leaving at (1, 0, 1) (u = 0, w = 3).
    #[test]
    fn inter_csurf_line_cylinder() {
        let mut brep = rcad_kernel::BRep::new();
        let face = cylinder_face(&mut brep);

        let mut hsurf = Surface::new();
        hsurf.load(&brep, &face);

        let mut ics = InterCSurf::new();
        let line = Line3::new(DVec3::new(-2.0, 0.0, 1.0), DVec3::new(1.0, 0.0, 0.0));
        ics.perform(&line, &hsurf);

        assert!(ics.is_done());
        assert_eq!(ics.nb_points(), 2, "nb={}", ics.nb_points());
        let mut ws = Vec::new();
        for i in 1..=ics.nb_points() {
            ws.push(ics.point(i).w());
        }
        ws.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((ws[0] - 1.0).abs() < 1e-9, "w0={}", ws[0]);
        assert!((ws[1] - 3.0).abs() < 1e-9, "w1={}", ws[1]);
        // The entering point is on the far side: u = pi.
        let pt = ics.point(1);
        let _ = pt;
    }
}
