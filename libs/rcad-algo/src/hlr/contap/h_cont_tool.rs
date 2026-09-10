// OCCT Contap_HContTool (TKHLR) — the sampling/projection tool for the
// contour engine (`TheSOBTool` of Contap_TheSearch, `TheSITool` of
// Contap_TheSearchInside).
//
// Contap_HContTool.hxx L29-140 + Contap_HContTool.cxx L24-371.
// The file-level mutable statics `uinf/vinf/usup/vsup` (cxx L24) map to a
// thread_local (the §5 rule of docs/tkhlr-port-plan.md).
// The Project() body uses Extrema_EPCOfExtPC2d (cxx L269-304).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::Point3;

use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Curve2dType};
use crate::geomalgo::extrema_gen_ext_pc2d::EPCOfExtPC2d;
use crate::geomalgo::int_patch::GeomAbsSurfaceType;

use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};
use crate::topalgo::adaptor3d::hvertex::{HVertex, HVertexBehavior};

thread_local! {
    /// OCCT Contap_HContTool.cxx L24: `static double uinf, vinf, usup, vsup;`
    /// — file-level mutable statics, thread_local per the port rules.
    static SAMPLE_BOUNDS: std::cell::RefCell<[f64; 4]> = const { std::cell::RefCell::new([0.0; 4]) };
}

/// OCCT NbSamplesV (cxx L26-68).
pub fn nb_samples_v(s: &dyn SurfaceAdapter, _v1: f64, _v2: f64) -> i32 {
    let nbs;
    let typ_s = s.get_type();
    match typ_s {
        GeomAbsSurfaceType::Plane => {
            nbs = 2;
        }
        GeomAbsSurfaceType::BezierSurface => {
            nbs = (3 + s.nb_v_poles()) as i32;
        }
        GeomAbsSurfaceType::BSplineSurface => {
            let mut n = (s.nb_v_knots() * s.v_degree()) as i32;
            if n < 2 {
                n = 2;
            }
            nbs = n;
        }
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::Torus
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => {
            nbs = 15;
        }
        _ => {
            nbs = 10;
        }
    }
    nbs
}

/// OCCT NbSamplesU (cxx L70-115).
pub fn nb_samples_u(s: &dyn SurfaceAdapter, _u1: f64, _u2: f64) -> i32 {
    let nbs;
    let typ_s = s.get_type();
    match typ_s {
        GeomAbsSurfaceType::Plane => {
            nbs = 2;
        }
        GeomAbsSurfaceType::BezierSurface => {
            nbs = (3 + s.nb_u_poles()) as i32;
        }
        GeomAbsSurfaceType::BSplineSurface => {
            let mut n = (s.nb_u_knots() * s.u_degree()) as i32;
            if n < 2 {
                n = 2;
            }
            nbs = n;
        }
        GeomAbsSurfaceType::Torus => {
            nbs = 20;
        }
        GeomAbsSurfaceType::Cylinder
        | GeomAbsSurfaceType::Cone
        | GeomAbsSurfaceType::Sphere
        | GeomAbsSurfaceType::SurfaceOfRevolution
        | GeomAbsSurfaceType::SurfaceOfExtrusion => {
            nbs = 10;
        }
        _ => {
            nbs = 10;
        }
    }
    nbs
}

/// OCCT NbSamplePoints (cxx L117-179) — normalizes the sample window and
/// stores it in the file statics.
pub fn nb_sample_points(s: &dyn SurfaceAdapter) -> i32 {
    let mut uinf = s.first_u_parameter();
    let mut usup = s.last_u_parameter();
    let mut vinf = s.first_v_parameter();
    let mut vsup = s.last_v_parameter();

    if usup < uinf {
        let temp = uinf;
        uinf = usup;
        usup = temp;
    }
    if vsup < vinf {
        let temp = vinf;
        vinf = vsup;
        vsup = temp;
    }
    if uinf == f64::MIN && usup == f64::MAX {
        uinf = -1.0e5;
        usup = 1.0e5;
    } else if uinf == f64::MIN {
        uinf = usup - 2.0e5;
    } else if usup == f64::MAX {
        usup = uinf + 2.0e5;
    }

    if vinf == f64::MIN && vsup == f64::MAX {
        vinf = -1.0e5;
        vsup = 1.0e5;
    } else if vinf == f64::MIN {
        vinf = vsup - 2.0e5;
    } else if vsup == f64::MAX {
        vsup = vinf + 2.0e5;
    }
    SAMPLE_BOUNDS.with(|b| *b.borrow_mut() = [uinf, vinf, usup, vsup]);
    if s.get_type() == GeomAbsSurfaceType::BSplineSurface {
        let m = (nb_samples_u(s, uinf, usup) / 3) * (nb_samples_v(s, vinf, vsup) / 3);
        if m > 5 {
            m
        } else {
            5
        }
    } else {
        5
    }
}

/// OCCT SamplePoint (cxx L181-228) — 1-based Index.
pub fn sample_point(s: &dyn SurfaceAdapter, index: i32) -> (f64, f64) {
    let [uinf, vinf, usup, vsup] = SAMPLE_BOUNDS.with(|b| *b.borrow());
    if s.get_type() == GeomAbsSurfaceType::BSplineSurface {
        let nb_int_u = nb_samples_u(s, uinf, usup) / 3;
        let nb_int_v = nb_samples_v(s, vinf, vsup) / 3;
        if nb_int_u * nb_int_v > 5 {
            let ind_u = (index - 1) / nb_int_u; //----   0 --> nbIntV
            let ind_v = (index - 1) - ind_u * nb_int_u; //----   0 --> nbIntU

            let u = uinf + ((usup - uinf) / (nb_int_u as f64 + 1.0)) * (ind_u + 1) as f64;
            let v = vinf + ((vsup - vinf) / (nb_int_v as f64 + 2.0)) * (ind_v + 1) as f64;
            return (u, v);
        }
    }

    match index {
        1 => (0.75 * uinf + 0.25 * usup, 0.75 * vinf + 0.25 * vsup), // 0.25;
        2 => (0.75 * uinf + 0.25 * usup, 0.25 * vinf + 0.75 * vsup), // 0.75;
        3 => (0.25 * uinf + 0.75 * usup, 0.75 * vinf + 0.25 * vsup), // 0.25;
        4 => (0.25 * uinf + 0.75 * usup, 0.25 * vinf + 0.75 * vsup), // 0.75;
        _ => (0.5 * (uinf + usup), 0.5 * (vinf + vsup)),             // 0.5;
    }
}

/// OCCT NbSamplesOnArc (cxx L230-259).
pub fn nb_samples_on_arc(a: &dyn Curve2dAdaptor) -> i32 {
    let curve_type = a.get_type();

    let nbs_on_c = match curve_type {
        Curve2dType::Line => 2.0,
        Curve2dType::Circle | Curve2dType::Ellipse | Curve2dType::Hyperbola | Curve2dType::Parabola => 10.0,
        Curve2dType::BezierCurve => a.nb_poles() as f64,
        Curve2dType::BSplineCurve => (2 + a.nb_knots() * a.degree()) as f64,
        _ => 10.0,
    };
    nbs_on_c as i32
}

/// OCCT Bounds (cxx L261-267).
pub fn bounds(a: &dyn Curve2dAdaptor) -> (f64, f64) {
    (a.first_parameter(), a.last_parameter())
}

/// OCCT Project (cxx L269-304) — projects the 2D point P on the arc C via
/// Extrema_EPCOfExtPC2d(P, C, Nbu=20, epsX=1.0e-8, Tol=1.0e-5); returns
/// (Paramproj, Ptproj) when successful.
pub fn project(c: &dyn Curve2dAdaptor, p: DVec2) -> Option<(f64, DVec2)> {
    let eps_x = 1.0e-8;
    let nbu = 20;
    let tol = 1.0e-5;

    let extrema = EPCOfExtPC2d::with_full_domain(p, c, nbu, eps_x, tol);
    if !extrema.is_done() {
        return None;
    }
    let nbext = extrema.nb_ext();
    if nbext == 0 {
        return None;
    }
    let mut indexmin = 1usize;
    let mut dist2 = extrema.square_distance(1);
    for i in 2..=nbext {
        if extrema.square_distance(i) < dist2 {
            indexmin = i;
            dist2 = extrema.square_distance(i);
        }
    }
    Some((extrema.point(indexmin).parameter(), extrema.point(indexmin).value()))
}

/// OCCT Tolerance (cxx L306-311) — V->Resolution(C).
pub fn tolerance(v: &dyn HVertexBehavior, c: &dyn Curve2dAdaptor) -> f64 {
    v.resolution(c)
}

/// OCCT Parameter (cxx L313-318) — V->Parameter(C).
pub fn parameter(v: &dyn HVertexBehavior, c: &dyn Curve2dAdaptor) -> f64 {
    v.parameter(c)
}

/// OCCT HasBeenSeen (cxx L320-323).
pub fn has_been_seen(_c: &dyn Curve2dAdaptor) -> bool {
    false
}

/// OCCT NbPoints (cxx L325-328).
pub fn nb_points(_c: &dyn Curve2dAdaptor) -> i32 {
    0
}

/// OCCT Value (cxx L330-337) — raises Standard_OutOfRange.
pub fn value(_c: &dyn Curve2dAdaptor, _index: i32) -> (Point3, f64, f64) {
    panic!("Standard_OutOfRange: Contap_HContTool::Value")
}

/// OCCT IsVertex (cxx L339-342).
pub fn is_vertex(_c: &dyn Curve2dAdaptor, _index: i32) -> bool {
    false
}

/// OCCT Vertex (cxx L344-349) — raises Standard_OutOfRange.
pub fn vertex(_c: &dyn Curve2dAdaptor, _index: i32) -> HVertex {
    panic!("Standard_OutOfRange: Contap_HContTool::Vertex")
}

/// OCCT NbSegments (cxx L351-354).
pub fn nb_segments(_c: &dyn Curve2dAdaptor) -> i32 {
    0
}

/// OCCT HasFirstPoint (cxx L356-359) — raises Standard_OutOfRange.
pub fn has_first_point(_c: &dyn Curve2dAdaptor, _index: i32) -> (bool, i32) {
    panic!("Standard_OutOfRange: Contap_HContTool::HasFirstPoint")
}

/// OCCT HasLastPoint (cxx L361-364) — raises Standard_OutOfRange.
pub fn has_last_point(_c: &dyn Curve2dAdaptor, _index: i32) -> (bool, i32) {
    panic!("Standard_OutOfRange: Contap_HContTool::HasLastPoint")
}

/// OCCT IsAllSolution (cxx L366-370).
pub fn is_all_solution(_c: &dyn Curve2dAdaptor) -> bool {
    false
}

// Silence the unused DVec3 import if Point3 aliases to it.
#[allow(unused)]
fn _dvec3_shape() -> DVec3 {
    DVec3::ZERO
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle2d, CylindricalSurface, Surface3};

    /// OCCT anchor: NbSamplesU/V dispatch for cylinder (10 / 15) and plane
    /// (2 / 2) (cxx L26-115).
    #[test]
    fn hcont_tool_sample_counts() {
        let cyl = GeomSurfaceAdapter::new(Surface3::Cylinder(CylindricalSurface {
            origin: DVec3::ZERO,
            axis: DVec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
            ref_dir: DVec3::new(1.0, 0.0, 0.0),
            y_dir: None,
        }));
        assert_eq!(nb_samples_u(&cyl, 0.0, 0.0), 10);
        assert_eq!(nb_samples_v(&cyl, 0.0, 0.0), 15);

        let pln = GeomSurfaceAdapter::new(Surface3::Plane(rcad_kernel::geom::Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(0.0, 0.0, 1.0),
            u_dir: DVec3::new(1.0, 0.0, 0.0),
            v_dir: DVec3::new(0.0, 1.0, 0.0),
        }));
        assert_eq!(nb_samples_u(&pln, 0.0, 0.0), 2);
        assert_eq!(nb_samples_v(&pln, 0.0, 0.0), 2);
    }

    /// OCCT anchor: NbSamplePoints + SamplePoint on a plane with finite
    /// bounds produce the 5 fixed barycentric samples (cxx L117-228).
    #[test]
    fn hcont_tool_sample_points() {
        let pln = GeomSurfaceAdapter::with_domain(
            Surface3::Plane(rcad_kernel::geom::Plane {
                origin: DVec3::ZERO,
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            }),
            [0.0, 4.0, 0.0, 4.0],
        );
        assert_eq!(nb_sample_points(&pln), 5);
        // cxx L209-210: U = 0.75*uinf + 0.25*usup = 1, V likewise.
        let (u1, v1) = sample_point(&pln, 1);
        assert!((u1 - 1.0).abs() < 1e-12 && (v1 - 1.0).abs() < 1e-12);
        // cxx L217-218: U = 0.25*uinf + 0.75*usup = 3, V = 1.
        let (u3, v3) = sample_point(&pln, 3);
        assert!((u3 - 3.0).abs() < 1e-12 && (v3 - 1.0).abs() < 1e-12);
        // default: the center.
        let (u5, v5) = sample_point(&pln, 5);
        assert!((u5 - 2.0).abs() < 1e-12 && (v5 - 2.0).abs() < 1e-12);
    }

    /// OCCT anchor: Project on a circle — the point (5,0) projects to
    /// (2,0) at parameter 0 (cxx L269-304 with EPCOfExtPC2d).
    #[test]
    fn hcont_tool_project_circle() {
        let circle = rcad_kernel::geom::Curve2d::Circle(Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 2.0,
        });
        let (param, pt) = project(&circle, DVec2::new(5.0, 0.0)).expect("Project done");
        assert!((pt.x - 2.0).abs() < 1e-4 && pt.y.abs() < 1e-4);
        let pu = param.rem_euclid(std::f64::consts::TAU);
        assert!(pu.abs() < 1e-3 || (pu - std::f64::consts::TAU).abs() < 1e-3);
    }
}
