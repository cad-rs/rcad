// OCCT HLRBRep_BSurfaceTool (TKHLR) — the static tool over the face
// adaptor (`BRepAdaptor_Surface`).
//
// HLRBRep_BSurfaceTool.hxx L36-170 + .lxx (one-liner delegates to the
// adaptor) + .cxx L22-153 (NbSamplesU/V). The remaining hxx statics are
// one-line forwards to the BRepAdaptor_Surface — carried by the
// `BRepAdaptorSurface` inherent methods (see `topalgo/brep_adaptor/surface.rs`).

use crate::topalgo::brep_adaptor::surface::BRepAdaptorSurface;

/// OCCT NbSamplesU(S) (cxx L22-64).
#[allow(clippy::match_like_matches_macro)]
pub fn nb_samples_u_total(s: &BRepAdaptorSurface<'_>) -> i32 {
    let nbs;
    let typ_s = s.get_type();
    match typ_s {
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane => {
            nbs = 2;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::BezierSurface => {
            nbs = (3 + s.nb_u_poles()) as i32;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::BSplineSurface => {
            let mut n = (s.nb_u_knots() * s.u_degree()) as i32;
            if n < 2 {
                n = 2;
            }
            nbs = n;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Torus => {
            nbs = 20;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfRevolution
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfExtrusion => {
            nbs = 10;
        }
        _ => {
            nbs = 10;
        }
    }
    nbs
}

/// OCCT NbSamplesV(S) (cxx L68-107).
#[allow(clippy::match_like_matches_macro)]
pub fn nb_samples_v_total(s: &BRepAdaptorSurface<'_>) -> i32 {
    let nbs;
    let typ_s = s.get_type();
    match typ_s {
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane => {
            nbs = 2;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::BezierSurface => {
            nbs = (3 + s.nb_v_poles()) as i32;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::BSplineSurface => {
            let mut n = (s.nb_v_knots() * s.v_degree()) as i32;
            if n < 2 {
                n = 2;
            }
            nbs = n;
        }
        crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::Torus
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfRevolution
        | crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfExtrusion => {
            nbs = 15;
        }
        _ => {
            nbs = 10;
        }
    }
    nbs
}

/// OCCT NbSamplesU(S, u1, u2) (cxx L111-130) — the `(u2 - u1) / (uf - ul)`
/// ratio uses the (negative) reversed full-domain denominator, verbatim.
pub fn nb_samples_u_range(s: &BRepAdaptorSurface<'_>, u1: f64, u2: f64) -> i32 {
    let nbs = nb_samples_u_total(s);
    let mut n = nbs;
    if nbs > 10 {
        let uf = s.first_u_parameter();
        let ul = s.last_u_parameter();
        n = (n as f64 * ((u2 - u1) / (uf - ul))) as i32;
        if n > nbs {
            n = nbs;
        }
        if n < 5 {
            n = 5;
        }
    }
    n
}

/// OCCT NbSamplesV(S, v1, v2) (cxx L134-153).
pub fn nb_samples_v_range(s: &BRepAdaptorSurface<'_>, v1: f64, v2: f64) -> i32 {
    let nbs = nb_samples_v_total(s);
    let mut n = nbs;
    if nbs > 10 {
        let vf = s.first_v_parameter();
        let vl = s.last_v_parameter();
        n = (n as f64 * ((v2 - v1) / (vf - vl))) as i32;
        if n > nbs {
            n = nbs;
        }
        if n < 5 {
            n = 5;
        }
    }
    n
}
