//! The `Geom_OffsetSurface`-level re-hosts of the OCCT offset surface class,
//! split out of [`crate::geom::offset_surface_utils`] to keep each file under
//! the repository's size guideline.  The `Geom_OffsetSurfaceUtils.pxx`
//! namespace itself (and its `EvaluateD0/D1/D2/DN` leaves) stays in
//! `offset_surface_utils.rs`; this module carries the class bodies:
//!
//!   - `Geom_OffsetSurface::SetBasisSurface` unwrap (Geom_OffsetSurface.cxx
//!     L150-170) — [`offset_basis_and_value`] — plus the `myOscSurf`
//!     construction condition (L258-267) — [`offset_surface_osculating`];
//!   - `Geom_OffsetSurface::Surface()` (cxx L867-993) — the equivalent
//!     non-offset surface — [`offset_equivalent_surface`];
//!   - `Geom_OffsetSurface::EvalD0/EvalD1/EvalD2/EvalDN` (cxx L338-358 /
//!     L362-389 / L393-424 / L466-497) — [`offset_surface_eval_d0`] /
//!     [`offset_surface_eval_d1`] (without the `myEvalRep` short circuit) and
//!     the rcad [`OffsetSurface`] payload re-hosts
//!     [`offset_payload_eval_d0`] / [`offset_payload_eval_d1`] /
//!     [`offset_payload_eval_d2`] / [`offset_payload_eval_dn`] (with it).
//!
//! Architecture differences (rcad value model): `occ::handle(Geom_Surface)`
//! maps to [`Surface3`]; the OCCT null handles map to `None`; the
//! `myEvalRep` / `myOscSurf` members are recomputed on demand (see the
//! `offset_surface_utils.rs` header).

use glam::DVec3;

use super::osculating_surface::OsculatingSurface;
use super::offset_surface_utils::{
    eval_dn, evaluate_d0, evaluate_d1, evaluate_d2, evaluate_dn, ResD1, ResD2,
};
use crate::core::precision::CONFUSION;
use crate::geom::{
    ConicalSurface, CylindricalSurface, OffsetSurface, Plane, SphericalSurface, Surface3,
    SurfaceEval, ToroidalSurface, TrimmedSurface,
};

/// OCCT `Geom_OffsetSurface::SetBasisSurface` unwrap (cxx L150-170): the
/// nested `Geom_RectangularTrimmedSurface` / `Geom_OffsetSurface` wrappers are
/// unwrapped and the offset values are accumulated into a single
/// `offsetValue`.  Returns `(aCheckingSurf, accumulated_offset)`.
pub fn offset_basis_and_value(the_surf: &Surface3, the_offset: f64) -> (Surface3, f64) {
    let mut a_checking_surf = the_surf.clone();
    let mut offset_value = the_offset;
    while matches!(
        a_checking_surf,
        Surface3::Trimmed(_) | Surface3::Offset(_)
    ) {
        if let Surface3::Trimmed(t) = &a_checking_surf {
            let basis = t.basis.as_ref().clone();
            a_checking_surf = basis;
        }
        if let Surface3::Offset(o) = &a_checking_surf {
            let basis = o.basis.as_ref().clone();
            offset_value += o.offset_distance;
            a_checking_surf = basis;
        }
    }
    (a_checking_surf, offset_value)
}

/// OCCT `Geom_OffsetSurface::Surface()` (cxx L867-993) — the equivalent
/// non-offset surface of this offset surface, or `None` when none exists.
///
/// The rcad counterpart of the OCCT `myEvalRep` member
/// (`GeomEval_RepSurfaceDesc::Full`, built by `makeFullSurfaceRep(Surface())`
/// in `SetBasisSurface` L255-256) is this recomputation; `directRepSurface`
/// (cxx L83-95) is therefore `offset_equivalent_surface(...) != None` for a
/// freshly constructed surface.
pub fn offset_equivalent_surface(the_surf: &Surface3, the_offset: f64) -> Option<Surface3> {
    if the_offset == 0.0 {
        // Direct case - no offset.
        return Some(the_surf.clone());
    }

    let tol = CONFUSION;
    let mut result: Option<Surface3> = None;

    // Handle trimmed surfaces - extract the basis surface and bounds.
    let (base, is_trimmed, u1, u2, v1, v2) = match the_surf {
        Surface3::Trimmed(t) => {
            let b = t.basis.as_ref();
            (b.clone(), true, t.trim[0], t.trim[1], t.trim[2], t.trim[3])
        }
        other => (other.clone(), false, 0., 0., 0., 0.),
    };

    // Handle canonical surfaces - compute the equivalent offset surface.
    // For direct orientation, offset is along the outward normal; for
    // indirect, it is reversed.
    match &base {
        Surface3::Plane(p) => {
            // Plane normal is already available as Position().Direction().
            let t = p.normal * the_offset;
            // OCCT Geom_Plane::Translated(T) — the translated geometry.
            result = Some(Surface3::Plane(Plane {
                origin: p.origin + t,
                ..p.clone()
            }));
        }
        Surface3::Cylinder(c) => {
            // OCCT: gp_Ax3 Axis = C->Position(); aSign = Axis.Direct() ? 1 : -1.
            let a_sign = ax3_direct(c.y_dir, c.ref_dir, c.axis);
            let radius = c.radius + a_sign * the_offset;
            if radius >= tol {
                result = Some(Surface3::Cylinder(CylindricalSurface {
                    radius,
                    ..c.clone()
                }));
            } else if radius <= -tol {
                // Negative radius: flip the X-axis to reverse the normal
                // orientation.
                result = Some(Surface3::Cylinder(CylindricalSurface {
                    radius: -radius,
                    ref_dir: -c.ref_dir,
                    y_dir: c.y_dir.map(|y| -y),
                    ..c.clone()
                }));
            }
            // else: degenerate surface - radius is too small.
        }
        Surface3::Cone(c) => {
            // OCCT: gp_Ax3 anAxis = C->Position(); aSign = anAxis.Direct() ? 1 : -1.
            let a_sign = ax3_direct(None, c.ref_dir, c.axis);
            let an_alpha = c.half_angle_rad;
            let a_cos = an_alpha.cos();
            let a_sin = an_alpha.sin();
            let a_radius = c.radius + a_sign * the_offset * a_cos;
            if a_radius >= 0. {
                // Translate the apex along the axis by the offset component
                // (anAxis.Translate(aZ) moves the gp_Ax3 location).
                let a_z = c.axis * (-a_sign * the_offset * a_sin);
                result = Some(Surface3::Cone(ConicalSurface {
                    apex: c.apex + a_z,
                    radius: a_radius,
                    ..c.clone()
                }));
            }
            // else: degenerate surface - radius is negative.
        }
        Surface3::Sphere(s) => {
            let a_sign = ax3_direct(None, s.ref_dir, s.axis);
            let radius = s.radius + a_sign * the_offset;
            if radius >= tol {
                result = Some(Surface3::Sphere(SphericalSurface {
                    radius,
                    ..s.clone()
                }));
            } else if radius <= -tol {
                // Negative radius: flip both X and Z axes to reverse the
                // normal orientation.
                result = Some(Surface3::Sphere(SphericalSurface {
                    radius: -radius,
                    axis: -s.axis,
                    ref_dir: -s.ref_dir,
                    ..s.clone()
                }));
            }
            // else: degenerate surface - radius is too small.
        }
        Surface3::Torus(t) => {
            let major_radius = t.major_radius;
            let a_sign = ax3_direct(None, t.ref_dir, t.axis);
            let minor_radius = t.minor_radius + a_sign * the_offset;
            // Only handle the non-self-intersecting torus
            // (MinorRadius <= MajorRadius).
            if minor_radius >= tol && minor_radius <= major_radius {
                result = Some(Surface3::Torus(ToroidalSurface {
                    minor_radius,
                    ..t.clone()
                }));
            } else if minor_radius <= -tol && -minor_radius <= major_radius {
                // Negative minor radius: flip the X-axis to reverse the
                // normal orientation.
                result = Some(Surface3::Torus(ToroidalSurface {
                    minor_radius: -minor_radius,
                    ref_dir: -t.ref_dir,
                    ..t.clone()
                }));
            }
            // else: degenerate or self-intersecting torus - no equivalent
            // surface.
        }
        _ => {}
    }

    // Trim the result if the basis surface was trimmed.
    if is_trimmed {
        if let Some(b) = result.take() {
            result = Some(Surface3::Trimmed(TrimmedSurface::new(b, u1, u2, v1, v2)));
        }
    }

    result
}

/// OCCT `gp_Ax3::Direct()` for the rcad surface frames: the frame
/// (X = `ref_dir`, Y, Z = `axis`) is right-handed (gp_Ax3.hxx L57-58
/// `Y = Z ^ X` for the direct form; the indirect form reverses Y,
/// gp_Ax3.cxx L44-46).
///
/// Architecture difference: only [`CylindricalSurface`] carries an explicit
/// frame Y direction in the rcad payload (`y_dir`), so a cylinder handedness
/// is decoded from `det(ref_dir, y_dir, axis)`; the
/// [`ConicalSurface`] / [`SphericalSurface`] / [`ToroidalSurface`] payloads
/// store `ref_dir` + `axis` only and can therefore represent the right-handed
/// frames exclusively — the OCCT `aSign = -1` arms of the cone / sphere /
/// torus are unreachable in rcad.
fn ax3_direct(y_dir: Option<DVec3>, ref_dir: DVec3, axis: DVec3) -> f64 {
    match y_dir {
        None => 1.0,
        Some(y) => {
            if ref_dir.dot(y.cross(axis)) > 0.0 {
                1.0
            } else {
                -1.0
            }
        }
    }
}

/// OCCT `Geom_OffsetSurface::SetBasisSurface` osculating-surface construction
/// (cxx L258-267): `myOscSurf` is built only for a BSpline / Bezier basis
/// (after the wrapper unwrap), with the hard-coded `Precision::Confusion()`
/// tolerance of `SetBasisSurface`.  `None` is the OCCT null `myOscSurf`.
pub fn offset_surface_osculating(the_basis: &Surface3) -> Option<OsculatingSurface> {
    // OCCT: constexpr double Tol = Precision::Confusion();
    let tol = CONFUSION;
    match the_basis {
        Surface3::BSpline(_) | Surface3::Bezier(_) => {
            Some(OsculatingSurface::with_surface(the_basis, tol))
        }
        _ => None,
    }
}

/// OCCT `Geom_OffsetSurface::EvalD0` (cxx L338-358) — the point of the offset
/// surface, without the `GeomEval_RepUtils::TryEvalSurfaceD0(myEvalRep, ...)`
/// short circuit (that member is recomputed by
/// [`offset_equivalent_surface`]; the caller decides whether to take it).
pub fn offset_surface_eval_d0(
    the_basis: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_u: f64,
    the_v: f64,
) -> Option<DVec3> {
    evaluate_d0(the_u, the_v, the_basis, the_offset, the_osc_query)
}

/// OCCT `Geom_OffsetSurface::EvalD1` (cxx L362-389) — the point + first
/// partials of the offset surface, without the
/// `GeomEval_RepUtils::TryEvalSurfaceD1(myEvalRep, ...)` short circuit.
pub fn offset_surface_eval_d1(
    the_basis: &Surface3,
    the_offset: f64,
    the_osc_query: Option<&OsculatingSurface>,
    the_u: f64,
    the_v: f64,
) -> Option<ResD1> {
    evaluate_d1(the_u, the_v, the_basis, the_offset, the_osc_query)
}

/// OCCT `Geom_OffsetSurface::Offset()`-driven payload convenience: the
/// rcad `Geom_OffsetSurface`-equivalent evaluation of a rcad
/// [`OffsetSurface`] payload (the `offset_basis_and_value` unwrap + the
/// `offset_surface_osculating` member rebuild).
pub fn offset_payload_osculating(of: &OffsetSurface) -> Option<OsculatingSurface> {
    let (a_checking_surf, _offset_value) =
        offset_basis_and_value(of.basis.as_ref(), of.offset_distance);
    offset_surface_osculating(&a_checking_surf)
}

/// The rcad [`OffsetSurface`] payload re-host of the prologue shared by
/// `Geom_OffsetSurface::EvalD0` / `EvalD1` / `EvalD2` / `EvalD3` / `EvalDN`
/// (cxx L338-497): the `Geom_OffsetSurface::SetBasisSurface` unwrap
/// ([`offset_basis_and_value`]) that produces `basisSurf` / `offsetValue`, then
/// the `GeomEval_RepUtils::TryEvalSurfaceD*` short circuit — `myEvalRep` is the
/// `makeFullSurfaceRep(Surface())` representation of the unwrapped basis, so it
/// is present exactly when [`offset_equivalent_surface`] yields a surface.
///
/// Returns `(basisSurf, offsetValue, equivalentSurface)`.  `myOscSurf` (the
/// member rebuilt by [`offset_surface_osculating`], present only for a BSpline
/// / Bezier basis) is recomputed by the callers that need it.
fn offset_payload_prologue(of: &OffsetSurface) -> (Surface3, f64, Option<Surface3>) {
    let (a_checking_surf, offset_value) =
        offset_basis_and_value(of.basis.as_ref(), of.offset_distance);
    let a_equiv = offset_equivalent_surface(&a_checking_surf, offset_value);
    (a_checking_surf, offset_value, a_equiv)
}

/// OCCT `Geom_OffsetSurface::EvalD0` (cxx L338-358) over the rcad
/// [`OffsetSurface`] payload.
pub fn offset_payload_eval_d0(of: &OffsetSurface, the_u: f64, the_v: f64) -> DVec3 {
    let (a_basis, offset_value, a_equiv) = offset_payload_prologue(of);
    if let Some(a_res) = a_equiv {
        // GeomEval_RepUtils::TryEvalSurfaceD0(myEvalRep, U, V, aRes).
        return SurfaceEval::point_at(&a_res, the_u, the_v);
    }
    let a_osc = offset_surface_osculating(&a_basis);
    offset_surface_eval_d0(&a_basis, offset_value, a_osc.as_ref(), the_u, the_v)
        .expect("Geom_UndefinedValue: Geom_OffsetSurface::EvalD0")
}

/// OCCT `Geom_OffsetSurface::EvalD1` (cxx L362-389) over the rcad
/// [`OffsetSurface`] payload.
pub fn offset_payload_eval_d1(of: &OffsetSurface, the_u: f64, the_v: f64) -> ResD1 {
    let (a_basis, offset_value, a_equiv) = offset_payload_prologue(of);
    if let Some(a_res) = a_equiv {
        // GeomEval_RepUtils::TryEvalSurfaceD1(myEvalRep, U, V, aRes).
        let (point, d1u, d1v) = SurfaceEval::derivatives(&a_res, the_u, the_v);
        return ResD1 { point, d1u, d1v };
    }
    let a_osc = offset_surface_osculating(&a_basis);
    offset_surface_eval_d1(&a_basis, offset_value, a_osc.as_ref(), the_u, the_v)
        .expect("Geom_UndefinedDerivative: Geom_OffsetSurface::EvalD1")
}

/// OCCT `Geom_OffsetSurface::EvalD2` (cxx L393-424) over the rcad
/// [`OffsetSurface`] payload.
pub fn offset_payload_eval_d2(of: &OffsetSurface, the_u: f64, the_v: f64) -> ResD2 {
    let (a_basis, offset_value, a_equiv) = offset_payload_prologue(of);
    if let Some(a_res) = a_equiv {
        // GeomEval_RepUtils::TryEvalSurfaceD2(myEvalRep, U, V, aRes).
        let d2 = SurfaceEval::derivatives2(&a_res, the_u, the_v);
        return ResD2 {
            point: d2.0,
            d1u: d2.1,
            d1v: d2.2,
            d2u: d2.3,
            d2v: d2.5,
            d2uv: d2.4,
        };
    }
    let a_osc = offset_surface_osculating(&a_basis);
    evaluate_d2(the_u, the_v, &a_basis, offset_value, a_osc.as_ref())
        .expect("Geom_UndefinedDerivative: Geom_OffsetSurface::EvalD2")
}

/// OCCT `Geom_OffsetSurface::EvalDN` (cxx L466-497) over the rcad
/// [`OffsetSurface`] payload.
pub fn offset_payload_eval_dn(of: &OffsetSurface, the_u: f64, the_v: f64, nu: i32, nv: i32) -> DVec3 {
    // OCCT: if (Nu + Nv < 1 || Nu < 0 || Nv < 0) throw Geom_UndefinedDerivative.
    assert!(
        nu + nv >= 1 && nu >= 0 && nv >= 0,
        "Geom_UndefinedDerivative: Geom_OffsetSurface::EvalDN"
    );
    let (a_basis, offset_value, a_equiv) = offset_payload_prologue(of);
    if let Some(a_res) = a_equiv {
        // GeomEval_RepUtils::TryEvalSurfaceDN(myEvalRep, U, V, Nu, Nv, aRes).
        return eval_dn(&a_res, the_u, the_v, nu, nv);
    }
    let a_osc = offset_surface_osculating(&a_basis);
    evaluate_dn(
        the_u,
        the_v,
        nu,
        nv,
        &a_basis,
        offset_value,
        a_osc.as_ref(),
    )
    .expect("Geom_UndefinedDerivative: Geom_OffsetSurface::EvalDN")
}
