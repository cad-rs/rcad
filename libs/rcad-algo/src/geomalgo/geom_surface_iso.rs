//! OCCT Geom_RectangularTrimmedSurface (ModelingData/TKG3d/Geom) — the rcad
//! translation of the trimmed-surface construction, shared by the two isoline
//! consumers (`geomalgo::geom_lib_iso_line` = GeomLib::buildC3dOnIsoLine and
//! `geomalgo::approx_curve_on_surface` = Approx_CurveOnSurface).
//!
//! History: this file used to carry a second body for the
//! `Geom_Surface::UIso` / `VIso` virtual dispatch (one copy of the dispatch;
//! another partial copy lived in `offset::brep_offset_tool_iso`).  OCCT carries
//! that dispatch as one virtual function per concrete surface class, so the
//! duplicate bodies violated the duplicate-implementation convergence rule
//! (see `docs/module-map.md`): the complete, 1:1 line-by-line translation with
//! OCCT anchors is the canonical body in
//! `crate::brep_fill::brep_fill_sweep::{surface_uiso, surface_viso}`, and the
//! copies were deleted (consumers redirected onto the canonical body).

use rcad_kernel::geom::{Surface3, TrimmedSurface};

/// OCCT Geom_RectangularTrimmedSurface(S, U1, U2, V1, V2, USense, VSense)
/// (Geom_RectangularTrimmedSurface.cxx L67-112) — the nested trimmed basis is
/// killed and the resulting surface carries `isutrimmed = isvtrimmed = true`.
pub(crate) fn surface_rectangular_trimmed(
    surf: &Surface3,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Surface3 {
    // OCCT: kill trimmed basis surfaces.
    let basis = match surf {
        Surface3::Trimmed(t) => (*t.basis).clone(),
        other => other.clone(),
    };
    Surface3::Trimmed(TrimmedSurface::new(basis, u1, u2, v1, v2))
}
