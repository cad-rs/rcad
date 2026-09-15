//! OCCT BRepBuilderAPI_Sewing.cxx L170-327 — the closed-surface analysis
//! helpers: the file-scope statics `IsClosedShape` (L170-215) and
//! `IsClosedByIsos` (L217-259), and the virtual methods
//! `IsUClosedSurface` (L261-293) / `IsVClosedSurface` (L295-327).
//!
//! Architecture differences:
//! - OCCT `BRep_Tool::CurveOnSurface(edge, surf, loc, f2d, l2d)` matches the
//!   pcurve BY SURFACE; the rcad pcurves are keyed by the owning face
//!   (mod.rs #1), so the rcad `IsUClosed/IsVClosedSurface` take the owning
//!   face in addition to the surface (every OCCT call site has it).
//! - `IsClosedByIsos` evaluates the iso sample points through the surface
//!   evaluation directly (the OCCT `surf->UIso/VIso` handles are not
//!   materialized; an iso point IS the surface evaluation at the fixed iso
//!   parameter, and the iso curve's parameter range is the surface's
//!   other-direction bounds — identical D0 values).

use glam::DVec3;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, Surface3};
use rcad_kernel::topo::topods::BRep;
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

use super::{brep_tool_degenerated, BRepBuilderAPISewing};

/// OCCT static IsClosedShape(theshape, v1, v2) (cxx L170-215).
pub(crate) fn is_closed_shape(
    _brep: &BRep,
    theshape: &Shape,
    v1: &Shape,
    v2: &Shape,
) -> bool {
    // OCCT L171: double TotLength = 0.0;
    let mut tot_length = 0.0f64;
    // OCCT L172-174: TopExp_Explorer aexp; for (aexp.Init(theshape, TopAbs_EDGE); ...).
    for aexp in bat::explorer(theshape, rcad_kernel::topo::topods::ShapeType::Edge, rcad_kernel::topo::topods::ShapeType::Shape) {
        // OCCT L176: TopoDS_Edge aedge = TopoDS::Edge(aexp.Current());
        let aedge = aexp;
        // OCCT L177-180: if (aedge.IsNull()) continue;
        if aedge.is_null() {
            continue;
        }
        // OCCT L182: TopExp::Vertices(aedge, ve1, ve2).
        let (ve1, ve2) = bat::top_exp_vertices_raw(&aedge);
        let (ve1, ve2) = match (ve1, ve2) {
            (Some(a), Some(b)) => (a, b),
            // OCCT proceeds with null ve1/ve2 only through the IsSame tests;
            // a null vertex is same to nothing, so the guard below falls
            // through to `continue` exactly like OCCT.
            (a, b) => (
                a.unwrap_or_else(Shape::null),
                b.unwrap_or_else(Shape::null),
            ),
        };
        // OCCT L183-186: if (!ve1.IsSame(v1) && !ve1.IsSame(v2)) continue;
        if !ve1.is_same(v1) && !ve1.is_same(v2) {
            continue;
        }
        // OCCT L187-190: if (BRep_Tool::Degenerated(aedge)) continue;
        if brep_tool_degenerated(&aedge) {
            continue;
        }
        // OCCT L192-193: BRep_Tool::Curve(aedge, first, last).
        let c3d = match bat::brep_tool_curve(&aedge) {
            Some((c, f, l)) => (c, f, l),
            None => continue,
        };
        // OCCT L194-205: accumulate the edge length; break when ve2 closes
        // on v1/v2.
        tot_length += super::gcpnts_abscissa_length(&c3d.0, c3d.1, c3d.2);
        if ve2.is_same(v1) || ve2.is_same(v2) {
            break;
        }
    }
    // OCCT L207-213.
    if tot_length > 0.0 {
        let p1 = bat::brep_tool_pnt(v1).unwrap_or(DVec3::ZERO);
        let p2 = bat::brep_tool_pnt(v2).unwrap_or(DVec3::ZERO);
        return p1.distance(p2) < tot_length / (1.2 * std::f64::consts::PI);
    }
    false
}

/// OCCT static IsClosedByIsos(thesurf, acrv2d, f2d, l2d, isUIsos)
/// (cxx L217-259).
pub(crate) fn is_closed_by_isos(
    thesurf: &Surface3,
    acrv2d: &rcad_kernel::geom::Curve2d,
    f2d: f64,
    l2d: f64,
    is_u_isos: bool,
) -> bool {
    // OCCT L221: bool isClosed = false;
    let mut is_closed = false;

    // OCCT L223-227: psurf1/psurf2 — the iso parameters at the pcurve ends
    // (periodic curves keep f2d/l2d, the others clamp onto the pcurve
    // domain: max(f2d, FirstParameter) / min(l2d, LastParameter)).
    let [pc_first, pc_last] = acrv2d.default_domain();
    let psurf1 = if acrv2d.is_periodic() {
        acrv2d.point_at(f2d)
    } else {
        acrv2d.point_at(f2d.max(pc_first))
    };
    let psurf2 = if acrv2d.is_periodic() {
        acrv2d.point_at(l2d)
    } else {
        acrv2d.point_at(l2d.min(pc_last))
    };
    // OCCT L229-238: aCrv1/aCrv2 — the iso curves at psurf1/psurf2.  The iso
    // parameter is the fixed direction coordinate; the iso curve's own
    // parameter range is the surface's other-direction bounds.
    let (iso1, iso2) = if is_u_isos {
        (psurf1.x, psurf2.x)
    } else {
        (psurf1.y, psurf2.y)
    };
    // The other-direction bounds (the iso curves' [FirstParameter,
    // LastParameter], cxx L239-242).
    let [u1, u2, v1, v2] = thesurf.default_domain();
    let (af1, al1, af2, al2) = if is_u_isos { (v1, v2, v1, v2) } else { (u1, u2, u1, u2) };
    // OCCT L243-248: p11/p1m/p12 and p21/p2m/p22 — the iso endpoints and
    // midpoints; each iso point is the surface evaluation at (fixed iso
    // parameter, iso curve parameter).
    let iso_point = |par: f64, iso: f64| -> DVec3 {
        if is_u_isos {
            thesurf.point_at(iso, par)
        } else {
            thesurf.point_at(par, iso)
        }
    };
    let p11 = iso_point(af1, iso1);
    let p1m = iso_point((af1 + al1) * 0.5, iso1);
    let p12 = iso_point(al1, iso1);
    let p21 = iso_point(af2, iso2);
    let p2m = iso_point((af2 + al2) * 0.5, iso2);
    let p22 = iso_point(al2, iso2);
    // OCCT L249-253.
    is_closed = ((p11 - p12).length() < (p11 - p1m).length() - rcad_kernel::core::precision::CONFUSION)
        && ((p21 - p22).length() < (p21 - p2m).length() - rcad_kernel::core::precision::CONFUSION);
    // OCCT L254.
    is_closed
}

impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing::IsUClosedSurface(surf, theEdge, theloc)
    /// (cxx L261-293) — defines if surface is U closed.
    pub(crate) fn is_u_closed_surface(
        &self,
        brep: &BRep,
        surf: &Surface3,
        the_edge: &Shape,
        the_face: &Shape,
        the_loc: u32,
    ) -> bool {
        let mut tmpsurf = surf.clone();
        // OCCT L263-274: unwrap the rectangular trimmed / offset basis.
        if let Surface3::Trimmed(ts) = &tmpsurf {
            tmpsurf = (*ts.basis).clone();
            return self.is_u_closed_surface(brep, &tmpsurf, the_edge, the_face, the_loc);
        } else if let Surface3::Offset(os) = &tmpsurf {
            tmpsurf = (*os.basis).clone();
            return self.is_u_closed_surface(brep, &tmpsurf, the_edge, the_face, the_loc);
        }
        // OCCT L276-290 (the else branch).
        let mut is_closed = tmpsurf.is_u_closed();
        if !is_closed {
            // OCCT L280-282: BRep_Tool::CurveOnSurface(edge, surf, theloc,
            // f2d, l2d) — the rcad read is keyed by the owning face.
            if let Some((acrv2d, f2d, l2d)) = bat::brep_tool_curve_on_surface(the_edge, the_face) {
                is_closed = is_closed_by_isos(&tmpsurf, &acrv2d, f2d, l2d, false);
            }
        }
        let _ = (the_loc, brep);
        // OCCT L291: return isClosed.
        is_closed
    }

    /// OCCT BRepBuilderAPI_Sewing::IsVClosedSurface(surf, theEdge, theloc)
    /// (cxx L295-327) — defines if surface is V closed.
    pub(crate) fn is_v_closed_surface(
        &self,
        brep: &BRep,
        surf: &Surface3,
        the_edge: &Shape,
        the_face: &Shape,
        the_loc: u32,
    ) -> bool {
        let mut tmpsurf = surf.clone();
        // OCCT L297-308: unwrap the rectangular trimmed / offset basis.
        if let Surface3::Trimmed(ts) = &tmpsurf {
            tmpsurf = (*ts.basis).clone();
            return self.is_v_closed_surface(brep, &tmpsurf, the_edge, the_face, the_loc);
        } else if let Surface3::Offset(os) = &tmpsurf {
            tmpsurf = (*os.basis).clone();
            return self.is_v_closed_surface(brep, &tmpsurf, the_edge, the_face, the_loc);
        }
        // OCCT L310-324 (the else branch).
        let mut is_closed = tmpsurf.is_v_closed();
        if !is_closed {
            if let Some((acrv2d, f2d, l2d)) = bat::brep_tool_curve_on_surface(the_edge, the_face) {
                is_closed = is_closed_by_isos(&tmpsurf, &acrv2d, f2d, l2d, true);
            }
        }
        let _ = (the_loc, brep);
        // OCCT L325: return isClosed.
        return is_closed;
    }
}

