//! OCCT BRepLib::EncodeRegularity family (TKTopAlgo/BRepLib —
//! BRepLib.cxx L2104-2587) — the 1:1 translation:
//! - file-scope class `SurfaceProperties` (cxx L2104-2197),
//! - `BRepLib::ContinuityOfFaces` (cxx L2196-2421),
//! - file-scope static `EncodeRegularity` (cxx L2428-2526),
//! - `BRepLib::EncodeRegularity(S, Tol)` (cxx L2533-2537),
//! - `BRepLib::EncodeRegularity(E, F1, F2, Tol)` (cxx L2570-2587).
//!
//! The `EncodeRegularity(S, LE, Tol)` list overload (cxx L2545-2567) has no
//! rcad consumer yet; its only distinct statement — the pure-edge set that
//! feeds the static's `theEdgesToEncode` filter — is carried by the static's
//! optional filter parameter with the same OCCT anchor.
//!
//! Architecture notes:
//! - the OCCT `gp_Trsf` location transforms (`aLoc1.Transformation()`,
//!   `Transformed(...)`, `IsNegative()`) are recorded no-ops: the rcad
//!   surfaces/pcurves travel world-space (locations applied at read), the
//!   same convention as `brep_lib::same_parameter`;
//! - the OCCT `try { OCC_CATCH_SIGNALS ... } catch (Standard_Failure)`
//!   shells (cxx L2461-2525 / L2574-2585) have no rcad counterpart — the
//!   translated bodies raise through panics like every Standard_Failure
//!   translation in the codebase;
//! - `BRep_Tool::Continuity(E, F1, F2)` (BRep_Tool.cxx L1180-1190 ->
//!   L1223-1252) is re-hosted here with the pool/pool-free dual read (the
//!   `fillet::fillet_surf` local precedent), matching the regularity
//!   records by surface value (`BRep_CurveRepresentation::IsRegularity`
//!   reduces to the value match — the rcad surfaces travel as values).

use std::sync::Arc;

use glam::{DVec2, DVec3};

use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_locate_ext_pc::LocateExtPC;
use rcad_kernel::base::geom_lprop::SLProps;
use rcad_kernel::base::proj_lib::adaptor::{Adaptor3dCurve, CurveOnSurface};
use rcad_kernel::base::proj_lib::{CurveType, Geom2dCurveAdaptor, GeomSurfaceAdaptor};
use rcad_kernel::core::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Surface3, SurfaceEval};
use rcad_kernel::topo::topods::{
    edge_data_pool_free, face_surface_value, shape_is_in_pool, surface_same, BRep, BRepBuilder,
    CurveRepresentation, GeomAbsShape, Orientation, ShapeType, TEdgeData, TShape,
};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

// =========================================================================
// class SurfaceProperties (BRepLib.cxx L2104-2197).
// =========================================================================

/// OCCT BRepLib.cxx L2104-2197 — the surface-properties helper feeding the
/// continuity decision.  The `gp_Trsf` member and its `Transformed` /
/// `IsNegative` uses are the recorded no-ops (architecture note above).
struct SurfaceProperties<'a> {
    /// OCCT mySurfaceProps — `GeomLProp_SLProps(theSurface, 2,
    /// Precision::Confusion())` (cxx L2113).
    surface_props: SLProps<'a>,
    /// OCCT myCurve2d.
    curve2d: &'a Curve2d,
    /// OCCT myIsReversed — the face based on the surface is reversed.
    is_reversed: bool,
    /// OCCT myCurveTangent — the tangent vector to the pcurve in UV.
    curve_tangent: DVec2,
}

impl<'a> SurfaceProperties<'a> {
    /// OCCT ctor (cxx L2107-2118): `SurfaceProperties(theSurface,
    /// theSurfaceTrsf, theCurve2D, theReversed)`.
    fn new(
        the_surface: &'a Surface3,
        _the_surface_trsf: (),
        the_curve2d: &'a Curve2d,
        the_reversed: bool,
    ) -> Self {
        SurfaceProperties {
            surface_props: SLProps::from_surface(the_surface, 2, CONFUSION),
            curve2d: the_curve2d,
            is_reversed: the_reversed,
            curve_tangent: DVec2::ZERO,
        }
    }

    /// OCCT Calculate (cxx L2121-2126): `myCurve2d->D1(theParamOnCurve, aUV,
    /// myCurveTangent); mySurfaceProps.SetParameters(aUV.X(), aUV.Y())`.
    fn calculate(&mut self, the_param_on_curve: f64) {
        let a_uv = Curve2dEval::point_at(self.curve2d, the_param_on_curve);
        self.curve_tangent = Curve2dEval::derivative_at(self.curve2d, the_param_on_curve);
        self.surface_props.set_parameters(a_uv.x, a_uv.y);
    }

    /// OCCT Value (cxx L2129-2131): `mySurfaceProps.Value().Transformed(...)` —
    /// the transform is the recorded no-op.
    fn value(&self) -> DVec3 {
        self.surface_props.value()
    }

    /// OCCT Derivative (cxx L2134-2153): the surface derivative orthogonal
    /// to the curve tangent vector.
    fn derivative(&mut self) -> DVec3 {
        // OCCT L2136-2141: the direction orthogonal to the tangent vector
        // of the curve; the zero vector comes back under Confusion.
        let mut an_ortho = DVec2::new(-self.curve_tangent.y, self.curve_tangent.x);
        let a_len = an_ortho.length();
        if a_len < CONFUSION {
            return DVec3::ZERO;
        }
        an_ortho /= a_len;
        if self.is_reversed {
            // OCCT L2143-2144: anOrtho.Reverse().
            an_ortho = -an_ortho;
        }
        // OCCT L2146-2147: SetLinearForm(anOrtho.X(), D1U, anOrtho.Y(), D1V).
        self.surface_props.d1u() * an_ortho.x + self.surface_props.d1v() * an_ortho.y
    }

    /// OCCT Normal (cxx L2156-2159).
    fn normal(&mut self) -> DVec3 {
        self.surface_props
            .normal()
            .expect("StdFail_UndefinedValue: SurfaceProperties::Normal")
    }

    /// OCCT Curvature (cxx L2162-2187): the principal curvatures and
    /// directions, negated for the reversed face and the negative
    /// transformation (the transform member is the recorded no-op).
    fn curvature(
        &mut self,
        the_principal_dir1: &mut DVec3,
        the_curvature1: &mut f64,
        the_principal_dir2: &mut DVec3,
        the_curvature2: &mut f64,
    ) {
        // OCCT L2164: CurvatureDirections(thePrincipalDir1, thePrincipalDir2).
        let (a_dir1, a_dir2) = self
            .surface_props
            .curvature_directions()
            .expect("StdFail_UndefinedValue: GeomLProp_SLProps::CurvatureDirections");
        *the_principal_dir1 = a_dir1;
        *the_principal_dir2 = a_dir2;
        // OCCT L2165-2166: MaxCurvature / MinCurvature.
        *the_curvature1 = self.surface_props.max_curvature();
        *the_curvature2 = self.surface_props.min_curvature();
        if self.is_reversed {
            // OCCT L2167-2170.
            *the_curvature1 = -*the_curvature1;
            *the_curvature2 = -*the_curvature2;
        }
        // OCCT L2172-2175: if (mySurfaceTrsf.IsNegative()) — the transform
        // member is the recorded no-op (world-space convention).
        // OCCT L2177-2178: the direction transforms — recorded no-ops.
    }
}

// =========================================================================
// BRepLib::ContinuityOfFaces (BRepLib.cxx L2196-2421).
// =========================================================================

/// OCCT BRep_Tool::CurveOnSurface(E, F, First, Last) (BRep_Tool.cxx
/// L301-312 -> L327-374) — the pcurve of the edge in the face, with the
/// closed-surface PCurve2 answer for a REVERSED local edge.  The
/// pool/pool-free dual read follows the `same_parameter::edge_data`
/// convention; the pcurves map is the rcad store fallback for
/// projection-built edges.
fn curve_on_surface_in_face(
    the_brep: &BRep,
    the_edg: &Shape,
    the_face: &Shape,
) -> Option<(Curve2d, f64, f64)> {
    // OCCT L334: Eisreversed = (E.Orientation() == TopAbs_REVERSED) — the
    // local edge orientation (the call sites pass FORWARD faces, so the
    // L303-310 face-orientation reversal branch is structural).
    let a_e_is_reversed = the_edg.orientation == Orientation::Reversed;
    let a_face_key = (the_face.ptr_id(), the_face.location);

    if shape_is_in_pool(the_brep, the_edg) {
        let a_te = the_brep.edge(the_edg.clone());
        curve_from_representations(&a_te, a_e_is_reversed, &a_face_key)
    } else {
        let a_te = edge_data_pool_free(the_edg)?;
        curve_from_representations(a_te, a_e_is_reversed, &a_face_key)
    }
}

/// The representation scan of OCCT BRep_Tool.cxx L339-364 (the
/// `IsCurveOnSurface(S, loc)` match answers the face-key match — the rcad
/// representations key the pcurves by owning face).
fn curve_from_representations(
    a_te: &TEdgeData,
    the_e_is_reversed: bool,
    the_face_key: &(u64, u32),
) -> Option<(Curve2d, f64, f64)> {
    for a_cr in &a_te.representations {
        match a_cr {
            CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } if face == the_face_key => {
                // OCCT L356-361: GC->Range(First, Last); return GC->PCurve().
                return Some((pcurve.clone(), range[0], range[1]));
            }
            CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                pcurve2,
                range,
            } if face == the_face_key => {
                // OCCT L354: IsCurveOnClosedSurface() && Eisreversed ->
                // PCurve2(); else PCurve().
                let a_c = if the_e_is_reversed { pcurve2 } else { pcurve1 };
                return Some((a_c.clone(), range[0], range[1]));
            }
            _ => {}
        }
    }
    // The rcad pcurves map (the projection-store encoding consumed by
    // bat::brep_tool_curve_on_surface).
    a_te
        .pcurves
        .get(the_face_key)
        .map(|(a_c, a_f, a_l)| (a_c.clone(), *a_f, *a_l))
}

/// OCCT BRepLib::ContinuityOfFaces(theEdge, theFace1, theFace2, theAngleTol)
/// (BRepLib.cxx L2196-2421).
pub fn continuity_of_faces(
    the_brep: &BRep,
    the_edge: &Shape,
    the_face1: &Shape,
    the_face2: &Shape,
    the_angle_tol: f64,
) -> GeomAbsShape {
    // OCCT L2199: bool isSeam = theFace1.IsEqual(theFace2).
    let a_is_seam = the_face1.is_equal(the_face2);

    // OCCT L2213-2221: the closed-edge branch — find the edge in the
    // FORWARD face 1 (this edge will have correct orientation).
    let (a_curve1, a_curve2, an_edge_in_face2) = if !the_face1.is_same(the_face2)
        && bat::brep_tool_is_closed_on_surface(the_edge, the_face1)
        && bat::brep_tool_is_closed_on_surface(the_edge, the_face2)
    {
        // OCCT L2213-2214: aFace1.Orientation(TopAbs_FORWARD).
        let a_face1_fwd = bat::oriented(the_face1, Orientation::Forward);
        let mut an_edge_in_face1 = Shape::null();
        for a_e in bat::explorer(&a_face1_fwd, ShapeType::Edge, ShapeType::Shape) {
            if a_e.is_same(the_edge) {
                an_edge_in_face1 = a_e;
                break;
            }
        }
        if an_edge_in_face1.is_null() {
            // OCCT L2222-2224.
            return GeomAbsShape::C0;
        }
        // OCCT L2227: aCurve1 = BRep_Tool::CurveOnSurface(anEdgeInFace1,
        // aFace1, aFirst, aLast).
        let a_c1 = curve_on_surface_in_face(the_brep, &an_edge_in_face1, &a_face1_fwd);
        // OCCT L2228-2231: aFace2 = FORWARD; anEdgeInFace2 =
        // anEdgeInFace1.Reverse(); aCurve2 = CurveOnSurface(anEdgeInFace2,
        // aFace2, ...).
        let a_face2_fwd = bat::oriented(the_face2, Orientation::Forward);
        let an_edge_in_face2 = bat::reversed(&an_edge_in_face1);
        let a_c2 = curve_on_surface_in_face(the_brep, &an_edge_in_face2, &a_face2_fwd);
        (a_c1, a_c2, an_edge_in_face2)
    } else {
        // OCCT L2233-2243: the default branch — the pcurves of the edge on
        // the two faces (the seam case reverses the local edge).
        let a_c1 = curve_on_surface_in_face(the_brep, the_edge, the_face1);
        let an_edge_in_face2 = if the_face1.is_same(the_face2) {
            bat::reversed(the_edge)
        } else {
            the_edge.clone()
        };
        let a_c2 = curve_on_surface_in_face(the_brep, &an_edge_in_face2, the_face2);
        (a_c1, a_c2, an_edge_in_face2)
    };

    // OCCT L2246-2249: a null pcurve answers GeomAbs_C0.
    let (Some(a_curve1), Some(a_curve2)) = (a_curve1, a_curve2) else {
        return GeomAbsShape::C0;
    };

    // OCCT L2253-2255: the face surfaces (the locations are the recorded
    // no-ops — world-space convention).
    let (Some(a_surface1_raw), Some(a_surface2_raw)) = (
        face_surface_value(the_brep, the_face1),
        face_surface_value(the_brep, the_face2),
    ) else {
        // The OCCT null surface handle makes the elementary-surface checks
        // false below; the properties construction would throw — the rcad
        // encoding carries no surface value at all, the C0 early-out.
        return GeomAbsShape::C0;
    };

    // OCCT L2258-2264: unwrap Geom_RectangularTrimmedSurface to the basis.
    let a_surface1_basis;
    let a_surface1 = match &a_surface1_raw {
        Surface3::Trimmed(a_ts) => {
            a_surface1_basis = a_ts.basis.as_ref().clone();
            &a_surface1_basis
        }
        _ => &a_surface1_raw,
    };
    let a_surface2_basis;
    let a_surface2 = match &a_surface2_raw {
        Surface3::Trimmed(a_ts) => {
            a_surface2_basis = a_ts.basis.as_ref().clone();
            &a_surface2_basis
        }
        _ => &a_surface2_raw,
    };

    // OCCT L2267-2272: the seam edge on elementary surfaces is always CN.
    let a_is_elementary = matches!(
        a_surface1,
        Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_)
    ) && matches!(
        a_surface2,
        Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_)
    );
    if a_is_seam && a_is_elementary {
        return GeomAbsShape::CN;
    }

    // OCCT L2274-2275: SurfaceProperties(theSurface, theSurfaceTrsf,
    // theCurve2D, theFace.Orientation() == TopAbs_REVERSED).
    let mut a_sp1 = SurfaceProperties::new(
        a_surface1,
        (),
        &a_curve1.0,
        the_face1.orientation == Orientation::Reversed,
    );
    let mut a_sp2 = SurfaceProperties::new(
        a_surface2,
        (),
        &a_curve2.0,
        the_face2.orientation == Orientation::Reversed,
    );

    // OCCT L2278-2285: the edge 3d range, shrunk by eps on both ends.
    let (mut a_f, mut a_l) = bat::brep_tool_range(the_edge);
    let a_eps = (a_l - a_f) / 100.0;
    a_f += a_eps; // to avoid calculations on
    a_l -= a_eps; // points of pointed squares.

    // OCCT L2287: const double anAngleTol2 = theAngleTol * theAngleTol.
    let an_angle_tol2 = the_angle_tol * the_angle_tol;

    // OCCT L2289-2297: the derivative carriers.
    let mut a_der1 = DVec3::ZERO;
    let mut a_der2 = DVec3::ZERO;
    let mut a_sq_len1 = 0.0f64;
    let mut a_sq_len2 = 0.0f64;
    let mut a_crv_dir1 = [DVec3::ZERO; 2];
    let mut a_crv_dir2 = [DVec3::ZERO; 2];
    let mut a_crv_len1 = [0.0f64; 2];
    let mut a_crv_len2 = [0.0f64; 2];

    // OCCT L2299: GeomAbs_Shape aCont = (isElementary ? CN : C2).
    let mut a_cont = if a_is_elementary {
        GeomAbsShape::CN
    } else {
        GeomAbsShape::C2
    };
    let mut a_cur_cont;

    // OCCT L2301-2302: `Extrema_LocateExtPC ext; occ::handle<BRepAdaptor_Curve>
    // aHC2;` — the projector is built lazily on the first refinement pass
    // (the OCCT `aHC2.IsNull()` gate, L2338).  The construction is replayed
    // at each refinement pass over the same range — the
    // `section_placement` replay convention (the OCCT member handle
    // persists; the observable behavior is identical).
    let mut a_refined = false;

    // OCCT L2302: for (int i = 0; i <= 20 && aCont > GeomAbs_C0; i++).
    for a_i in 0..=20i32 {
        if a_cont <= GeomAbsShape::C0 {
            break;
        }
        // OCCT L2306: u = f + (l - f) * i / 20.
        let a_u = a_f + (a_l - a_f) * a_i as f64 / 20.0;

        // First suppose that this is sameParameter (cxx L2308).
        a_cur_cont = GeomAbsShape::C0;

        // OCCT L2315-2316: aSP1.Calculate(u); aSP2.Calculate(u).
        a_sp1.calculate(a_u);
        a_sp2.calculate(a_u);

        // OCCT L2318-2322.
        a_der1 = a_sp1.derivative();
        a_sq_len1 = a_der1.length_squared();
        a_der2 = a_sp2.derivative();
        a_sq_len2 = a_der2.length_squared();
        let mut a_is_smooth_suspect =
            a_der1.cross(a_der2).length_squared() <= an_angle_tol2 * a_sq_len1 * a_sq_len2;
        if a_is_smooth_suspect {
            // OCCT L2324-2333: the normals with the face-orientation
            // reversal; opposing normals answer C0.
            let mut a_normal1 = a_sp1.normal();
            if the_face1.orientation == Orientation::Reversed {
                a_normal1 = -a_normal1;
            }
            let mut a_normal2 = a_sp2.normal();
            if the_face2.orientation == Orientation::Reversed {
                a_normal2 = -a_normal2;
            }
            if a_normal1.dot(a_normal2) < 0.0 {
                return GeomAbsShape::C0;
            }
        }

        if !a_is_smooth_suspect {
            // OCCT L2336-2351: refine by projection — the pcurve adaptor on
            // the second surface, built on the first pass (the aHC2 null
            // gate).  The construction is replayed at each refinement pass
            // over the same range — the section_placement replay
            // convention (the OCCT member handle persists; the observable
            // behavior is identical).
            if !a_refined {
                a_refined = true;
            }
            // OCCT L2340: aHC2 = new BRepAdaptor_Curve(anEdgeInFace2,
            // theFace2) — the CurveOnSurface adaptor over the pcurve.
            let (a_pc2, a_f2, a_l2) = match curve_on_surface_in_face(
                the_brep,
                &an_edge_in_face2,
                &the_face2,
            ) {
                Some(a_c) => a_c,
                None => return GeomAbsShape::C0,
            };
            let a_cs = Box::new(CurveOnSurface::new(
                Arc::new(Geom2dCurveAdaptor::with_range(a_pc2, a_f2, a_l2)),
                Arc::new(GeomSurfaceAdaptor::new(a_surface2_raw.clone())),
            ));
            // OCCT L2341: ext.Initialize(*aHC2, f, l,
            // Precision::PConfusion()) — over the shrunk window.
            let a_cth = Box::new(CurveToolHandle::new(
                a_cs.as_ref(),
                CurveType::Other,
                a_cs.is_periodic(),
                a_cs.period(),
                a_cs.resolution(CONFUSION),
                a_cs.is_closed(),
            ));
            let mut a_ext = Box::new(LocateExtPC::new());
            a_ext.initialize(a_cth.as_ref(), a_f, a_l, PCONFUSION);
            // OCCT L2342: ext.Perform(aSP1.Value(), u).
            a_ext.perform(a_sp1.value(), a_u);
            if a_ext.is_done() && a_ext.is_min() {
                // OCCT L2344-2348: aSP2.Calculate(poc.Parameter()).
                let a_param = a_ext.point().param;
                a_sp2.calculate(a_param);
                a_der2 = a_sp2.derivative();
                a_sq_len2 = a_der2.length_squared();
            }
            a_is_smooth_suspect =
                a_der1.cross(a_der2).length_squared() <= an_angle_tol2 * a_sq_len1 * a_sq_len2;
        }

        if a_is_smooth_suspect {
            // OCCT L2353-2361: the G1/C1 classification.
            a_cur_cont = GeomAbsShape::G1;
            if (a_sq_len1.sqrt() - a_sq_len2.sqrt()).abs() < CONFUSION
                && a_der1.dot(a_der2) > SQUARE_CONFUSION
            {
                // <= check vectors are codirectional.
                a_cur_cont = GeomAbsShape::C1;
            }
        } else {
            // OCCT L2362-2365.
            return GeomAbsShape::C0;
        }

        // OCCT L2367-2370: no need further processing, because maximal
        // continuity is less than G2.
        if a_cont < GeomAbsShape::G2 {
            continue;
        }

        // OCCT L2372-2377: the principal curvatures on each surface (the
        // carriers are locals — the two mutable out-parameters cannot alias
        // the array slots).
        let (mut a_d10, mut a_l10, mut a_d11, mut a_l11) =
            (DVec3::ZERO, 0.0f64, DVec3::ZERO, 0.0f64);
        a_sp1.curvature(&mut a_d10, &mut a_l10, &mut a_d11, &mut a_l11);
        a_crv_dir1 = [a_d10, a_d11];
        a_crv_len1 = [a_l10, a_l11];
        let (mut a_d20, mut a_l20, mut a_d21, mut a_l21) =
            (DVec3::ZERO, 0.0f64, DVec3::ZERO, 0.0f64);
        a_sp2.curvature(&mut a_d20, &mut a_l20, &mut a_d21, &mut a_l21);
        a_crv_dir2 = [a_d20, a_d21];
        a_crv_len2 = [a_l20, a_l21];

        // OCCT L2378-2404: the G2/C2 classification — the directions of
        // principal curvatures parallel and the curvature values equal.
        for a_step in 0..=1i32 {
            let a_other = (1 - a_step) as usize;
            if a_crv_dir1[0].cross(a_crv_dir2[a_step as usize]).length_squared()
                <= SQUARE_CONFUSION
                && (a_crv_len1[0] - a_crv_len2[a_step as usize]).abs() < CONFUSION
                && a_crv_dir1[1].cross(a_crv_dir2[a_other]).length_squared() <= SQUARE_CONFUSION
                && (a_crv_len1[1] - a_crv_len2[a_other]).abs() < CONFUSION
            {
                if a_cur_cont == GeomAbsShape::C1
                    && a_crv_dir1[0].dot(a_crv_dir2[a_step as usize]) > CONFUSION
                    && a_crv_dir1[1].dot(a_crv_dir2[a_other]) > CONFUSION
                {
                    a_cur_cont = GeomAbsShape::C2;
                } else {
                    a_cur_cont = GeomAbsShape::G2;
                }
                break;
            }
        }

        // OCCT L2406-2409.
        if a_cur_cont < a_cont {
            a_cont = a_cur_cont;
        }
    }

    // OCCT L2412-2418: the elementary-surface C2 is totally CN.
    if a_is_elementary && a_cont == GeomAbsShape::C2 {
        a_cont = GeomAbsShape::CN;
    }
    a_cont
}

// =========================================================================
// OCCT BRep_Tool::Continuity(E, F1, F2) (BRep_Tool.cxx L1180-1190 ->
// L1223-1252) — the regularity-record reader (the canonical re-host; the
// fillet_surf / loc_ope_gluer local copies predate it).
// =========================================================================

/// The dual-path edge-data read (the `same_parameter::edge_data`
/// convention).
fn edge_data_of<'a>(the_brep: &'a BRep, the_e: &'a Shape) -> Option<&'a TEdgeData> {
    if shape_is_in_pool(the_brep, the_e) {
        Some(the_brep.edge(the_e.clone()))
    } else {
        edge_data_pool_free(the_e)
    }
}

/// OCCT BRep_Tool::Continuity(E, F1, F2) — the regularity record
/// (BRep_CurveOn2Surfaces) matching the two face surfaces;
/// GeomAbs_C0 when no record matches.  The OCCT location composition
/// (`L.Predivided(E.Location())`) reduces to the surface-value match (the
/// rcad surfaces travel as values).
fn brep_tool_continuity(
    the_brep: &BRep,
    the_e: &Shape,
    the_f1: &Shape,
    the_f2: &Shape,
) -> GeomAbsShape {
    // OCCT L1183-1186: S1 = Surface(F1, l1); S2 = Surface(F2, l2).
    let (Some(a_s1), Some(a_s2)) = (
        face_surface_value(the_brep, the_f1),
        face_surface_value(the_brep, the_f2),
    ) else {
        return GeomAbsShape::C0;
    };
    let Some(a_te) = edge_data_of(the_brep, the_e) else {
        return GeomAbsShape::C0;
    };
    // OCCT L1240-1248: the IsRegularity(S1, S2, l1, l2) walk.
    for a_cr in &a_te.representations {
        if let CurveRepresentation::CurveOn2Surfaces {
            surface1,
            surface2,
            continuity,
            ..
        } = a_cr
        {
            if (surface_same(surface1, &a_s1) && surface_same(surface2, &a_s2))
                || (surface_same(surface1, &a_s2) && surface_same(surface2, &a_s1))
            {
                // OCCT L1246: return cr->Continuity().
                return *continuity;
            }
        }
    }
    // OCCT L1250: return GeomAbs_C0.
    GeomAbsShape::C0
}

// =========================================================================
// static EncodeRegularity (BRepLib.cxx L2428-2526) and the two public
// overloads (L2533-2537 / L2570-2587).
// =========================================================================

/// OCCT static EncodeRegularity(theShape, theTolAng, theMap,
/// theEdgesToEncode = empty) (cxx L2428-2526) — codes the regularities on
/// all edges of the shape, boundary of two faces that do not have it.
fn encode_regularity_impl(
    the_brep: &mut BRep,
    the_shape: &Shape,
    the_tol_ang: f64,
    the_map: &mut std::collections::HashSet<(u64, u32)>,
    the_edges_to_encode: Option<&std::collections::HashSet<(u64, u32)>>,
) {
    // OCCT L2433-2439: nullify the location (the world-space convention —
    // the recorded no-op) and the twice-processing guard.
    let a_shape_key = (the_shape.ptr_id(), the_shape.location);
    if !the_map.insert(a_shape_key) {
        return; // do not need to process shape twice
    }

    // OCCT L2441-2450: the compound / compsolid recursion.
    if the_shape.shape_type() == ShapeType::Compound
        || the_shape.shape_type() == ShapeType::CompSolid
    {
        for a_child in bat::sub_shapes(the_shape) {
            encode_regularity_impl(the_brep, &a_child, the_tol_ang, the_map, the_edges_to_encode);
        }
        return;
    }

    // OCCT L2452-2456 (the try shell is the architecture note): the
    // EDGE/FACE ancestors map.
    let mut a_m: indexmap::IndexMap<(u64, u32), (Shape, Vec<Shape>)> = indexmap::IndexMap::new();
    crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
        the_shape,
        ShapeType::Edge,
        ShapeType::Face,
        &mut a_m,
    );

    // OCCT L2465-2516: the per-edge walk.
    for (_a_key, (a_e, a_faces)) in a_m.iter() {
        // OCCT L2468-2478: process only the edges from the list when the
        // filter is not empty (the LE-overload path).
        if let Some(a_filter) = the_edges_to_encode {
            if !a_filter.contains(&(_a_key.0, _a_key.1)) {
                continue;
            }
        }

        // OCCT L2481-2500: F1 = the first ancestor face, F2 = the first
        // ancestor face different from F1.
        let mut a_found = false;
        let mut a_f1 = Shape::null();
        let mut a_f2 = Shape::null();
        for a_it in a_faces {
            if a_f1.is_null() {
                a_f1 = a_it.clone();
            } else if !a_f1.is_same(a_it) {
                a_found = true;
                a_f2 = a_it.clone();
                break;
            }
        }

        // OCCT L2501-2512: the seam-edge detection — the same edge twice
        // in the face with different orientations.
        if !a_found && !a_f1.is_null() {
            let a_or_e = a_e.orientation;
            let a_f1_edges = bat::explorer(&a_f1, ShapeType::Edge, ShapeType::Shape);
            for a_cur_e in a_f1_edges {
                if a_e.is_same(&a_cur_e) && a_or_e != a_cur_e.orientation {
                    a_found = true;
                    a_f2 = a_f1.clone();
                    break;
                }
            }
        }

        // OCCT L2513-2515.
        if a_found {
            encode_regularity_edge(the_brep, a_e, &a_f1, &a_f2, the_tol_ang);
        }
    }
}

/// OCCT BRepLib::EncodeRegularity(const TopoDS_Shape& S, const double Tol)
/// (cxx L2533-2537).
pub fn encode_regularity(the_brep: &mut BRep, the_s: &Shape, the_tol: f64) {
    let mut a_map = std::collections::HashSet::new();
    encode_regularity_impl(the_brep, the_s, the_tol, &mut a_map, None);
}

/// OCCT BRepLib::EncodeRegularity(E, F1, F2, Tol) (cxx L2570-2587) — the
/// regularity between 2 faces connected by an edge.
pub fn encode_regularity_edge(
    the_brep: &mut BRep,
    the_e: &Shape,
    the_f1: &Shape,
    the_f2: &Shape,
    the_tol: f64,
) {
    // OCCT L2573: the Continuity gate.
    if brep_tool_continuity(the_brep, the_e, the_f1, the_f2) <= GeomAbsShape::C0 {
        // OCCT L2575-2577 (the catch shell is the architecture note).
        let a_cont = continuity_of_faces(the_brep, the_e, the_f1, the_f2, the_tol);
        // OCCT L2578: B.Continuity(E, F1, F2, aCont).
        let mut a_b = BRepBuilder::new();
        a_b.continuity(the_brep, the_e, the_f1, the_f2, a_cont);
    }
}

// =========================================================================
// Tests (hand-derived closed forms).
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Curve3, Line2d, Line3, Plane};

    /// One rectangular boundary edge: the 3d line from (x0, y0) along
    /// `dir`, with the plane pcurve (origin UV + the same direction).
    fn add_rect_edge(
        the_brep: &mut BRep,
        the_v1: Shape,
        the_v2: Shape,
        the_x0: f64,
        the_y0: f64,
        the_dir: DVec3,
        the_len: f64,
    ) -> Shape {
        the_brep.add_tedge(
            Some(Curve3::Line(Line3::new(
                DVec3::new(the_x0, the_y0, 0.0),
                the_dir,
            ))),
            the_v1,
            the_v2,
            [0.0, the_len],
        )
    }

    /// The plane pcurve representation pushed on the edge for the face
    /// (the global plane UV of the z = 0 plane).
    fn push_plane_pcurve(the_brep: &mut BRep, the_e: &Shape, the_face: &Shape) {
        let a_te = the_brep.edge(the_e.clone());
        let Some(Curve3::Line(a_ln)) = &a_te.curve else {
            panic!("expected a line edge");
        };
        let a_pc = Curve2d::Line(Line2d::new(
            DVec2::new(a_ln.origin.x, a_ln.origin.y),
            DVec2::new(a_ln.direction.x, a_ln.direction.y),
        ));
        let a_te2 = the_brep.edge(the_e.clone());
        let a_r = (a_te2.range[0], a_te2.range[1]);
        drop(a_te2);
        drop(a_te);
        the_brep
            .edge_mut_inplace(the_e.clone())
            .representations
            .push(CurveRepresentation::CurveOnSurface {
                face: (the_face.ptr_id(), the_face.location),
                pcurve: a_pc,
                range: [a_r.0, a_r.1],
            });
    }

    /// The shared construction: two coplanar rectangular faces on z = 0 —
    /// face 1 over [0,2]x[0,2], face 2 over [2,4]x[0,2] — sharing the ONE
    /// edge TShape at x = 2 (the same Shape in both wires, two
    /// CurveOnSurface representations).  Returns the compound of the two
    /// faces.
    fn make_coplanar_pair(the_brep: &mut BRep) -> Shape {
        let a_plane = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let mut a_b = BRepBuilder::new();

        let a_p00 = the_brep.add_tvertex(DVec3::new(0.0, 0.0, 0.0));
        let a_vs = the_brep.add_tvertex(DVec3::new(2.0, 0.0, 0.0));
        let a_p20 = the_brep.add_tvertex(DVec3::new(4.0, 0.0, 0.0));
        let a_p01 = the_brep.add_tvertex(DVec3::new(0.0, 2.0, 0.0));
        let a_ve = the_brep.add_tvertex(DVec3::new(2.0, 2.0, 0.0));
        let a_p21 = the_brep.add_tvertex(DVec3::new(4.0, 2.0, 0.0));

        // Face 1: bottom / shared / top / left.
        let a_bottom1 =
            add_rect_edge(the_brep, a_p00.clone(), a_vs.clone(), 0.0, 0.0, DVec3::X, 2.0);
        let a_shared = add_rect_edge(the_brep, a_vs.clone(), a_ve.clone(), 2.0, 0.0, DVec3::Y, 2.0);
        let a_top1 = add_rect_edge(the_brep, a_p01.clone(), a_ve.clone(), 0.0, 2.0, DVec3::X, 2.0);
        let a_left = add_rect_edge(the_brep, a_p00.clone(), a_p01.clone(), 0.0, 0.0, DVec3::Y, 2.0);
        let a_wire1 = a_b.build_wire(
            the_brep,
            vec![
                a_bottom1.clone(),
                a_shared.clone(),
                a_top1.clone(),
                a_left.clone(),
            ],
        );
        let a_face1 = a_b.make_face(the_brep, Some(a_plane.clone()), a_wire1);

        // Face 2: bottom / right / top / the SAME shared edge.
        let a_bottom2 =
            add_rect_edge(the_brep, a_vs.clone(), a_p20.clone(), 2.0, 0.0, DVec3::X, 2.0);
        let a_right = add_rect_edge(the_brep, a_p20.clone(), a_p21.clone(), 4.0, 0.0, DVec3::Y, 2.0);
        let a_top2 = add_rect_edge(the_brep, a_ve.clone(), a_p21.clone(), 2.0, 2.0, DVec3::X, 2.0);
        let a_wire2 = a_b.build_wire(
            the_brep,
            vec![
                a_bottom2.clone(),
                a_right.clone(),
                a_top2.clone(),
                a_shared.clone(),
            ],
        );
        let a_face2 = a_b.make_face(the_brep, Some(a_plane.clone()), a_wire2);

        // The pcurves (after the faces exist — the representations key on
        // the face shapes).
        for (a_e, a_f) in [
            (&a_bottom1, &a_face1),
            (&a_shared, &a_face1),
            (&a_top1, &a_face1),
            (&a_left, &a_face1),
            (&a_bottom2, &a_face2),
            (&a_right, &a_face2),
            (&a_top2, &a_face2),
            (&a_shared, &a_face2),
        ] {
            push_plane_pcurve(the_brep, a_e, a_f);
        }

        // The SHELL of the two faces — the OCCT static recurses through
        // COMPOUND/COMPSOLID only (cxx L2442-2450), so a compound of faces
        // would map each face alone with EMPTY ancestor lists; the shell is
        // the OCCT input form where the EDGE/FACE ancestors map sees both
        // faces.
        let a_shell = a_b.make_shell(the_brep);
        a_b.add_to_shell(the_brep, a_shell.clone(), a_face1.clone());
        a_b.add_to_shell(the_brep, a_shell.clone(), a_face2.clone());
        a_shell
    }

    /// The probe: the edge carrying two CurveOnSurface representations
    /// (the shared x = 2 edge).
    fn find_shared_edge(the_brep: &BRep, the_comp: &Shape) -> Shape {
        for a_f in bat::sub_shapes(the_comp) {
            for a_e in bat::explorer(&a_f, ShapeType::Edge, ShapeType::Shape) {
                if shape_is_in_pool(the_brep, &a_e) {
                    let a_te = the_brep.edge(a_e.clone());
                    let a_pcs = a_te
                        .representations
                        .iter()
                        .filter(|a_cr| {
                            matches!(a_cr, CurveRepresentation::CurveOnSurface { .. })
                        })
                        .count();
                    if a_pcs > 1 {
                        return a_e;
                    }
                }
            }
        }
        Shape::null()
    }

    /// Two coplanar forward faces (elementary surfaces): the orthogonal
    /// derivatives are equal (C1); the principal directions of the SAME
    /// plane are identical (dot = 1 > Confusion), so the codirectionality
    /// promotion (cxx L2391-2397) raises C1 to C2; aCont = min(CN, C2) = C2
    /// through the loop, and the elementary tail (cxx L2412-2418) promotes
    /// the C2 to CN.  (First expectation was G2 — wrong: the zero-curvature
    /// directions of one plane are NOT zero vectors, they coincide.)
    #[test]
    fn coplanar_forward_faces_code_cn() {
        let mut a_brep = BRep::new();
        let a_comp = make_coplanar_pair(&mut a_brep);
        encode_regularity(&mut a_brep, &a_comp, CONFUSION);

        let a_shared = find_shared_edge(&a_brep, &a_comp);
        assert!(
            !a_shared.is_null(),
            "the shared edge must carry two pcurve representations"
        );
        let a_te = a_brep.edge(a_shared.clone());
        let a_rec = a_te.representations.iter().find_map(|a_cr| match a_cr {
            CurveRepresentation::CurveOn2Surfaces { continuity, .. } => Some(*continuity),
            _ => None,
        });
        assert_eq!(a_rec, Some(GeomAbsShape::CN), "the coplanar pair codes CN");
    }

    /// The reversed second face: opposing normals answer GeomAbs_C0
    /// directly (cxx L2324-2335).
    #[test]
    fn reversed_second_face_answers_c0() {
        let mut a_brep = BRep::new();
        let a_comp = make_coplanar_pair(&mut a_brep);
        let a_children = bat::sub_shapes(&a_comp);
        assert_eq!(a_children.len(), 2);
        let a_face1 = a_children[0].clone();
        let a_face2 = bat::oriented(&a_children[1], Orientation::Reversed);
        let a_shared = find_shared_edge(&a_brep, &a_comp);
        let a_cont = continuity_of_faces(&a_brep, &a_shared, &a_face1, &a_face2, CONFUSION);
        assert_eq!(a_cont, GeomAbsShape::C0);
    }

    /// The reader: after the coding, BRep_Tool::Continuity(E, F1, F2)
    /// answers CN (the elementary-surface coplanar outcome, see the walk
    /// test above).
    #[test]
    fn continuity_reader_answers_the_coded_record() {
        let mut a_brep = BRep::new();
        let a_comp = make_coplanar_pair(&mut a_brep);
        let a_children = bat::sub_shapes(&a_comp);
        let a_shared = find_shared_edge(&a_brep, &a_comp);
        encode_regularity_edge(&mut a_brep, &a_shared, &a_children[0], &a_children[1], CONFUSION);
        assert_eq!(
            brep_tool_continuity(&a_brep, &a_shared, &a_children[0], &a_children[1]),
            GeomAbsShape::CN
        );
    }
}
