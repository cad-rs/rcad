//! OCCT BRepTools_Modification (TKBRep/BRepTools — TKTopAlgo landing zone per
//! docs/module-map.md) — 1:1 translation of
//! BRepTools_Modification.hxx (L39-145) + BRepTools_Modification.cxx
//! (L24-38), together with the gp_Trsf subclass
//! BRepTools_TrsfModification (BRepTools_TrsfModification.hxx L35-148 +
//! BRepTools_TrsfModification.cxx L43-449).
//!
//! Source: $OCCT_SRC/src/ModelingData/TKBRep/BRepTools/
//!
//! The abstract interface is driven by BRepTools_Modifier — the sibling
//! module brep_tools_modifier.rs.
//!
//! Architecture differences:
//! 1. The OCCT virtual out-parameters are `occ::handle<Geom_Curve>` /
//!    `occ::handle<Geom_Surface>` / `occ::handle<Geom2d_Curve>` by reference;
//!    rcad carries the geometry by value, so the out-parameters are
//!    `&mut Option<Curve3>` / `&mut Option<Surface3>` / `&mut Option<Curve2d>`
//!    (the `None` value is the OCCT null handle).
//! 2. `TopLoc_Location` maps to the u32 location id (0 = identity) — the same
//!    carrier draft_modification.rs uses (its arch. diff. #4); the OCCT
//!    `L.Identity()` re-assignment of Draft_Modification.cxx L251/L283 is the
//!    id 0 in that carrier.
//! 3. `occ::handle<Poly_Triangulation>` maps to
//!    `rcad_kernel::math::poly::Triangulation` (the kernel's Poly_Triangulation
//!    translation).  `Poly_Polygon3D` / `Poly_PolygonOnTriangulation` have no
//!    rcad carrier at all (rcad BRep edges carry no polygon representation —
//!    see brep_check/brep_check_edge.rs L956-968), so those two
//!    out-parameters are dropped; the OCCT default bodies
//!    (BRepTools_Modification.cxx L29-38, `return false`) are preserved
//!    verbatim, which is the entire observable behaviour of the untranslated
//!    out-parameter.
//! 4. `TopoDS_Face` / `TopoDS_Edge` / `TopoDS_Vertex` arguments are `&Shape`.
//! 5. The virtuals are `&mut self` (the OCCT overrides are non-const and cache
//!    state — Draft_Modification.cxx L310-427 writes myCurves1d/myTol1d).

use glam::DVec3;

use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3};
use rcad_kernel::math::gp::Trsf;
use rcad_kernel::precision::{is_positive_infinite_value, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{GeomAbsShape, TShape};

use crate::brep_algo::tool::{
    brep_tool_curve_on_surface, brep_tool_parameter, brep_tool_pnt, brep_tool_range,
    brep_tool_tolerance, top_exp_vertices_raw,
};

// ---------------------------------------------------------------------------
// OCCT BRepTools_Modification (BRepTools_Modification.hxx L39-145).
// ---------------------------------------------------------------------------

/// OCCT BRepTools_Modification (hxx L41) — defines the geometric modifications
/// applied to a shape, i.e. the changes to its faces, edges and vertices.
pub trait BRepToolsModification {
    /// OCCT BRepTools_Modification::NewSurface (hxx L58-63, pure virtual).
    fn new_surface(
        &mut self,
        the_f: &Shape,
        the_s: &mut Option<Surface3>,
        the_l: &mut u32,
        the_tol: &mut f64,
        the_rev_wires: &mut bool,
        the_rev_face: &mut bool,
    ) -> bool;

    /// OCCT BRepTools_Modification::NewTriangulation (hxx L68-69) —
    /// BRepTools_Modification.cxx L24-27 default body `return false`.
    /// The OCCT out-parameter `Handle(Poly_Triangulation)` is the kernel
    /// Poly_Triangulation translation (arch. diff. #3).
    fn new_triangulation(
        &mut self,
        _the_f: &Shape,
        _the_t: &mut Option<rcad_kernel::math::poly::Triangulation>,
    ) -> bool {
        // OCCT L26: return false;
        false
    }

    /// OCCT BRepTools_Modification::NewCurve (hxx L78-81, pure virtual).
    fn new_curve(
        &mut self,
        the_e: &Shape,
        the_c: &mut Option<Curve3>,
        the_l: &mut u32,
        the_tol: &mut f64,
    ) -> bool;

    /// OCCT BRepTools_Modification::NewPolygon (hxx L86) —
    /// BRepTools_Modification.cxx L29-32 default body `return false`.
    /// The OCCT out-parameter `Handle(Poly_Polygon3D)` has no rcad carrier
    /// (arch. diff. #3), so it is dropped; the false result — the entire
    /// behaviour of the default body — is preserved.
    fn new_polygon(&mut self, _the_e: &Shape) -> bool {
        // OCCT L31: return false;
        false
    }

    /// OCCT BRepTools_Modification::NewPolygonOnTriangulation (hxx L91-94) —
    /// BRepTools_Modification.cxx L34-38 default body `return false`.
    /// The OCCT out-parameter `Handle(Poly_PolygonOnTriangulation)` has no
    /// rcad carrier (arch. diff. #3), so it is dropped.
    fn new_polygon_on_triangulation(&mut self, _the_e: &Shape, _the_f: &Shape) -> bool {
        // OCCT L37: return false;
        false
    }

    /// OCCT BRepTools_Modification::NewPoint (hxx L102, pure virtual).
    fn new_point(&mut self, the_v: &Shape, the_p: &mut DVec3, the_tol: &mut f64) -> bool;

    /// OCCT BRepTools_Modification::NewCurve2d (hxx L114-119, pure virtual).
    /// NewE is `&mut Shape` (the OCCT NewE handle is mutated in place by the
    /// BRep_Builder calls of the implementations — Draft_Modification
    /// arch. diff. #3).
    fn new_curve2d(
        &mut self,
        the_e: &Shape,
        the_f: &Shape,
        the_new_e: &mut Shape,
        the_new_f: &Shape,
        the_c: &mut Option<Curve2d>,
        the_tol: &mut f64,
    ) -> bool;

    /// OCCT BRepTools_Modification::NewParameter (hxx L127-130, pure virtual).
    fn new_parameter(
        &mut self,
        the_v: &Shape,
        the_e: &Shape,
        the_p: &mut f64,
        the_tol: &mut f64,
    ) -> bool;

    /// OCCT BRepTools_Modification::Continuity (hxx L137-142, pure virtual).
    fn continuity(
        &mut self,
        the_e: &Shape,
        the_f1: &Shape,
        the_f2: &Shape,
        the_new_e: &Shape,
        the_new_f1: &Shape,
        the_new_f2: &Shape,
    ) -> GeomAbsShape;
}

// ---------------------------------------------------------------------------
// OCCT BRepTools_TrsfModification
// (BRepTools_TrsfModification.hxx L35-148 + .cxx L43-449).
// ---------------------------------------------------------------------------

/// OCCT BRepTools_TrsfModification (hxx L38) — the gp_Trsf driven
/// modification: every function returns true and transforms the geometry.
pub struct BRepToolsTrsfModification {
    /// OCCT: myTrsf.
    my_trsf: Trsf,
    /// OCCT: myCopyMesh.
    my_copy_mesh: bool,
}

impl BRepToolsTrsfModification {
    /// OCCT BRepTools_TrsfModification::BRepTools_TrsfModification(T)
    /// (cxx L43-47).
    pub fn new(the_trsf: Trsf) -> Self {
        BRepToolsTrsfModification {
            my_trsf: the_trsf,
            my_copy_mesh: false,
        }
    }

    /// OCCT BRepTools_TrsfModification::Trsf() (cxx L51-54) — the OCCT
    /// accessor returns `gp_Trsf&`; the rcad read accessor returns `&Trsf`.
    pub fn trsf(&self) -> &Trsf {
        &self.my_trsf
    }

    /// OCCT BRepTools_TrsfModification::IsCopyMesh() (cxx L58-61) — the OCCT
    /// accessor returns `bool&`; the rcad mutable accessor is the same write
    /// path (the callers assign through it).
    pub fn is_copy_mesh_mut(&mut self) -> &mut bool {
        &mut self.my_copy_mesh
    }
}

impl BRepToolsModification for BRepToolsTrsfModification {
    /// OCCT BRepTools_TrsfModification::NewSurface (cxx L65-92).
    fn new_surface(
        &mut self,
        the_f: &Shape,
        the_s: &mut Option<Surface3>,
        the_l: &mut u32,
        the_tol: &mut f64,
        the_rev_wires: &mut bool,
        the_rev_face: &mut bool,
    ) -> bool {
        // OCCT L72: S = BRep_Tool::Surface(F, L);
        let (a_surf, a_loc) = brep_tool_surface_and_location(the_f);
        *the_s = a_surf;
        *the_l = a_loc;
        // OCCT L73-77: if (S.IsNull()) return false;
        let Some(a_surface) = the_s.clone() else {
            // processing cases when there is no geometry
            return false;
        };

        // OCCT L79: Tol = BRep_Tool::Tolerance(F);
        // OCCT L80: Tol *= std::abs(myTrsf.ScaleFactor());
        *the_tol = brep_tool_tolerance(the_f) * self.my_trsf.scale.abs();
        // OCCT L81-82.
        *the_rev_wires = false;
        *the_rev_face = self.my_trsf.is_negative();

        // OCCT L84-87: gp_Trsf LT = L.Transformation(); LT.Invert();
        // LT.Multiply(myTrsf); LT.Multiply(L.Transformation());
        let lt = location_transformation(*the_l)
            * self.my_trsf.to_daffine3()
            * location_transformation(*the_l);

        // OCCT L89: S = occ::down_cast<Geom_Surface>(S->Transformed(LT));
        *the_s = Some(rcad_kernel::geom::transform_surface(&a_surface, &lt));

        // OCCT L91: return true;
        true
    }

    /// OCCT BRepTools_TrsfModification::NewTriangulation (cxx L96-165).
    fn new_triangulation(
        &mut self,
        the_f: &Shape,
        the_t: &mut Option<rcad_kernel::math::poly::Triangulation>,
    ) -> bool {
        // OCCT L99-102: if (!myCopyMesh) return false;
        if !self.my_copy_mesh {
            return false;
        }
        let _ = (the_f, the_t);
        // OCCT L104-105: TopLoc_Location aLoc;
        // theTriangulation = BRep_Tool::Triangulation(theFace, aLoc);
        panic!(
            "GAP: BRep_Tool::Triangulation (BRep_Tool.cxx L1206-1230) not \
             translated — rcad BRep faces carry no Poly_Triangulation \
             (TFaceData has no triangulation slot)"
        );
    }

    /// OCCT BRepTools_TrsfModification::NewPolygon (cxx L169-218).
    fn new_polygon(&mut self, the_e: &Shape) -> bool {
        // OCCT L172-175: if (!myCopyMesh) return false;
        if !self.my_copy_mesh {
            return false;
        }
        let _ = the_e;
        // OCCT L177-178: TopLoc_Location aLoc;
        // theP = BRep_Tool::Polygon3D(theE, aLoc);
        panic!(
            "GAP: BRep_Tool::Polygon3D (BRep_Tool.cxx L1300-1330) not \
             translated — rcad BRep edges carry no Poly_Polygon3D \
             representation"
        );
    }

    /// OCCT BRepTools_TrsfModification::NewPolygonOnTriangulation
    /// (cxx L222-271).
    fn new_polygon_on_triangulation(&mut self, the_e: &Shape, the_f: &Shape) -> bool {
        // OCCT L227-230: if (!myCopyMesh) return false;
        if !self.my_copy_mesh {
            return false;
        }
        let _ = (the_e, the_f);
        // OCCT L232-233: TopLoc_Location aLoc;
        // occ::handle<Poly_Triangulation> aT = BRep_Tool::Triangulation(theF, aLoc);
        panic!(
            "GAP: BRep_Tool::Triangulation (BRep_Tool.cxx L1206-1230) not \
             translated — rcad BRep faces carry no Poly_Triangulation \
             (TFaceData has no triangulation slot)"
        );
    }

    /// OCCT BRepTools_TrsfModification::NewCurve (cxx L275-297).
    fn new_curve(
        &mut self,
        the_e: &Shape,
        the_c: &mut Option<Curve3>,
        the_l: &mut u32,
        the_tol: &mut f64,
    ) -> bool {
        // OCCT L280-281: double f, l;
        // C = BRep_Tool::Curve(E, L, f, l);
        let (a_curve, a_loc, _a_first, _a_last) = brep_tool_curve_and_location(the_e);
        *the_c = a_curve;
        *the_l = a_loc;

        // OCCT L283-284: Tol = BRep_Tool::Tolerance(E);
        // Tol *= std::abs(myTrsf.ScaleFactor());
        *the_tol = brep_tool_tolerance(the_e) * self.my_trsf.scale.abs();

        // OCCT L286-289: gp_Trsf LT = L.Transformation(); LT.Invert();
        // LT.Multiply(myTrsf); LT.Multiply(L.Transformation());
        let lt = location_transformation(*the_l)
            * self.my_trsf.to_daffine3()
            * location_transformation(*the_l);

        // OCCT L291-294: if (!C.IsNull()) C = down_cast(C->Transformed(LT));
        if let Some(a_c) = the_c.clone() {
            *the_c = Some(rcad_kernel::geom::transform_curve(&a_c, &lt));
        }

        // OCCT L296: return true;
        true
    }

    /// OCCT BRepTools_TrsfModification::NewPoint (cxx L301-309).
    fn new_point(&mut self, the_v: &Shape, the_p: &mut DVec3, the_tol: &mut f64) -> bool {
        // OCCT L303: P = BRep_Tool::Pnt(V);
        *the_p = brep_tool_pnt(the_v).unwrap_or(DVec3::ZERO);
        // OCCT L304-305: Tol = BRep_Tool::Tolerance(V);
        // Tol *= std::abs(myTrsf.ScaleFactor());
        *the_tol = brep_tool_tolerance(the_v) * self.my_trsf.scale.abs();
        // OCCT L306: P.Transform(myTrsf);
        *the_p = self.my_trsf.to_daffine3().transform_point3(*the_p);

        // OCCT L308: return true;
        true
    }

    /// OCCT BRepTools_TrsfModification::NewCurve2d (cxx L313-409).
    fn new_curve2d(
        &mut self,
        the_e: &Shape,
        the_f: &Shape,
        _the_new_e: &mut Shape,
        _the_new_f: &Shape,
        the_c: &mut Option<Curve2d>,
        the_tol: &mut f64,
    ) -> bool {
        // OCCT L320-321: TopLoc_Location loc;
        // OCCT L322: Tol = BRep_Tool::Tolerance(E);
        // OCCT L323: double scale = myTrsf.ScaleFactor();
        // OCCT L324: Tol *= std::abs(scale);
        let scale = self.my_trsf.scale;
        *the_tol = brep_tool_tolerance(the_e) * scale.abs();
        // OCCT L325: const occ::handle<Geom_Surface>& S = BRep_Tool::Surface(F, loc);
        let (a_surf, _a_loc) = brep_tool_surface_and_location(the_f);

        // OCCT L327-331: if (S.IsNull()) return false;
        let Some(a_surface) = a_surf else {
            // processing the case when the surface (geometry) is deleted
            return false;
        };
        // OCCT L332-335: GeomAdaptor_Surface GAsurf(S);
        // if (GAsurf.GetType() == GeomAbs_Plane) return false;
        let (a_basis, _a_bounds) = rcad_kernel::topods::surface_adaptor_basis_and_bounds(&a_surface);
        if matches!(a_basis, Surface3::Plane(_)) {
            return false;
        }

        // OCCT L337-342: double f, l;
        // occ::handle<Geom2d_Curve> NewC = BRep_Tool::CurveOnSurface(E, F, f, l);
        // if (NewC.IsNull()) return false;
        let Some((mut new_c, mut f, mut l)) = brep_tool_curve_on_surface(the_e, the_f) else {
            return false;
        };

        // OCCT L344-352: the Geom2d_TrimmedCurve basis unwrap.
        if let Curve2d::Trimmed(tc) = &new_c {
            new_c = (*tc.curve).clone();
        }

        // OCCT L354: double fc = NewC->FirstParameter(), lc = NewC->LastParameter();
        let (fc, lc) = curve2d_range(&new_c);

        // OCCT L356-377: if (!NewC->IsPeriodic()) { the range readjustment }.
        if !new_c.is_periodic() {
            if fc - f > PCONFUSION {
                f = fc;
            }
            if l - lc > PCONFUSION {
                l = lc;
            }
            if (l - f).abs() < PCONFUSION {
                if (f - fc).abs() < PCONFUSION && !is_positive_infinite_value(lc) {
                    l = lc;
                } else if !is_positive_infinite_value(fc) && !is_negative_infinite(fc) {
                    f = fc;
                }
            }
        }

        // OCCT L379-380: newf = f; newl = l;
        // (the OCCT reassignment of newf / newl lives in the L381-397 branch
        // below, whose GTransform step is untranslated.)
        let new_f = f;
        let new_l = l;
        // OCCT L381-397: if (std::abs(scale) != 1.) { the GTransform branch }.
        if scale.abs() != 1.0 {
            // OCCT L384: NewC = new Geom2d_TrimmedCurve(NewC, f, l);
            new_c = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                curve: Box::new(new_c),
                t_min: f,
                t_max: l,
            });
            // OCCT L385: gp_GTrsf2d gtrsf = S->ParametricTransformation(myTrsf);
            let _ = &new_c;
            panic!(
                "GAP: Geom_Surface::ParametricTransformation (Geom_Surface.cxx \
                 L300-330) and GeomLib::GTransform (GeomLib.cxx L1001-1100) \
                 not translated — the BRepTools_TrsfModification::NewCurve2d \
                 L381-397 branch is unreachable in rcad"
            );
        }

        // OCCT L398-405: the 3D / 2D range readjustment.
        // TopoDS_Vertex V1, V2; TopExp::Vertices(E, V1, V2);
        let (a_v1, a_v2) = top_exp_vertices_raw(the_e);
        let a_v1 = a_v1.unwrap_or_else(Shape::null);
        let a_v2 = a_v2.unwrap_or_else(Shape::null);
        // TopoDS_Edge EFOR = TopoDS::Edge(E.Oriented(TopAbs_FORWARD));
        let a_efor = {
            let mut s = the_e.clone();
            s.orientation = rcad_kernel::topods::Orientation::Forward;
            s
        };
        // double aTolV; NewParameter(V1, EFOR, f, aTolV); NewParameter(V2, EFOR, l, aTolV);
        let mut a_tol_v = 0.0;
        self.new_parameter(&a_v1, &a_efor, &mut f, &mut a_tol_v);
        self.new_parameter(&a_v2, &a_efor, &mut l, &mut a_tol_v);
        // OCCT L406: GeomLib::SameRange(Precision::PConfusion(), NewC, newf,
        // newl, f, l, C);
        *the_c = Some(rcad_kernel::geom::same_range_2d(
            PCONFUSION, new_c, new_f, new_l, f, l,
        )
        .expect("Standard_NullObject: GeomLib::SameRange (GeomLib.cxx L842-970)"));

        // OCCT L408: return true;
        true
    }

    /// OCCT BRepTools_TrsfModification::NewParameter (cxx L413-437).
    fn new_parameter(
        &mut self,
        the_v: &Shape,
        the_e: &Shape,
        the_p: &mut f64,
        the_tol: &mut f64,
    ) -> bool {
        // OCCT L418-421: if (V.IsNull()) return false; // infinite edge may
        // have Null vertex
        if the_v.is_null() {
            return false;
        }

        // OCCT L423-426: TopLoc_Location loc;
        // Tol = BRep_Tool::Tolerance(V); Tol *= std::abs(myTrsf.ScaleFactor());
        // P = BRep_Tool::Parameter(V, E);
        *the_tol = brep_tool_tolerance(the_v) * self.my_trsf.scale.abs();
        *the_p = brep_tool_parameter(the_v, the_e);

        // OCCT L428-434: double f, l;
        // occ::handle<Geom_Curve> C = BRep_Tool::Curve(E, loc, f, l);
        // if (!C.IsNull()) P = C->TransformedParameter(P, myTrsf);
        let (a_curve, _a_loc, _a_f, _a_l) = brep_tool_curve_and_location(the_e);
        if let Some(a_c) = a_curve {
            *the_p = a_c.transformed_parameter(*the_p);
        }

        // OCCT L436: return true;
        true
    }

    /// OCCT BRepTools_TrsfModification::Continuity (cxx L441-449).
    fn continuity(
        &mut self,
        the_e: &Shape,
        the_f1: &Shape,
        the_f2: &Shape,
        _the_new_e: &Shape,
        _the_new_f1: &Shape,
        _the_new_f2: &Shape,
    ) -> GeomAbsShape {
        // OCCT L448: return BRep_Tool::Continuity(E, F1, F2);
        crate::offset::draft_modification::brep_tool_continuity(the_e, the_f1, the_f2)
    }
}

// ---------------------------------------------------------------------------
// Re-hosts of the TKBRep readers the translation needs.
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Surface(F, L) (BRep_Tool.cxx L594-620) — the face surface
/// in its LOCAL frame plus the location that places it in world coordinates
/// (rcad: the face TShape's stored surface and its location id; the OCCT
/// location is the face's own `TopLoc_Location`).
fn brep_tool_surface_and_location(the_f: &Shape) -> (Option<Surface3>, u32) {
    match the_f.data.as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fd.surface_location),
        _ => (None, 0),
    }
}

/// OCCT BRep_Tool::Curve(E, L, f, l) (BRep_Tool.cxx L184-215) — the 3D curve
/// with its location and parameter range.
fn brep_tool_curve_and_location(the_e: &Shape) -> (Option<Curve3>, u32, f64, f64) {
    let curve = crate::brep_algo::tool::brep_tool_curve(the_e).map(|(c, _, _)| c);
    let (f, l) = brep_tool_range(the_e);
    (curve, the_e.location, f, l)
}

/// OCCT Geom2d_Curve::FirstParameter() / LastParameter() (Geom2d_Curve.hxx
/// L120-137) — the natural parameter domain of the 2d curve (the kernel
/// Curve2dEval::default_domain).
fn curve2d_range(the_c: &Curve2d) -> (f64, f64) {
    let d = the_c.default_domain();
    (d[0], d[1])
}

/// OCCT TopLoc_Location::Transformation() (TopLoc_Location.hxx L139-142) for
/// the rcad u32 location id.  The table an id indexes is the owning pool's
/// `BRep::locations`; the modification carrier holds no table, and the rcad
/// modification implementors carry the identity location
/// (Draft_Modification.cxx L251/L283 `L.Identity()` is the id 0 in this
/// carrier — draft_modification.rs arch. diff. #4), for which the OCCT
/// transformation is the identity matrix.
fn location_transformation(_the_loc: u32) -> glam::DAffine3 {
    glam::DAffine3::IDENTITY
}

/// OCCT Precision::IsNegativeInfinite (Precision.hxx L118-124) applied to the
/// doubled value (the gp 1d value lives in the rcad f64 parameter space).
fn is_negative_infinite(r: f64) -> bool {
    r <= -rcad_kernel::precision::REAL_LAST
}
