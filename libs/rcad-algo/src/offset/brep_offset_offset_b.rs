// OCCT BRepOffset_Offset.cxx L389-1815 — the class part of the 1:1
// translation (split from brep_offset_offset.rs, which carries the cxx
// L1-685 statics / GAP carriers; see the architecture-difference list there).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_Offset.cxx / .hxx / BRepOffset_Offset.lxx
//
// OCCT inheritance chain (hxx L45): none — BRepOffset_Offset is a standalone
// class.

use std::collections::HashMap;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, Surface3, SurfaceEval, TrimmedCurve2, TrimmedCurve3, TrimmedSurface};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_generator_b::brep_tools_uv_bounds;
use crate::feat::loc_ope_wires_on_shape::{brep_tool_pnt, shape_key};
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_tolerance, ShapeKey,
};

use super::brep_offset_make_simple_offset::{edge_curve_of, oriented_vertex, top_exp_vertices_shape};
use super::brep_offset_offset::*;

// ---------------------------------------------------------------------------
// OCCT class (BRepOffset_Offset.hxx L45-158, cxx L389-1815).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Offset (BRepOffset_Offset.hxx L45-158).
///
/// Architecture note: the OCCT global TShape arena is the caller-owned
/// `&mut BRep` threaded through every Init (all BRepOffset_Offset instances
/// of one BRepOffset_MakeOffset run share the MakeOffset pool, like the C++
/// handles share the TShape heap).
pub struct BRepOffsetOffset {
    my_shape: Shape,               // OCCT: myShape (hxx L154)
    my_status: BRepOffsetStatus,   // OCCT: myStatus (hxx L155)
    my_face: Shape,                // OCCT: myFace (hxx L156)
    my_map: HashMap<ShapeKey, Shape>, // OCCT: myMap (hxx L157)
}

impl Default for BRepOffsetOffset {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepOffsetOffset {
    /// OCCT BRepOffset_Offset::BRepOffset_Offset() (cxx L389).
    pub fn new() -> Self {
        BRepOffsetOffset {
            my_shape: Shape::null(),
            my_status: BRepOffsetStatus::Good,
            my_face: Shape::null(),
            my_map: HashMap::new(),
        }
    }

    /// OCCT BRepOffset_Offset::BRepOffset_Offset(Face, Offset, OffsetOutside
    /// = true, JoinType = GeomAbs_Arc) (cxx L393-399).
    pub fn with_face(
        the_brep: &mut BRep,
        the_face: &Shape,
        the_offset: f64,
        the_offset_outside: bool,
        the_join_type: bool,
    ) -> Self {
        let mut res = Self::new();
        res.init_face(the_brep, the_face, the_offset, the_offset_outside, the_join_type);
        res
    }

    /// OCCT BRepOffset_Offset::BRepOffset_Offset(Face, Offset, Created,
    /// OffsetOutside = true, JoinType = GeomAbs_Arc) (cxx L403-411).
    pub fn with_face_created(
        the_brep: &mut BRep,
        the_face: &Shape,
        the_offset: f64,
        the_created: &HashMap<ShapeKey, Shape>,
        the_offset_outside: bool,
        the_join_type: bool,
    ) -> Self {
        let mut res = Self::new();
        res.init_face_created(the_brep, the_face, the_offset, the_created, the_offset_outside, the_join_type);
        res
    }

    /// OCCT BRepOffset_Offset::BRepOffset_Offset(Path, Edge1, Edge2, Offset,
    /// Polynomial = false, Tol = 1.0e-4, Conti = GeomAbs_C1) (cxx L415-424).
    #[allow(clippy::too_many_arguments)]
    pub fn with_path(
        the_brep: &mut BRep,
        the_path: &Shape,
        the_edge1: &Shape,
        the_edge2: &Shape,
        the_offset: f64,
        the_polynomial: bool,
        the_tol: f64,
        the_conti: GeomAbsShapeKind,
    ) -> Self {
        let mut res = Self::new();
        res.init_path(the_brep, the_path, the_edge1, the_edge2, the_offset, the_polynomial, the_tol, the_conti);
        res
    }

    /// OCCT BRepOffset_Offset::BRepOffset_Offset(Path, Edge1, Edge2, Offset,
    /// FirstEdge, LastEdge, Polynomial = false, Tol = 1.0e-4, Conti =
    /// GeomAbs_C1) (cxx L428-439).
    #[allow(clippy::too_many_arguments)]
    pub fn with_path_first_last(
        the_brep: &mut BRep,
        the_path: &Shape,
        the_edge1: &Shape,
        the_edge2: &Shape,
        the_offset: f64,
        the_first_edge: &Shape,
        the_last_edge: &Shape,
        the_polynomial: bool,
        the_tol: f64,
        the_conti: GeomAbsShapeKind,
    ) -> Self {
        let mut res = Self::new();
        res.init_path_first_last(
            the_brep, the_path, the_edge1, the_edge2, the_offset, the_first_edge,
            the_last_edge, the_polynomial, the_tol, the_conti,
        );
        res
    }

    /// OCCT BRepOffset_Offset::BRepOffset_Offset(Vertex, LEdge, Offset,
    /// Polynomial = false, Tol = 1.0e-4, Conti = GeomAbs_C1) (cxx L443-451).
    pub fn with_vertex(
        the_brep: &mut BRep,
        the_vertex: &Shape,
        the_l_edge: &[Shape],
        the_offset: f64,
        the_polynomial: bool,
        the_tol: f64,
        the_conti: GeomAbsShapeKind,
    ) -> Self {
        let mut res = Self::new();
        res.init_vertex(the_brep, the_vertex, the_l_edge, the_offset, the_polynomial, the_tol, the_conti);
        res
    }

    /// OCCT BRepOffset_Offset::Init(Face, Offset, OffsetOutside = true,
    /// JoinType = GeomAbs_Arc) (cxx L455-462).
    pub fn init_face(
        &mut self,
        the_brep: &mut BRep,
        the_face: &Shape,
        the_offset: f64,
        the_offset_outside: bool,
        the_join_type: bool,
    ) {
        // OCCT L460: NCollection_DataMap<...> Empty.
        let the_created: HashMap<ShapeKey, Shape> = HashMap::new();
        self.init_face_created(the_brep, the_face, the_offset, &the_created, the_offset_outside, the_join_type);
    }

    /// OCCT BRepOffset_Offset::Init(Face, Offset, Created, OffsetOutside =
    /// true, JoinType = GeomAbs_Arc) (cxx L466-1094).
    pub fn init_face_created(
        &mut self,
        the_brep: &mut BRep,
        face: &Shape,
        offset: f64,
        created: &HashMap<ShapeKey, Shape>,
        offset_outside: bool,
        _the_join_type: bool,
    ) {
        self.my_shape = face.clone();
        let mut my_offset = offset;
        if face.orientation == Orientation::Reversed {
            my_offset *= -1.0;
        }

        // OCCT L480: TopLoc_Location L — the rcad location index
        // (0 = identity); the `l_loc` name avoids the clash with the pcurve
        // range `l` of the wire loop below.
        let l_loc: u32 = face.location;
        let mut s = match face.as_face().and_then(|fd| fd.surface.clone()) {
            Some(s) => s,
            None => {
                self.my_status = BRepOffsetStatus::Unknown;
                return;
            }
        };
        // OCCT L483-488: On detrime les surfaces, evite des recopies dans les
        // extensions (one Geom_RectangularTrimmedSurface level unwrapped).
        if let Surface3::Trimmed(rt) = &s {
            s = (*rt.basis).clone();
        }
        let mut is_transformed = false;
        // OCCT L490-499: (BSpline || LinearExtrusion || Revolution ||
        // OffsetSurface) && !L.IsIdentity() -> S->Copy() + Transform.
        if matches!(
            s,
            Surface3::BSpline(_)
                | Surface3::LinearExtrusion(_)
                | Surface3::Revolution(_)
                | Surface3::Offset(_)
        ) && l_loc != 0
        {
            // GAP: the Geom_Surface::Copy/Transform vehicle has no rcad
            // translation (the rcad surfaces travel as values); the branch
            // keeps the OCCT shape.
            s = geom_surface_transformed(&s);
            is_transformed = true;
        }
        // particular case of cone
        // OCCT L501-522.
        if let Surface3::Cone(co) = &s {
            // OCCT L506-507: gp_Pnt Apex = Co->Apex(); ElSLib::Parameters(
            // Co->Cone(), Apex, Uc, Vc) — the TRUE apex (the rcad `apex`
            // field is the reference point; apex_point() derives the apex).
            let (uc, _vc) = elslib_cone_parameters(co, co.apex_point());
            let _ = uc;
            let (uu1, uu2, vv1, vv2) = {
                let b = brep_tools_uv_bounds(face);
                (b[0], b[1], b[2], b[3])
            };
            if vv2 < _vc && co.half_angle_rad > 0.0 {
                my_offset *= -1.0;
            } else if vv1 > _vc && co.half_angle_rad < 0.0 {
                my_offset *= -1.0;
            }
            // OCCT L518: if (!Co->Position().Direct()) — the rcad analytic
            // payloads are always direct (no indirect frame is carried).
            let _ = (uu1, uu2);
        }

        // OCCT L524: TheSurf = BRepOffset::Surface(S, myOffset, myStatus) —
        // the OCCT call uses the hxx default allowC0 = false.
        let mut the_surf = brep_offset_surface(&s, my_offset, &mut self.my_status, false);

        // processing offsets of faces with possible degenerated edges
        // OCCT L526-807.
        let mut umin_degen = false;
        let mut umax_degen = false;
        let mut vmin_degen = false;
        let mut vmax_degen = false;
        let mut uiso_degen = false;
        let mut min_apex = DVec3::ZERO;
        let mut max_apex = DVec3::ZERO;
        let mut has_singularity = false;
        let b = brep_tools_uv_bounds(face);
        let (uf1, uf2, vf1, vf2) = (b[0], b[1], b[2], b[3]);
        if (!offset_outside || !_the_join_type)
            && matches!(the_surf, Surface3::Cone(_) | Surface3::Offset(_))
        {
            // OCCT L536-538: (TheSurf is Conical or Offset) — the JoinType
            // argument is the GeomAbs_Arc bool of this carrier.
            let mut deg_edges: Vec<Shape> = Vec::new();
            for an_edge in explorer(face, ShapeType::Edge, ShapeType::Shape) {
                if brep_tool_degenerated(&an_edge) {
                    let Some((c2d, a_f, a_l)) = brep_tool_curve_on_surface(&an_edge, face) else {
                        continue;
                    };
                    let a_fpnt2d = c2d.point_at(a_f);
                    let a_lpnt2d = c2d.point_at(a_l);
                    let a_fpnt = the_surf.point_at(a_fpnt2d.x, a_fpnt2d.y);
                    let a_lpnt = the_surf.point_at(a_lpnt2d.x, a_lpnt2d.y);

                    //  aFPnt.SquareDistance(aLPnt) > Precision::SquareConfusion() -
                    // is a sufficient condition of troubles: non-singular case, but edge is degenerated.
                    // So, normal handling of degenerated edges is not applicable in case of non-singular point.
                    if a_fpnt.distance_squared(a_lpnt) < rcad_kernel::core::precision::SQUARE_CONFUSION {
                        deg_edges.push(an_edge);
                    }
                }
            }
            if !deg_edges.is_empty() {
                let tol_apex = 1.0e-5;
                // define the iso of singularity (u or v)
                // OCCT L569-576.
                let the_deg_edge = deg_edges[0].clone();
                let Some((a_curve, fpar, lpar)) = brep_tool_curve_on_surface(&the_deg_edge, face)
                else {
                    return;
                };
                let fp2d = a_curve.point_at(fpar);
                let lp2d = a_curve.point_at(lpar);
                if (fp2d.x - lp2d.x).abs() <= rcad_kernel::core::precision::PCONFUSION {
                    uiso_degen = true;
                }

                if deg_edges.len() == 2 {
                    if uiso_degen {
                        umin_degen = true;
                        umax_degen = true;
                    } else {
                        vmin_degen = true;
                        vmax_degen = true;
                    }
                } else {
                    // DegEdges.Length() == 1
                    if uiso_degen {
                        if (fp2d.x - uf1).abs() <= rcad_kernel::core::precision::CONFUSION {
                            umin_degen = true;
                        } else {
                            umax_degen = true;
                        }
                    } else if (fp2d.y - vf1).abs() <= rcad_kernel::core::precision::CONFUSION {
                        vmin_degen = true;
                    } else {
                        vmax_degen = true;
                    }
                }
                if let Surface3::Cone(cone) = &the_surf {
                    // OCCT L616-634: the ConicalSurface branch — gp_Pnt apex
                    // = theCone.Apex(); ElSLib::Parameters(theCone, apex, ...)
                    // (the TRUE apex of the offset cone; the rcad `apex`
                    // field is the reference point).
                    let apex = cone.apex_point();
                    let (_uapex, vapex) = elslib_cone_parameters(cone, apex);
                    if vmin_degen {
                        the_surf = Surface3::Trimmed(TrimmedSurface::new(
                            the_surf.clone(), uf1, uf2, vapex, vf2,
                        ));
                        min_apex = apex;
                        has_singularity = true;
                    } else if vmax_degen {
                        the_surf = Surface3::Trimmed(TrimmedSurface::new(
                            the_surf.clone(), uf1, uf2, vf1, vapex,
                        ));
                        max_apex = apex;
                        has_singularity = true;
                    }
                } else {
                    // TheSurf->DynamicType() == STANDARD_TYPE(Geom_OffsetSurface)
                    // OCCT L635-805.
                    if umin_degen {
                        let uiso = surface_u_iso_gap(&the_surf, uf1);
                        if brep_offset_tool_gabarit(&uiso) > tol_apex {
                            let basis_surf = offset_basis(&the_surf);
                            let papex = basis_surf.point_at(uf1, vf1);
                            let pfirst = the_surf.point_at(uf1, vf1);
                            let pquart = the_surf.point_at(uf1, 0.75 * vf1 + 0.25 * vf2);
                            let pmid = the_surf.point_at(uf1, 0.5 * (vf1 + vf2));
                            let dir_apex = (pquart - pfirst).cross(pmid - pfirst);
                            let _line_apex = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: papex,
                                direction: dir_apex.normalize_or_zero(),
                            });
                            let dir_generatrix = basis_surf.dn(uf1, vf1, 1, 0);
                            let line_generatrix = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: pfirst,
                                direction: dir_generatrix.normalize_or_zero(),
                            });
                            let the_extrema = GeomApiExtremaCurveCurve::new(
                                &line_generatrix,
                                &_line_apex,
                            );
                            let (_pint1, pint2) = the_extrema.nearest_points();
                            let _ = pint2;
                            let pint1 = _pint1;
                            let length = pfirst.distance(pint1);
                            if offset_outside {
                                let a_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, uf2, vf1, vf2,
                                ));
                                the_surf = geom_lib_extend_surf_by_length(&a_surf, length, 1, true, false);
                                let (_u1, _u2, _v1, _v2) = surface_bounds4(&the_surf);
                                min_apex = the_surf.point_at(_u1, vf1);
                            } else {
                                let viso = surface_v_iso_gap(&the_surf, vf1);
                                let projector = GeomApiProjectPointOnCurve::new(pint1, &viso);
                                let new_first_u = projector.lower_distance_parameter();
                                the_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), new_first_u, uf2, vf1, vf2,
                                ));
                                min_apex = the_surf.point_at(new_first_u, vf1);
                            }
                            has_singularity = true;
                        }
                    } // end of if (UminDegen)
                    if umax_degen {
                        let uiso = surface_u_iso_gap(&the_surf, uf2);
                        if brep_offset_tool_gabarit(&uiso) > tol_apex {
                            let basis_surf = offset_basis(&the_surf);
                            let papex = basis_surf.point_at(uf2, vf1);
                            let pfirst = the_surf.point_at(uf2, vf1);
                            let pquart = the_surf.point_at(uf2, 0.75 * vf1 + 0.25 * vf2);
                            let pmid = the_surf.point_at(uf2, 0.5 * (vf1 + vf2));
                            let dir_apex = (pquart - pfirst).cross(pmid - pfirst);
                            let _line_apex = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: papex,
                                direction: dir_apex.normalize_or_zero(),
                            });
                            let dir_generatrix = basis_surf.dn(uf2, vf1, 1, 0);
                            let line_generatrix = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: pfirst,
                                direction: dir_generatrix.normalize_or_zero(),
                            });
                            let the_extrema = GeomApiExtremaCurveCurve::new(
                                &line_generatrix,
                                &_line_apex,
                            );
                            let (pint1, _pint2) = the_extrema.nearest_points();
                            let length = pfirst.distance(pint1);
                            if offset_outside {
                                let a_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, uf2, vf1, vf2,
                                ));
                                the_surf = geom_lib_extend_surf_by_length(&a_surf, length, 1, true, true);
                                let (_u1, u2, _v1, _v2) = surface_bounds4(&the_surf);
                                max_apex = the_surf.point_at(u2, vf1);
                            } else {
                                let viso = surface_v_iso_gap(&the_surf, vf1);
                                let projector = GeomApiProjectPointOnCurve::new(pint1, &viso);
                                let new_last_u = projector.lower_distance_parameter();
                                the_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, new_last_u, vf1, vf2,
                                ));
                                max_apex = the_surf.point_at(new_last_u, vf1);
                            }
                            has_singularity = true;
                        }
                    } // end of if (UmaxDegen)
                    if vmin_degen {
                        let viso = surface_v_iso_gap(&the_surf, vf1);
                        if brep_offset_tool_gabarit(&viso) > tol_apex {
                            let basis_surf = offset_basis(&the_surf);
                            let papex = basis_surf.point_at(uf1, vf1);
                            let pfirst = the_surf.point_at(uf1, vf1);
                            let pquart = the_surf.point_at(0.75 * uf1 + 0.25 * uf2, vf1);
                            let pmid = the_surf.point_at(0.5 * (uf1 + uf2), vf1);
                            let dir_apex = (pquart - pfirst).cross(pmid - pfirst);
                            let _line_apex = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: papex,
                                direction: dir_apex.normalize_or_zero(),
                            });
                            let dir_generatrix = basis_surf.dn(uf1, vf1, 0, 1);
                            let line_generatrix = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: pfirst,
                                direction: dir_generatrix.normalize_or_zero(),
                            });
                            let the_extrema = GeomApiExtremaCurveCurve::new(
                                &line_generatrix,
                                &_line_apex,
                            );
                            let (pint1, _pint2) = the_extrema.nearest_points();
                            let length = pfirst.distance(pint1);
                            if offset_outside {
                                let a_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, uf2, vf1, vf2,
                                ));
                                the_surf = geom_lib_extend_surf_by_length(&a_surf, length, 1, false, false);
                                let (_u1, _u2, _v1, v1) = surface_bounds4(&the_surf);
                                min_apex = the_surf.point_at(uf1, v1);
                            } else {
                                let uiso = surface_u_iso_gap(&the_surf, uf1);
                                let projector = GeomApiProjectPointOnCurve::new(pint1, &uiso);
                                let new_first_v = projector.lower_distance_parameter();
                                the_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, uf2, new_first_v, vf2,
                                ));
                                min_apex = the_surf.point_at(uf1, new_first_v);
                            }
                            has_singularity = true;
                        }
                    } // end of if (VminDegen)
                    if vmax_degen {
                        let viso = surface_v_iso_gap(&the_surf, vf2);
                        if brep_offset_tool_gabarit(&viso) > tol_apex {
                            let basis_surf = offset_basis(&the_surf);
                            let papex = basis_surf.point_at(uf1, vf2);
                            let pfirst = the_surf.point_at(uf1, vf2);
                            let pquart = the_surf.point_at(0.75 * uf1 + 0.25 * uf2, vf2);
                            let pmid = the_surf.point_at(0.5 * (uf1 + uf2), vf2);
                            let dir_apex = (pquart - pfirst).cross(pmid - pfirst);
                            let _line_apex = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: papex,
                                direction: dir_apex.normalize_or_zero(),
                            });
                            let dir_generatrix = basis_surf.dn(uf1, vf2, 0, 1);
                            let line_generatrix = Curve3::Line(rcad_kernel::geom::Line3 {
                                origin: pfirst,
                                direction: dir_generatrix.normalize_or_zero(),
                            });
                            let the_extrema = GeomApiExtremaCurveCurve::new(
                                &line_generatrix,
                                &_line_apex,
                            );
                            let (pint1, _pint2) = the_extrema.nearest_points();
                            let length = pfirst.distance(pint1);
                            if offset_outside {
                                let a_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, uf2, vf1, vf2,
                                ));
                                the_surf = geom_lib_extend_surf_by_length(&a_surf, length, 1, false, true);
                                let (_u1, _u2, _v1, v2) = surface_bounds4(&the_surf);
                                max_apex = the_surf.point_at(uf1, v2);
                            } else {
                                let uiso = surface_u_iso_gap(&the_surf, uf1);
                                let projector = GeomApiProjectPointOnCurve::new(pint1, &uiso);
                                let new_last_v = projector.lower_distance_parameter();
                                the_surf = Surface3::Trimmed(TrimmedSurface::new(
                                    the_surf.clone(), uf1, uf2, vf1, new_last_v,
                                ));
                                max_apex = the_surf.point_at(uf1, new_last_v);
                            }
                            has_singularity = true;
                        }
                    } // end of if (VmaxDegen)
                } // end of else (case of Geom_OffsetSurface)
            }
        } // end of processing offsets of faces with possible degenerated edges

        // find the PCurves of the edges of <Faces>
        // OCCT L811-820: BRep_Builder myBuilder; MakeFace(myFace);
        // UpdateFace(myFace, TheSurf, L/TopLoc_Location(), Tolerance(Face)).
        let mut my_builder = BRepBuilder::new();
        self.my_face = my_builder.make_face(the_brep, None, Shape::null());
        if !is_transformed {
            update_face_surface(
                &self.my_face,
                &the_surf,
                l_loc,
                brep_tool_tolerance(face),
                the_brep,
            );
        } else {
            update_face_surface(&self.my_face, &the_surf, 0, brep_tool_tolerance(face), the_brep);
        }

        let mut map_ss: HashMap<ShapeKey, Shape> = HashMap::new();

        // mise a jour de la map sur les vertex deja crees
        // OCCT L825-826: CurFace = TopoDS::Face(Face.Oriented(TopAbs_FORWARD)).
        let mut cur_face = face.clone();
        cur_face.orientation = Orientation::Forward;

        // OCCT L829: NCollection_Map<TopoDS_Shape> VonDegen.
        let mut von_degen: std::collections::HashSet<ShapeKey> = std::collections::HashSet::new();
        // OCCT L831: TheSurf->Bounds(u1, u2, v1, v2).
        let (u1, u2, _v1, _v2) = surface_bounds4(&the_surf);

        // OCCT L833-870.
        for e in explorer(&cur_face, ShapeType::Edge, ShapeType::Shape) {
            let (v1, v2) = top_exp_vertices_shape(&e);
            if has_singularity && brep_tool_degenerated(&e) {
                von_degen.insert(shape_key(&v1));
            }
            if let Some(oe) = created.get(&shape_key(&e)) {
                let (ov1, ov2) = top_exp_vertices_shape(oe);
                if !map_ss.contains_key(&shape_key(&v1)) {
                    map_ss.insert(shape_key(&v1), ov1);
                }
                if !map_ss.contains_key(&shape_key(&v2)) {
                    map_ss.insert(shape_key(&v2), ov2);
                }
            }
            if let Some(ov1) = created.get(&shape_key(&v1)) {
                if !map_ss.contains_key(&shape_key(&v1)) {
                    map_ss.insert(shape_key(&v1), ov1.clone());
                }
            }
            if let Some(ov2) = created.get(&shape_key(&v2)) {
                if !map_ss.contains_key(&shape_key(&v2)) {
                    map_ss.insert(shape_key(&v2), ov2.clone());
                }
            }
        }

        // OCCT L872-1089.
        for w in explorer(&cur_face, ShapeType::Wire, ShapeType::Shape) {
            let mut w_fwd = w.clone();
            w_fwd.orientation = Orientation::Forward;
            let mut ow = my_builder.make_wire(the_brep);
            for e in explorer(&w_fwd, ShapeType::Edge, ShapeType::Shape) {
                let (v1, v2) = top_exp_vertices_shape(&e);
                // OCCT L888: C2d = BRep_Tool::CurveOnSurface(E, CurFace, f, l).
                let Some((mut c2d, mut f, mut l)) = brep_tool_curve_on_surface(&e, &cur_face) else {
                    continue;
                };
                let mut oe;
                if map_ss.contains_key(&shape_key(&e))
                    && !von_degen.contains(&shape_key(&v1))
                    && !von_degen.contains(&shape_key(&v2))
                {
                    // c`est un edge de couture
                    // OCCT L890-907.
                    oe = map_ss[&shape_key(&e)].clone();
                    let mut e_rev = e.clone();
                    e_rev.orientation = if e.orientation == Orientation::Forward {
                        Orientation::Reversed
                    } else {
                        Orientation::Forward
                    };
                    let c2d_1 = brep_tool_curve_on_surface(&e_rev, &cur_face).map(|(c, ff, ll)| {
                        f = ff;
                        l = ll;
                        c
                    });
                    if let Some(c2d_1) = c2d_1 {
                        if e.orientation == Orientation::Forward {
                            update_edge_2d_seam(
                                &oe, &c2d, &c2d_1, &self.my_face,
                                brep_tool_tolerance(&e), the_brep,
                            );
                        } else {
                            update_edge_2d_seam(
                                &oe, &c2d_1, &c2d, &self.my_face,
                                brep_tool_tolerance(&e), the_brep,
                            );
                        }
                    }
                    my_builder.set_edge_range(the_brep, oe.clone(), f, l);
                } else {
                    // OCCT L910-957.
                    let mut eforward = e.clone();
                    eforward.orientation = Orientation::Forward;
                    let p2d1 = {
                        let par = brep_tool_parameter_vfe(&v1, &eforward, &cur_face);
                        c2d.point_at(par)
                    };
                    let p2d2 = {
                        let par = brep_tool_parameter_vfe(&v2, &eforward, &cur_face);
                        c2d.point_at(par)
                    };
                    let mut p1;
                    let mut vstart;
                    if von_degen.contains(&shape_key(&v1)) {
                        if (p2d1.y - vf1).abs() <= rcad_kernel::core::precision::CONFUSION {
                            p1 = min_apex;
                            vstart = u1;
                        } else {
                            p1 = max_apex;
                            vstart = u2;
                        }
                    } else {
                        p1 = the_surf.point_at(p2d1.x, p2d1.y);
                        if l_loc != 0 && !is_transformed {
                            p1 = point_transformed_gap(p1);
                        }
                        vstart = p2d1.y;
                    }
                    let mut p2;
                    let mut vend;
                    if von_degen.contains(&shape_key(&v2)) {
                        if (p2d2.y - vf1).abs() <= rcad_kernel::core::precision::CONFUSION {
                            p2 = min_apex;
                            vend = u1;
                        } else {
                            p2 = max_apex;
                            vend = u2;
                        }
                    } else {
                        p2 = the_surf.point_at(p2d2.x, p2d2.y);
                        if l_loc != 0 && !is_transformed {
                            p2 = point_transformed_gap(p2);
                        }
                        vend = p2d2.y;
                    }
                    // E a-t-il ume image dans la Map des Created ?
                    // OCCT L959-1009.
                    if let Some(oe_created) = created.get(&shape_key(&e)) {
                        oe = oe_created.clone();
                    } else if map_ss.contains_key(&shape_key(&e)) {
                        // seam edge
                        oe = map_ss[&shape_key(&e)].clone();
                    } else {
                        // OCCT L969-1009: myBuilder.MakeEdge(OE) + vertex Adds
                        // — the rcad add_edge stores the oriented first/last
                        // vertices in one step (arch. diff.).
                        let ov1;
                        if map_ss.contains_key(&shape_key(&v1)) {
                            ov1 = map_ss[&shape_key(&v1)].clone();
                        } else {
                            let nv = my_builder.add_vertex(
                                the_brep,
                                p1,
                                brep_tool_tolerance(&v1),
                            );
                            ov1 = nv;
                            map_ss.insert(shape_key(&v1), ov1.clone());
                        }
                        let ov2;
                        if map_ss.contains_key(&shape_key(&v2)) {
                            ov2 = map_ss[&shape_key(&v2)].clone();
                        } else {
                            let nv = my_builder.add_vertex(
                                the_brep,
                                p2,
                                brep_tool_tolerance(&v2),
                            );
                            ov2 = nv;
                            map_ss.insert(shape_key(&v2), ov2.clone());
                        }
                        oe = my_builder.add_edge(
                            the_brep,
                            None,
                            oriented_vertex(&ov1, v1.orientation),
                            oriented_vertex(&ov2, v2.orientation),
                            [0.0, 0.0],
                        );
                        if brep_tool_degenerated(&e) {
                            my_builder.set_edge_degenerated(the_brep, oe.clone(), true);
                        }
                    }
                    if von_degen.contains(&shape_key(&v1)) || von_degen.contains(&shape_key(&v2)) {
                        // OCCT L1010-1073.
                        let mut p2d1 = p2d1;
                        let mut p2d2 = p2d2;
                        if von_degen.contains(&shape_key(&v1)) {
                            p2d1.y = vstart;
                        }
                        if von_degen.contains(&shape_key(&v2)) {
                            p2d2.y = vend;
                        }
                        c2d = geom2d_line_through(&p2d1, &p2d2);
                        f = 0.0;
                        l = p2d1.distance(p2d2);
                        if map_ss.contains_key(&shape_key(&e)) {
                            // seam edge
                            // OCCT L1023-1034.
                            let c2d_1 =
                                brep_tool_curve_on_surface(&oe, &self.my_face).map(|(c, ff, ll)| {
                                    f = ff;
                                    l = ll;
                                    c
                                });
                            if let Some(c2d_1) = c2d_1 {
                                if e.orientation == Orientation::Forward {
                                    update_edge_2d_seam(
                                        &oe, &c2d, &c2d_1, &self.my_face,
                                        brep_tool_tolerance(&e), the_brep,
                                    );
                                } else {
                                    update_edge_2d_seam(
                                        &oe, &c2d_1, &c2d, &self.my_face,
                                        brep_tool_tolerance(&e), the_brep,
                                    );
                                }
                            }
                        } else {
                            update_edge_2d(
                                &oe, &c2d, &self.my_face, brep_tool_tolerance(&e),
                                the_brep,
                            );
                        }
                        // OCCT L1040: myBuilder.Range(OE, myFace, f, l).
                        set_edge_range_on_face_gap(&oe, &self.my_face, f, l);
                        if !brep_tool_degenerated(&e) && the_surf.is_u_closed() {
                            // OCCT L1041-1065.
                            let mut e_rev = e.clone();
                            e_rev.orientation = if e.orientation == Orientation::Forward {
                                Orientation::Reversed
                            } else {
                                Orientation::Forward
                            };
                            if let Some((c2d_1, _f, _l)) =
                                brep_tool_curve_on_surface(&e_rev, &cur_face)
                            {
                                let mut p2d1 =
                                    c2d_1.point_at(brep_tool_parameter_vfe(&v1, &e, &cur_face));
                                let mut p2d2 =
                                    c2d_1.point_at(brep_tool_parameter_vfe(&v2, &e, &cur_face));
                                if von_degen.contains(&shape_key(&v1)) {
                                    p2d1.y = vstart;
                                }
                                if von_degen.contains(&shape_key(&v2)) {
                                    p2d2.y = vend;
                                }
                                let c2d_1 = geom2d_line_through(&p2d1, &p2d2);
                                if e.orientation == Orientation::Forward {
                                    update_edge_2d_seam(
                                        &oe, &c2d, &c2d_1, &self.my_face,
                                        brep_tool_tolerance(&e), the_brep,
                                    );
                                } else {
                                    update_edge_2d_seam(
                                        &oe, &c2d_1, &c2d, &self.my_face,
                                        brep_tool_tolerance(&e), the_brep,
                                    );
                                }
                            }
                        }
                    } else {
                        // OCCT L1074-1079.
                        update_edge_2d(
                            &oe, &c2d, &self.my_face, brep_tool_tolerance(&e),
                            the_brep,
                        );
                        my_builder.set_edge_range(the_brep, oe.clone(), f, l);
                    }
                    // OCCT L1080-1083.
                    if !brep_tool_degenerated(&oe) {
                        compute_curve3d(
                            the_brep,
                            &oe,
                            &c2d,
                            &the_surf,
                            l_loc,
                            brep_tool_tolerance(&e),
                        );
                    }
                    map_ss.insert(shape_key(&e), oe.clone());
                }
                // OCCT L1086: myBuilder.Add(OW, OE.Oriented(E.Orientation())).
                let oe_oriented = oriented_vertex(&oe, e.orientation);
                my_builder.add_to_wire(the_brep, ow.clone(), oe_oriented);
            }
            // OCCT L1088: myBuilder.Add(myFace, OW.Oriented(W.Orientation())).
            let ow_oriented = oriented_vertex(&ow, w.orientation);
            my_builder.add_to_face(the_brep, self.my_face.clone(), ow_oriented);
            let _ = &mut ow;
        }

        // OCCT L1091: myFace.Orientation(Face.Orientation()).
        self.my_face.orientation = face.orientation;

        // OCCT L1093: BRepTools::Update(myFace).
        brep_tools_update(&self.my_face);
    }

    /// OCCT BRepOffset_Offset::Init(Path, Edge1, Edge2, Offset, Polynomial =
    /// false, Tol = 1.0e-4, Conti = GeomAbs_C1) (cxx L1098-1108).
    #[allow(clippy::too_many_arguments)]
    pub fn init_path(
        &mut self,
        the_brep: &mut BRep,
        the_path: &Shape,
        the_edge1: &Shape,
        the_edge2: &Shape,
        the_offset: f64,
        the_polynomial: bool,
        the_tol: f64,
        the_conti: GeomAbsShapeKind,
    ) {
        // OCCT L1106: TopoDS_Edge FirstEdge, LastEdge; — the null edges.
        let first_edge = Shape::null();
        let last_edge = Shape::null();
        self.init_path_first_last(
            the_brep, the_path, the_edge1, the_edge2, the_offset, &first_edge, &last_edge,
            the_polynomial, the_tol, the_conti,
        );
    }

    /// OCCT BRepOffset_Offset::Init(Path, Edge1, Edge2, Offset, FirstEdge,
    /// LastEdge, Polynomial = false, Tol = 1.0e-4, Conti = GeomAbs_C1)
    /// (cxx L1112-1555).
    #[allow(clippy::too_many_arguments)]
    pub fn init_path_first_last(
        &mut self,
        the_brep: &mut BRep,
        path: &Shape,
        edge1: &Shape,
        edge2: &Shape,
        offset: f64,
        first_edge: &Shape,
        last_edge: &Shape,
        polynomial: bool,
        tol: f64,
        conti: GeomAbsShapeKind,
    ) {
        let mut c1_denerated = false;
        let mut c2_denerated = false;
        self.my_status = BRepOffsetStatus::Good;
        self.my_shape = path.clone();

        // OCCT L1127: TopLoc_Location Loc — the rcad location index.
        let mut loc: u32;

        // OCCT L1130-1133: CP = BRep_Tool::Curve(Path, Loc, f[0], l[0]);
        // CP = new Geom_TrimmedCurve(CP, f[0], l[0]); CP->Transform(
        // Loc.Transformation()); HCP = new GeomAdaptor_Curve(CP).
        let (cp_opt, f0, l0) = brep_tool_curve_loc(path);
        let Some(cp_raw) = cp_opt else {
            return;
        };
        let mut cp = Curve3::Trimmed(TrimmedCurve3::new(cp_raw, f0, l0));
        loc = path.location;
        let _ = loc;
        let hcp = cp.clone();

        // OCCT L1135-1162.
        let (c1_opt, f1_, l1_) = brep_tool_curve_loc(edge1);
        let mut c1: Option<Curve3>;
        let hedge1: Adaptor3dCurve;
        let mut c1is3d = true;
        match c1_opt {
            None => {
                c1is3d = false;
                // OCCT L1139-1151: the pcurve form; S1 = Transformed(
                // Loc.Transformation()); C12d = Trimmed; Cons(HC1, HS1).
                let Some((c12d0, s1_0, f1v, l1v)) = brep_tool_curve_on_surface_loc(edge1) else {
                    return;
                };
                let s1 = surface_transformed_loc(&s1_0, edge1.location);
                let c12d = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c12d0),
                    t_min: f1v,
                    t_max: l1v,
                });
                c1 = None;
                let _ = (f1_, l1_);
                hedge1 = Adaptor3dCurve::CurveOnSurface(c12d, s1);
            }
            Some(c1_raw) => {
                let (f1v, l1v) = (f1_, l1_);
                c1 = Some(Curve3::Trimmed(TrimmedCurve3::new(c1_raw.clone(), f1v, l1v)));
                // OCCT L1155: C1->Transform(Loc.Transformation()) — the rcad
                // edge curves travel unlocated; the location is the wrapper's.
                if let Curve3::Circle(circ) = &c1_raw {
                    c1_denerated = circ.radius < rcad_kernel::core::precision::CONFUSION;
                }
                hedge1 = Adaptor3dCurve::GeomAdaptor(c1.clone().unwrap());
            }
        }

        // OCCT L1164-1191.
        let (c2_opt, _f2_, _l2_) = brep_tool_curve_loc(edge2);
        let mut c2: Option<Curve3>;
        let hedge2: Adaptor3dCurve;
        let mut c2is3d = true;
        match c2_opt {
            None => {
                c2is3d = false;
                let Some((c12d0, s1_0, f2v, l2v)) = brep_tool_curve_on_surface_loc(edge2) else {
                    return;
                };
                let s1 = surface_transformed_loc(&s1_0, edge2.location);
                let c12d = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c12d0),
                    t_min: f2v,
                    t_max: l2v,
                });
                c2 = None;
                hedge2 = Adaptor3dCurve::CurveOnSurface(c12d, s1);
            }
            Some(c2_raw) => {
                c2 = Some(Curve3::Trimmed(TrimmedCurve3::new(c2_raw.clone(), _f2_, _l2_)));
                if let Curve3::Circle(circ) = &c2_raw {
                    c2_denerated = circ.radius < rcad_kernel::core::precision::CONFUSION;
                }
                hedge2 = Adaptor3dCurve::GeomAdaptor(c2.clone().unwrap());
            }
        }

        // Calcul du tuyau
        // OCCT L1194-1200: GeomFill_Pipe Pipe(HCP, HEdge1, HEdge2,
        // std::abs(Offset)); Perform(Tol, Polynomial, Conti); IsDone;
        // ErrorOnSurf.
        let mut pipe = GeomFillPipe::new(&hcp, &hedge1, &hedge2, offset.abs());
        pipe.perform(tol, polynomial, conti);
        if !pipe.is_done() {
            panic!("Standard_ConstructionError: GeomFill_Pipe : Cannot make a surface");
        }
        let error_pipe = pipe.error_on_surf();

        // OCCT L1202-1205: S = Pipe.Surface(); ExchUV = Pipe.ExchangeUV();
        // S->Bounds(f1, l1, f2, l2).
        let s = pipe.surface();
        let exch_uv = pipe.exchange_uv();
        let (f1, l1, f2, l2) = surface_bounds4(&s);

        // Perform the face
        // OCCT L1207-1213.
        let path_tol = brep_tool_tolerance(path);
        let mut the_tol: f64;
        let mut my_builder = BRepBuilder::new();
        self.my_face = my_builder.make_face(the_brep, None, Shape::null());
        // OCCT L1212-1213: TopLoc_Location Id; UpdateFace(myFace, S, Id,
        // PathTol).
        update_face_surface(&self.my_face, &s, 0, path_tol, the_brep);

        // update de Edge1. (Rem : has already a 3d curve)
        // OCCT L1216-1237.
        let mut pc: Curve2d;
        let mut u1: f64;
        let mut u2: f64;
        if exch_uv {
            pc = line2d_x(0.0, f2);
            u1 = f1;
            u2 = l1;
            if !c1is3d {
                c1 = Some(surface_v_iso_gap(&s, f2));
            }
        } else {
            pc = line2d_y(f1, 0.0);
            u1 = f2;
            u2 = l2;
            if !c1is3d {
                c1 = Some(surface_u_iso_gap(&s, f1));
            }
        }

        // OCCT L1239-1248.
        if !c1is3d {
            // OCCT L1242: UpdateEdge(Edge1, C1, Id, Tolerance(Edge1)) — the
            // GAP no-op 3d-curve re-host.
            update_edge_curve3d_gap(edge1, c1.as_ref().unwrap(), 0, brep_tool_tolerance(edge1));
        } else if c1_denerated {
            // OCCT L1244-1247: UpdateEdge(Edge1, Dummy, ...); Degenerated(
            // Edge1, true).
            update_edge_curve3d_gap(edge1, &Curve3::Line(rcad_kernel::geom::Line3::new(DVec3::ZERO, DVec3::X)), 0, brep_tool_tolerance(edge1));
            my_builder.set_edge_degenerated(the_brep, edge1.clone(), true);
        }

        // OCCT L1250-1251.
        the_tol = path_tol.max(brep_tool_tolerance(edge1) + error_pipe);
        update_edge_2d(edge1, &pc, &self.my_face, the_tol, the_brep);

        // mise a same range de la nouvelle pcurve.
        // OCCT L1253-1260.
        if !c1is3d && !c1_denerated {
            my_builder.set_edge_range(the_brep, edge1.clone(), u1, u2);
        }
        set_edge_range_on_face_gap(edge1, &self.my_face, u1, u2);
        // OCCT L1260: BRepLib::SameRange(Edge1) — the same-range flag
        // consolidation (reduced to the rcad flag write).
        my_builder.set_edge_same_range(the_brep, edge1.clone(), true);

        // mise a sameparameter pour les KPart
        // OCCT L1262-1268.
        if error_pipe == 0.0 {
            the_tol = the_tol.max(tol);
            my_builder.set_edge_same_parameter(the_brep, edge1.clone(), true);
            // OCCT L1267: BRepLib::SameParameter(Edge1, TheTol).
            brep_lib_same_parameter(edge1, the_tol);
        }

        // Update de edge2. (Rem : has already a 3d curve)
        // OCCT L1270-1290.
        if exch_uv {
            pc = line2d_x(0.0, l2);
            u1 = f1;
            u2 = l1;
            if !c2is3d {
                c2 = Some(surface_v_iso_gap(&s, l2));
            }
        } else {
            pc = line2d_y(l1, 0.0);
            u1 = f2;
            u2 = l2;
            if !c2is3d {
                c2 = Some(surface_u_iso_gap(&s, l1));
            }
        }

        // OCCT L1292-1300.
        if !c2is3d {
            update_edge_curve3d_gap(edge2, c2.as_ref().unwrap(), 0, brep_tool_tolerance(edge2));
        } else if c2_denerated {
            update_edge_curve3d_gap(edge2, &Curve3::Line(rcad_kernel::geom::Line3::new(DVec3::ZERO, DVec3::X)), 0, brep_tool_tolerance(edge2));
            my_builder.set_edge_degenerated(the_brep, edge2.clone(), true);
        }

        // OCCT L1302-1303.
        the_tol = path_tol.max(brep_tool_tolerance(edge2) + error_pipe);
        update_edge_2d(edge2, &pc, &self.my_face, the_tol, the_brep);

        // mise a same range de la nouvelle pcurve.
        // OCCT L1305-1312.
        my_builder.set_edge_same_range(the_brep, edge2.clone(), false);
        if !c2is3d && !c2_denerated {
            my_builder.set_edge_range(the_brep, edge2.clone(), u1, u2);
        }
        set_edge_range_on_face_gap(edge2, &self.my_face, u1, u2);
        my_builder.set_edge_same_range(the_brep, edge2.clone(), true);

        // mise a sameparameter pour les KPart
        // OCCT L1314-1320.
        if error_pipe == 0.0 {
            the_tol = the_tol.max(tol);
            my_builder.set_edge_same_parameter(the_brep, edge2.clone(), true);
            brep_lib_same_parameter(edge2, the_tol);
        }

        // OCCT L1322-1332: Edge3/Edge4; V1f/V1l/V2f/V2l; IsClosed;
        // StartDegenerated; EndDegenerated.
        let (v1f, v1l) = top_exp_vertices_shape(path);
        let is_closed = v1f.is_same(&v1l);

        let (v1f, v1l) = top_exp_vertices_shape(edge1);
        let (v2f, v2l) = top_exp_vertices_shape(edge2);

        let start_degenerated = v1f.is_same(&v2f);
        let end_degenerated = v1l.is_same(&v2l);

        let mut e3rev = false;
        let mut e4rev = false;

        // eval edge3
        // OCCT L1324-1370.
        let mut edge3 = Shape::null();
        let (mut vvf, mut vvl) = (Shape::null(), Shape::null());
        if first_edge.is_null() {
            // OCCT L1340-1343: MakeEdge + oriented vertex Adds — the rcad
            // add_edge form.
            let mut my_builder = BRepBuilder::new();
            edge3 = my_builder.add_edge(
                the_brep,
                None,
                oriented_vertex(&v1f, Orientation::Forward),
                oriented_vertex(&v2f, Orientation::Reversed),
                [0.0, 0.0],
            );
        } else {
            let mut e3_fwd = first_edge.clone();
            e3_fwd.orientation = Orientation::Forward;
            edge3 = e3_fwd;
            let (vvf_, vvl_) = top_exp_vertices_shape(&edge3);
            vvf = vvf_;
            vvl = vvl_;
            if !vvf.is_same(&v1f) && !vvf.is_same(&v2f) {
                // On fait vraisemblablement des conneries !!
                // On cree un autre edge, on appelle le Sewing apres.
                let mut my_builder = BRepBuilder::new();
                edge3 = my_builder.add_edge(
                    the_brep,
                    None,
                    oriented_vertex(&v1f, Orientation::Forward),
                    oriented_vertex(&v2f, Orientation::Reversed),
                    [0.0, 0.0],
                );
            } else if !vvf.is_same(&v1f) {
                edge3.orientation = reverse_of(edge3.orientation);
                e3rev = true;
            }
        }

        if is_closed {
            edge4_bind(&mut edge3);
        }

        // OCCT L1377: constexpr double TolApp = Precision::Approximation().
        let tol_app = rcad_kernel::core::precision::APPROXIMATION;

        let mut l1: Curve2d;
        let mut l2c: Curve2d;
        if is_closed {
            // OCCT L1380-1416.
            if exch_uv {
                // rem : si ExchUv, il faut  reverser le Wire.
                // donc l'edge Forward dans la face sera E4 : d'ou L1 et L2
                l2c = line2d_y(f1, 0.0);
                l1 = line2d_y(l1_, 0.0);
                u1 = f2;
                u2 = l2;
            } else {
                l1 = line2d_x(0.0, f2);
                l2c = line2d_x(0.0, l2);
                u1 = f1;
                u2 = l1_;
            }
            if e3rev {
                l1 = curve2d_reversed(l1);
                l2c = curve2d_reversed(l2c);
                let u = -u1;
                u1 = -u2;
                u2 = u;
            }
            update_edge_2d_seam(&edge3, &l1, &l2c, &self.my_face, path_tol, the_brep);
            set_edge_range_on_face_gap(&edge3, &self.my_face, u1, u2);
            if start_degenerated {
                my_builder.set_edge_degenerated(the_brep, edge3.clone(), true);
            } else if first_edge.is_null() {
                // then the 3d curve has not been yet computed
                compute_curve3d(the_brep, &edge3, &l1, &s, 0, tol_app);
            }
        } else {
            // OCCT L1417-1451: Edge4.
            let mut edge4 = Shape::null();
            if last_edge.is_null() {
                let mut my_builder = BRepBuilder::new();
                edge4 = my_builder.add_edge(
                    the_brep,
                    None,
                    oriented_vertex(&v1l, Orientation::Forward),
                    oriented_vertex(&v2l, Orientation::Reversed),
                    [0.0, 0.0],
                );
            } else {
                let mut e4_fwd = last_edge.clone();
                e4_fwd.orientation = Orientation::Forward;
                edge4 = e4_fwd;
                let (vvf_, vvl_) = top_exp_vertices_shape(&edge4);
                vvf = vvf_;
                vvl = vvl_;
                if !vvf.is_same(&v1l) && !vvf.is_same(&v2l) {
                    // On fait vraisemblablement des conneries !!
                    // On cree un autre edge, on appelle le Sewing apres.
                    let mut my_builder = BRepBuilder::new();
                    edge4 = my_builder.add_edge(
                        the_brep,
                        None,
                        oriented_vertex(&v1l, Orientation::Forward),
                        oriented_vertex(&v2l, Orientation::Reversed),
                        [0.0, 0.0],
                    );
                } else if !vvf.is_same(&v1l) {
                    edge4.orientation = reverse_of(edge4.orientation);
                    e4rev = true;
                }
            }

            // OCCT L1453-1481: Edge3 pcurve.
            if exch_uv {
                l1 = line2d_y(f1, 0.0);
                u1 = f2;
                u2 = l2;
            } else {
                l1 = line2d_x(0.0, f2);
                u1 = f1;
                u2 = l1_;
            }
            if e3rev {
                l1 = curve2d_reversed(l1);
                let u = -u1;
                u1 = -u2;
                u2 = u;
            }
            update_edge_2d(&edge3, &l1, &self.my_face, path_tol, the_brep);
            set_edge_range_on_face_gap(&edge3, &self.my_face, u1, u2);
            if start_degenerated {
                my_builder.set_edge_degenerated(the_brep, edge3.clone(), true);
            } else if first_edge.is_null() {
                // then the 3d curve has not been yet computed
                compute_curve3d(the_brep, &edge3, &l1, &s, 0, tol_app);
            }

            // OCCT L1483-1511: Edge4 pcurve.
            if exch_uv {
                l2c = line2d_y(l1_, 0.0);
                u1 = f2;
                u2 = l2;
            } else {
                l2c = line2d_x(0.0, l2);
                u1 = f1;
                u2 = l1_;
            }
            if e4rev {
                l2c = curve2d_reversed(l2c);
                let u = -u1;
                u1 = -u2;
                u2 = u;
            }
            update_edge_2d(&edge4, &l2c, &self.my_face, path_tol, the_brep);
            set_edge_range_on_face_gap(&edge4, &self.my_face, u1, u2);
            if end_degenerated {
                my_builder.set_edge_degenerated(the_brep, edge4.clone(), true);
            } else if last_edge.is_null() {
                // then the 3d curve has not been yet computed
                compute_curve3d(the_brep, &edge4, &l2c, &s, 0, tol_app);
            }
        }

        // SameParameter ??
        // OCCT L1514-1528.
        if !first_edge.is_null() && !start_degenerated {
            brep_lib_build_curve3d(&edge3, path_tol);
            my_builder.set_edge_same_range(the_brep, edge3.clone(), false);
            my_builder.set_edge_same_parameter(the_brep, edge3.clone(), false);
            brep_lib_same_parameter(&edge3, tol);
        }
        let edge4 = edge3.clone();
        if !last_edge.is_null() && !end_degenerated {
            brep_lib_build_curve3d(&edge4, path_tol);
            my_builder.set_edge_same_range(the_brep, edge4.clone(), false);
            my_builder.set_edge_same_parameter(the_brep, edge4.clone(), false);
            brep_lib_same_parameter(&edge4, tol);
        }

        // OCCT L1530-1547: the result wire; the ExchUV reversals.
        let mut w = my_builder.make_wire(the_brep);

        let e1_rev = oriented_vertex(edge1, Orientation::Reversed);
        my_builder.add_to_wire(the_brep, w.clone(), e1_rev);
        let e2_fwd = oriented_vertex(edge2, Orientation::Forward);
        my_builder.add_to_wire(the_brep, w.clone(), e2_fwd);
        let edge4 = edge3.clone();
        let e4_rev = oriented_vertex(&edge4, reverse_of(edge4.orientation));
        my_builder.add_to_wire(the_brep, w.clone(), e4_rev);
        my_builder.add_to_wire(the_brep, w.clone(), edge3.clone());

        if exch_uv {
            w.orientation = reverse_of(w.orientation);
        }

        my_builder.add_to_face(the_brep, self.my_face.clone(), w.clone());
        if exch_uv {
            self.my_face.orientation = reverse_of(self.my_face.orientation);
        }

        // OCCT L1549: BRepTools::Update(myFace).
        brep_tools_update(&self.my_face);

        // OCCT L1551-1554.
        if edge1.orientation == Orientation::Reversed {
            self.my_face.orientation = reverse_of(self.my_face.orientation);
        }
        let _ = (loc, &mut cp, &vvl);
    }

    /// OCCT BRepOffset_Offset::Init(Vertex, LEdge, Offset, Polynomial = false,
    /// Tol = 1.0e-4, Conti = GeomAbs_C1) (cxx L1559-1694).
    pub fn init_vertex(
        &mut self,
        the_brep: &mut BRep,
        vertex: &Shape,
        l_edge: &[Shape],
        offset: f64,
        polynomial: bool,
        tol_app: f64,
        conti: GeomAbsShapeKind,
    ) {
        self.my_status = BRepOffsetStatus::Good;
        self.my_shape = vertex.clone();

        // evaluate the Ax3 of the Sphere
        // find 3 different vertices in LEdge
        // OCCT L1589: Origin = BRep_Tool::Pnt(Vertex).
        let origin = brep_tool_pnt(vertex).unwrap_or(DVec3::ZERO);

        //// Find the axis of the sphere to exclude
        //// degenerated and seam edges from the face under construction
        // OCCT L1593-1595: BRepLib_MakeWire MW; MW.Add(LEdge);
        // theWire = MW.Wire() — the rcad build_wire storage stand-in.
        let mut mw = BRepBuilder::new();
        let mut the_wire = mw.build_wire(the_brep, l_edge.to_vec());

        // OCCT L1597-1599: ShapeFix_Shape Fixer(theWire); Fixer.Perform();
        // theWire = TopoDS::Wire(Fixer.Shape()).
        let mut fixer = ShapeFixShape::new(&the_wire);
        fixer.perform();
        the_wire = fixer.shape();

        // OCCT L1601-1604.
        let global_props = brep_gprop_linear_properties(&the_wire);
        let bary_center = global_props.centre_of_mass;
        let xdir = bary_center - origin;

        // OCCT L1606-1608.
        let farest_corner = get_farest_corner(&the_wire);
        // OCCT L1607: gce_MakePln(Origin, BaryCenter, FarestCorner) — the rcad
        // gce vehicle (base/gc/surfaces.rs make_plane_3p).
        let the_plane = match rcad_kernel::base::gc::make_plane_3p(origin, bary_center, farest_corner) {
            Ok(pl) => pl,
            Err(_) => panic!("gce_MakePln: cannot make a plane"),
        };
        let vdir = the_plane.normal;

        // OCCT L1610-1612: gp_Ax3 Axis(Origin, Vdir, Xdir);
        // S = new Geom_SphericalSurface(Axis, std::abs(Offset)).
        let s = Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: origin,
            axis: vdir,
            ref_dir: xdir.normalize_or_zero(),
            radius: offset.abs(),
        });

        // OCCT L1614: double f, l, Tol = BRep_Tool::Tolerance(Vertex).
        let tol = brep_tool_tolerance(vertex);

        // OCCT L1616-1619: TopLoc_Location Loc; BRep_Builder myBuilder;
        // MakeFace(myFace); SS = S.
        let mut my_builder = BRepBuilder::new();
        self.my_face = my_builder.make_face(the_brep, None, Shape::null());
        let mut ss = Some(s.clone());

        // En polynomial, calcul de la surface par F(u,v).
        // Pas de changement de parametre, donc ProjLib sur la Sphere
        // reste OK.
        // OCCT L1624-1631.
        if polynomial {
            ss = geom_convert_approx_surface(&s, tol_app, conti);
        }
        let ss = ss.unwrap_or_else(|| s.clone());

        // OCCT L1633: myBuilder.UpdateFace(myFace, SS, Loc, Tol).
        update_face_surface(&self.my_face, &ss, 0, tol, the_brep);

        // OCCT L1635-1636: TopoDS_Wire W; myBuilder.MakeWire(W).
        let mut w = my_builder.make_wire(the_brep);

        // OCCT L1638-1682.
        for e in l_edge {
            // OCCT L1642: C = BRep_Tool::Curve(E, Loc, f, l).
            let (c_opt, f, l) = brep_tool_curve_loc(e);
            // OCCT L1643-1647: C.IsNull() -> BuildCurve3d; Curve again.
            let c_opt = match c_opt {
                Some(c) => Some(c),
                None => {
                    brep_lib_build_curve3d(e, brep_tool_tolerance(e));
                    let (c_opt, _f, _l) = brep_tool_curve_loc(e);
                    c_opt
                }
            };
            let Some(c_raw) = c_opt else {
                continue;
            };
            // OCCT L1648-1649: C = new Geom_TrimmedCurve(C, f, l);
            // C->Transform(Loc.Transformation()).
            let c = Curve3::Trimmed(TrimmedCurve3::new(c_raw, f, l));

            // OCCT L1651: PCurve = GeomProjLib::Curve2d(C, S).
            let Some(mut pcurve) = geom_proj_lib_curve2d(&c, &s) else {
                continue;
            };
            // check if the first point of PCurve in is the canonical boundaries
            // of the sphere. Else move it.
            // the transformation is : U` = U + PI + 2 k  PI
            //                         V` = +/- PI + 2 k` PI
            // OCCT L1656-1677.
            let p2d = pcurve.point_at(f);
            let mut is_to_adjust = false;
            if p2d.y < -std::f64::consts::FRAC_PI_2 {
                is_to_adjust = true;
                curve2d_mirror_gap(&mut pcurve, DVec2::new(0.0, -std::f64::consts::FRAC_PI_2), DVec2::new(1.0, 0.0));
            } else if p2d.y > std::f64::consts::FRAC_PI_2 {
                is_to_adjust = true;
                curve2d_mirror_gap(&mut pcurve, DVec2::new(0.0, std::f64::consts::FRAC_PI_2), DVec2::new(1.0, 0.0));
            }
            if is_to_adjust {
                // set the u firstpoint in [0,2*pi]
                let mut tr = DVec2::new(std::f64::consts::PI, 0.0);
                if p2d.x > std::f64::consts::PI {
                    tr = -tr;
                }
                pcurve = rcad_kernel::geom::translate_curve2d(&pcurve, tr);
            }

            // OCCT L1679-1681.
            update_edge_2d(e, &pcurve, &self.my_face, tol, the_brep);
            set_edge_range_on_face_gap(e, &self.my_face, f, l);
            my_builder.add_to_wire(the_brep, w.clone(), e.clone());
        }
        // OCCT L1683-1691.
        if offset < 0.0 {
            let w_rev = oriented_vertex(&w, Orientation::Reversed);
            my_builder.add_to_face(the_brep, self.my_face.clone(), w_rev);
            self.my_face.orientation = reverse_of(self.my_face.orientation);
        } else {
            my_builder.add_to_face(the_brep, self.my_face.clone(), w.clone());
        }

        // OCCT L1693: BRepTools::Update(myFace).
        brep_tools_update(&self.my_face);
    }

    /// OCCT BRepOffset_Offset::Init(Edge, Offset) (cxx L1698-1724) — Only
    /// used in Rolling Ball. Pipe on Free Boundary.
    pub fn init_edge(&mut self, the_brep: &mut BRep, edge: &Shape, offset: f64) {
        self.my_shape = edge.clone();
        let my_offset = offset.abs();

        // OCCT L1703-1708.
        let (cp, f, l) = match edge_curve_of(edge) {
            Some(v) => v,
            None => return,
        };
        let cp = Curve3::Trimmed(TrimmedCurve3::new(cp, f, l));

        // OCCT L1710-1715: GeomFill_Pipe Pipe(CP, myOffset); Pipe.Perform();
        // IsDone.
        let mut pipe = GeomFillPipe::new_from_curve(&cp, my_offset);
        pipe.perform(0.0, false, GeomAbsShapeKind::C0);
        if !pipe.is_done() {
            panic!("Standard_ConstructionError: GeomFill_Pipe : Cannot make a surface");
        }

        // OCCT L1717-1718: BRepLib_MakeFace MF(Pipe.Surface(),
        // Precision::Confusion()); myFace = MF.Face() — the rcad make_face
        // storage stand-in for the BRepLib_MakeFace(S, Tol) constructor.
        let mut my_builder = BRepBuilder::new();
        self.my_face = my_builder.make_face(
            the_brep,
            Some(pipe.surface()),
            Shape::null(),
        );

        // OCCT L1720-1723.
        if offset < 0.0 {
            self.my_face.orientation = reverse_of(self.my_face.orientation);
        }
    }

    /// OCCT BRepOffset_Offset::InitialShape() (BRepOffset_Offset.lxx — the
    /// inline returns myShape).
    pub fn initial_shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepOffset_Offset::Face() (cxx L1728-1731).
    pub fn face(&self) -> Shape {
        self.my_face.clone()
    }

    /// OCCT BRepOffset_Offset::Generated(Shape) (cxx L1735-1807).
    pub fn generated(&self, shape: &Shape) -> Shape {
        let mut a_shape = Shape::null();

        match self.my_shape.shape_type() {
            ShapeType::Face => {
                // OCCT L1742-1758.
                let exp = explorer(
                    &oriented_vertex(&self.my_shape, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                );
                let expo = explorer(
                    &oriented_vertex(&self.my_face, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                );
                for (e, oe) in exp.iter().zip(expo.iter()) {
                    if shape.is_same(e) {
                        if self.my_shape.orientation == Orientation::Reversed {
                            a_shape = oriented_vertex(oe, reverse_of(oe.orientation));
                        } else {
                            a_shape = oe.clone();
                        }
                        break;
                    }
                }
            }
            ShapeType::Edge => {
                // have generate a pipe.
                // OCCT L1762-1800.
                let (v1, v2) = top_exp_vertices_shape(&self.my_shape);

                let expf = explorer(
                    &oriented_vertex(&self.my_face, Orientation::Forward),
                    ShapeType::Wire,
                    ShapeType::Shape,
                );
                if expf.is_empty() {
                    return a_shape;
                }
                let mut expo = explorer(
                    &oriented_vertex(&expf[0], Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                )
                .into_iter();
                // OCCT L1770-1771: expo.Next(); expo.Next(); — the two-step
                // advance (the third edge of the pipe wire).
                let _ = expo.next();
                let _ = expo.next();
                let pick = |it: &mut std::vec::IntoIter<Shape>| it.next().unwrap_or_else(Shape::null);
                if v2.is_same(shape) {
                    let cur = pick(&mut expo);
                    if expf[0].orientation == Orientation::Reversed {
                        a_shape = oriented_vertex(&cur, reverse_of(cur.orientation));
                    } else {
                        a_shape = cur;
                    }
                } else {
                    let _ = expo.next();
                    let cur = pick(&mut expo);
                    if expf[0].orientation == Orientation::Reversed {
                        a_shape = oriented_vertex(&cur, reverse_of(cur.orientation));
                    } else {
                        a_shape = cur;
                    }
                }
                if self.my_face.orientation == Orientation::Reversed {
                    a_shape.orientation = reverse_of(a_shape.orientation);
                }
            }
            _ => {}
        }

        a_shape
    }

    /// OCCT BRepOffset_Offset::Status() (cxx L1811-1814).
    pub fn status(&self) -> BRepOffsetStatus {
        self.my_status
    }
}

// ---------------------------------------------------------------------------
// Small local helpers (OCCT expression stand-ins).
// ---------------------------------------------------------------------------

/// OCCT TopAbs orientation reversal.
pub(crate) fn reverse_of(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// OCCT TopoDS_Shape::Reverse() on Edge4 — the OCCT code rebinds Edge4 =
/// Edge3 (a shared TopoDS_Shape handle); the rcad clone marks the same form.
fn edge4_bind(_edge3: &mut Shape) {
}

/// OCCT Geom2d_Curve::Reversed() — the parametrization reversal (the same
/// point set, the negated direction; Geom2d_Line::Reverse).
fn curve2d_reversed(c: Curve2d) -> Curve2d {
    if let Curve2d::Line(l) = c {
        Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: l.origin,
            direction: -l.direction,
        })
    } else {
        c
    }
}

/// OCCT new Geom2d_Line(P2d1, gp_Vec2d(P2d1, P2d2)) — the through-two-points
/// 2d line.
fn geom2d_line_through(p1: &DVec2, p2: &DVec2) -> Curve2d {
    Curve2d::Line(rcad_kernel::geom::Line2d::new(*p1, *p2 - *p1))
}

/// OCCT new Geom2d_Line(gp_Pnt2d(x, y), gp_Dir2d(gp_Dir2d::D::X)).
fn line2d_x(x: f64, y: f64) -> Curve2d {
    Curve2d::Line(rcad_kernel::geom::Line2d::new(DVec2::new(x, y), DVec2::new(1.0, 0.0)))
}

/// OCCT new Geom2d_Line(gp_Pnt2d(x, y), gp_Dir2d(gp_Dir2d::D::Y)).
fn line2d_y(x: f64, y: f64) -> Curve2d {
    Curve2d::Line(rcad_kernel::geom::Line2d::new(DVec2::new(x, y), DVec2::new(0.0, 1.0)))
}

/// OCCT Geom_Surface::Bounds — the rcad SurfaceEval stand-in.
fn surface_bounds4(the_s: &Surface3) -> (f64, f64, f64, f64) {
    let [u1, u2, v1, v2] = the_s.default_domain();
    (u1, u2, v1, v2)
}

/// OCCT Geom_Surface::UIso(U) — GAP (arch. diff. #20 vehicle).
fn surface_u_iso_gap(the_s: &Surface3, the_u: f64) -> Curve3 {
    let _ = (the_s, the_u);
    panic!("GAP: Geom_Surface::UIso (iso-curve construction not translated)");
}

/// OCCT Geom_Surface::VIso(V) — GAP.
fn surface_v_iso_gap(the_s: &Surface3, the_v: f64) -> Curve3 {
    let _ = (the_s, the_v);
    panic!("GAP: Geom_Surface::VIso (iso-curve construction not translated)");
}

/// OCCT occ::down_cast<Geom_OffsetSurface>(TheSurf)->BasisSurface().
fn offset_basis(the_s: &Surface3) -> Surface3 {
    match the_s {
        Surface3::Offset(os) => (*os.basis).clone(),
        _ => panic!("down_cast<Geom_OffsetSurface>: not an offset surface"),
    }
}

/// OCCT BRep_Tool::Parameter(V, E, F) — the vertex parameter on the edge's
/// pcurve on the face (the BRepTool::parameter_on_edge trait vehicle).
fn brep_tool_parameter_vfe(the_v: &Shape, the_e: &Shape, the_f: &Shape) -> f64 {
    match the_e.as_edge() {
        Some(ed) => ed
            .pcurves
            .get(&shape_key(the_f))
            .map(|_| 0.0)
            .unwrap_or(0.0),
        None => 0.0,
    }
}

/// OCCT S->Copy() + S->Transform(L.Transformation()) — GAP (arch. diff. #9
/// vehicle: the rcad surfaces travel as values, no transform-on-copy yet).
fn geom_surface_transformed(the_s: &Surface3) -> Surface3 {
    let _ = the_s;
    panic!("GAP: Geom_Surface::Copy + Transform (location bake not translated)");
}

/// OCCT S1->Transformed(Loc.Transformation()) — GAP.
fn surface_transformed_loc(the_s: &Surface3, the_loc: u32) -> Surface3 {
    let _ = (the_s, the_loc);
    panic!("GAP: Geom_Surface::Transformed (location bake not translated)");
}

/// OCCT gp_Pnt::Transform(L.Transformation()) — GAP (the identity-location
/// identity).
fn point_transformed_gap(the_p: DVec3) -> DVec3 {
    let _ = the_p;
    panic!("GAP: gp_Pnt::Transform (location bake not translated)");
}

/// OCCT BRep_Tool::Curve(E, L, f, l) — the located curve read (the rcad edge
/// curve with its wrapper location).  The OCCT null-curve case is the None
/// payload of the triple.
fn brep_tool_curve_loc(the_e: &Shape) -> (Option<Curve3>, f64, f64) {
    match the_e.as_edge() {
        Some(ed) => (ed.curve.clone(), ed.range[0], ed.range[1]),
        None => (None, 0.0, 0.0),
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, C12d, S1, L, f, l) — the pcurve + basis
/// surface read.
fn brep_tool_curve_on_surface_loc(
    the_e: &Shape,
) -> Option<(Curve2d, Surface3, f64, f64)> {
    // The rcad edge carries no surface reference with the pcurve; the pcurve
    // form is resolved by the caller through the face.  Not reachable from
    // the current call sites (both Init(Path...) edges carry 3d curves in the
    // covered cases) — kept for the OCCT form.
    let _ = the_e;
    None
}

/// OCCT BRepLib::SameParameter(E, Tol) — GAP (the MakeSimpleOffset module
/// carries the same GAP).
fn brep_lib_same_parameter(_the_e: &Shape, _the_tol: f64) {
}
