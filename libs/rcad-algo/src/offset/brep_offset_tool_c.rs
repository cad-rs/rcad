// OCCT BRepOffset_Tool.cxx L2580-L4659 — module c of the 1:1 translation
// (see the architecture-difference list in brep_offset_tool.rs; numbering
// continues there).
//
// Module c carries the statics MakeFace / EnlargeGeometry / UpdatePCurves /
// CompactUVBounds and the class methods CheckBounds / EnLargeFace /
// TryParameter(static) / MapVertexEdges / BuildNeighbour / ExtentFace /
// Deboucle3D / CorrectOrientation / CheckPlanesNormals plus the file
// statics PerformPlanes / UpdateVertexTolerances.
//
// Additional architecture notes for this module:
// 28. BRep_Builder::UpdateEdge(E, C2d1, C2d2, F, Tol) (the seam
//     two-pcurve form) — the rcad pcurve store carries one pcurve per face
//     key, so the second (reversed) pcurve overwrites the first (annotated
//     at the call sites; the store extension is staged).
// 29. BRepTopAdaptor_FClass2d -> the topalgo fclass2d over FaceShapeSource
//     (the loc_ope_wires_on_shape_b.rs #8 precedent);
//     GCPnts_QuasiUniformDeflection -> the reduced uniform sampling.
// 30. IntTools_FaceFace::Perform(theFace1, theFace2) (the bare-face form)
//     — the rcad IntTools FaceFace is DS-bound; the PerformPlanes carrier
//     takes the OCCT !IsDone() path behind the GAP annotation.
// 31. BRepLib::Update(S) (the Maj des UVPoints of ExtentFace) and
//     BRepTools::DetectClosedness — GAP no-ops / reduced re-hosts
//     (annotated at the call sites).

use std::collections::HashMap;

use glam::DVec2;
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, SurfaceEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::topo::topods::{Orientation, ShapeType, State, TShape};
use rcad_kernel::topo_shape::Shape;

use rcad_kernel::math::el::elslib_cone_parameters;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_range, brep_tool_tolerance,
};

use super::brep_offset_tool::*;
use super::brep_offset_tool_d::{BRepOffsetAnalyse, FClass2dCarrier};
use super::brep_offset_tool_b::{
    b_update_vertex_on_edge, extent_edge_tool, inter2d, inter3d, project_vertex_on_edge, try_project,
};

// ---------------------------------------------------------------------------
// Local builder re-hosts.
// ---------------------------------------------------------------------------

/// OCCT Precision::IsInfinite(theVal) — the rcad INFINITE_VALUE form.
fn precision_is_infinite(the_val: f64) -> bool {
    the_val >= rcad_kernel::precision::INFINITE_VALUE
        || the_val <= -rcad_kernel::precision::INFINITE_VALUE
}

/// OCCT BRep_Builder::UpdateEdge(E, Tol) — the edge tolerance SET form
/// (TE->Tolerance(TE->Tolerance() * 10.) writes the value directly).
fn b_set_edge_tolerance(e: &mut Shape, tol: f64) {
    if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut e.data) {
        ed.tolerance = tol;
    }
}

/// OCCT BRep_Builder::UpdateFace(F, S, L, Tol) — the face surface update
/// (the rcad face carries the surface directly; location-baked, annotated).
fn b_update_face_surface(f: &mut Shape, s: &Surface3, tol: f64) {
    if let TShape::Face(fd) = std::sync::Arc::make_mut(&mut f.data) {
        fd.surface = Some(s.clone());
        fd.tolerance = fd.tolerance.max(tol);
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, F, Tol) — remove the pcurve entry
/// (the NullPCurve form of UpdatePCurves).
fn b_remove_edge_pcurve(e: &mut Shape, f: &Shape) {
    let key = bat::shape_key(f);
    if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut e.data) {
        ed.pcurves.shift_remove(&key);
    }
}

/// OCCT BRep_Builder::MakeFace(F, S, Tol) — the bare surface face.
fn b_make_face_surface(s: &Surface3, tol: f64) -> Shape {
    Shape {
        data: std::sync::Arc::new(TShape::Face(
            rcad_kernel::topo::topods::TFaceData {
                my_shapes: Vec::new(),
                flags: rcad_kernel::topo::topods::tshape_flags::DEFAULT,
                surface: Some(s.clone()),
                surface_location: 0,
                outer_wire: Shape::null(),
                inner_wires: Vec::new(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: Vec::new(),
                tolerance: tol,
                natural_restriction: false,
            },
        )),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::MakeSolid(SS) — an empty solid TShape.
fn b_make_solid() -> Shape {
    Shape {
        data: std::sync::Arc::new(TShape::Solid(
            rcad_kernel::topo::topods::TSolidData {
                my_shapes: Vec::new(),
                flags: rcad_kernel::topo::topods::tshape_flags::DEFAULT,
                shells: Vec::new(),
                internal_vertices: Vec::new(),
                internal_edges: Vec::new(),
            },
        )),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::Add(Solid, Shell).
fn b_add_solid_shell(solid: &mut Shape, shell: &Shape) {
    if let TShape::Solid(sd) = std::sync::Arc::make_mut(&mut solid.data) {
        sd.my_shapes.push(shell.clone());
        sd.shells.push(shell.clone());
    }
}

// ---------------------------------------------------------------------------
// OCCT static MakeFace (cxx L2580-2919).
// ---------------------------------------------------------------------------

/// OCCT static MakeFace(S, Um, UM, Vm, VM, uclosed, vclosed, isVminDegen,
/// isVmaxDegen, F) (cxx L2580-2919) — the bounding face of the enlarged
/// surface.
#[allow(clippy::too_many_arguments)]
fn make_face(
    s: &Surface3,
    um: f64,
    um2: f64,
    vm: f64,
    vm2: f64,
    uclosed: bool,
    vclosed: bool,
    is_vmin_degen: bool,
    is_vmax_degen: bool,
    f: &mut Shape,
) {
    let u_min = um;
    let u_max = um2;
    let v_min = vm;
    let v_max = vm2;

    // compute infinite flags
    // OCCT L2597-2600: Precision::IsNegativeInfinite / IsPositiveInfinite.
    let neg_inf = -rcad_kernel::precision::INFINITE_VALUE;
    let pos_inf = rcad_kernel::precision::INFINITE_VALUE;
    let umininf = u_min <= neg_inf;
    let umaxinf = u_max >= pos_inf;
    let vmininf = v_min <= neg_inf;
    let vmaxinf = v_max >= pos_inf;

    // degenerated flags (for cones)
    // OCCT L2603-2624: the cone apex probe through ElSLib::Parameters.
    let mut vmindegen = is_vmin_degen;
    let mut vmaxdegen = is_vmax_degen;
    let mut the_surf: &Surface3 = s;
    if let Surface3::Trimmed(ts) = s {
        the_surf = ts.basis.as_ref();
    }
    if let Surface3::Cone(cone) = the_surf {
        let the_apex = cone.apex;
        let y_dir = cone.axis.cross(cone.ref_dir);
        let (_u_apex, v_apex) = elslib_cone_parameters(
            the_apex,
            cone.apex,
            cone.ref_dir,
            y_dir,
            cone.axis,
            cone.radius,
            cone.half_angle_rad,
        );
        if (v_min - v_apex).abs() <= rcad_kernel::precision::CONFUSION {
            vmindegen = true;
        }
        if (v_max - v_apex).abs() <= rcad_kernel::precision::CONFUSION {
            vmaxdegen = true;
        }
    }

    // compute vertices
    // OCCT L2627-2628: BRep_Builder B; tol = Precision::Confusion().
    let tol = rcad_kernel::precision::CONFUSION;

    let mut v00 = Shape::null();
    let mut v10 = Shape::null();
    let mut v11 = Shape::null();
    let mut v01 = Shape::null();

    // OCCT L2632-2653: the corner vertices (B.MakeVertex(V, S->Value, tol)).
    let make_vertex_at = |u: f64, v: f64| -> Shape {
        let mut vtx = bat::builder_make_vertex();
        let p = s.point_at(u, v);
        bat::builder_update_vertex_point_tol(&mut vtx, p, tol);
        vtx
    };
    if !umininf {
        if !vmininf {
            v00 = make_vertex_at(u_min, v_min);
        }
        if !vmaxinf {
            v01 = make_vertex_at(u_min, v_max);
        }
    }
    if !umaxinf {
        if !vmininf {
            v10 = make_vertex_at(u_max, v_min);
        }
        if !vmaxinf {
            v11 = make_vertex_at(u_max, v_max);
        }
    }

    // OCCT L2655-2674: the closed / degenerate vertex merges.
    if uclosed {
        v10 = v00.clone();
        v11 = v01.clone();
    }
    if vclosed {
        v01 = v00.clone();
        v11 = v10.clone();
    }
    if vmindegen {
        v10 = v00.clone();
    }
    if vmaxdegen {
        v11 = v01.clone();
    }

    // make the lines
    // OCCT L2677-2693: the Geom2d_Line iso forms.
    let y_line = |x: f64| Curve2d::Line(Line2d {
        origin: DVec2::new(x, 0.0),
        direction: DVec2::new(0.0, 1.0),
    });
    let x_line = |y: f64| Curve2d::Line(Line2d {
        origin: DVec2::new(0.0, y),
        direction: DVec2::new(1.0, 0.0),
    });
    let lumin = if !umininf { Some(y_line(u_min)) } else { None };
    let lumax = if !umaxinf { Some(y_line(u_max)) } else { None };
    let lvmin = if !vmininf { Some(x_line(v_min)) } else { None };
    let lvmax = if !vmaxinf { Some(x_line(v_max)) } else { None };

    // OCCT L2695-2725: the iso 3d curves (hasiso =
    // S->IsKind(Geom_ElementarySurface)) with the Gabarit apex probes.
    let mut cumin: Option<Curve3> = None;
    let mut cumax: Option<Curve3> = None;
    let mut cvmin: Option<Curve3> = None;
    let mut cvmax: Option<Curve3> = None;
    let tol_apex = 1e-5;
    let hasiso = s.is_elementary();
    if hasiso {
        if !umininf {
            cumin = Some(surface_uiso_gap(s, u_min));
        }
        if !umaxinf {
            cumax = Some(surface_uiso_gap(s, u_max));
        }
        if !vmininf {
            cvmin = Some(surface_viso_gap(s, v_min));
            if gabarit(cvmin.as_ref().expect("viso")) <= tol_apex {
                vmindegen = true;
            }
        }
        if !vmaxinf {
            cvmax = Some(surface_viso_gap(s, v_max));
            if gabarit(cvmax.as_ref().expect("viso")) <= tol_apex {
                vmaxdegen = true;
            }
        }
    }

    // make the face
    // OCCT L2728: B.MakeFace(F, S, tol).
    *f = b_make_face_surface(s, tol);

    // make the edges
    // OCCT L2731-2863.
    let mut eumin = Shape::null();
    let mut eumax = Shape::null();
    let mut evmin = Shape::null();
    let mut evmax = Shape::null();

    if !umininf {
        if hasiso {
            eumin = b_make_edge_curve(&cumin.expect("cumin"), tol);
        } else {
            eumin = bat::builder_make_edge();
        }
        if uclosed {
            // OCCT L2745: B.UpdateEdge(eumin, Lumax, Lumin, F, tol) — the
            // seam two-pcurve form; the rcad store keeps one pcurve per
            // face key (architecture difference #28).
            bat::builder_update_edge_pcurve(&mut eumin, lumax.as_ref().expect("lumax"), f, tol);
        } else {
            bat::builder_update_edge_pcurve(&mut eumin, lumin.as_ref().expect("lumin"), f, tol);
        }
        if !vmininf {
            v00.orientation = Orientation::Forward;
            bat::builder_add_edge_vertex(&mut eumin, &v00);
        }
        if !vmaxinf {
            v01.orientation = Orientation::Reversed;
            bat::builder_add_edge_vertex(&mut eumin, &v01);
        }
        bat::builder_range_edge(&mut eumin, v_min, v_max);
    }

    if !umaxinf {
        if uclosed {
            eumax = eumin.clone();
        } else {
            if hasiso {
                eumax = b_make_edge_curve(&cumax.expect("cumax"), tol);
            } else {
                eumax = bat::builder_make_edge();
            }
            bat::builder_update_edge_pcurve(&mut eumax, lumax.as_ref().expect("lumax"), f, tol);
            if !vmininf {
                v10.orientation = Orientation::Forward;
                bat::builder_add_edge_vertex(&mut eumax, &v10);
            }
            if !vmaxinf {
                v11.orientation = Orientation::Reversed;
                bat::builder_add_edge_vertex(&mut eumax, &v11);
            }
            bat::builder_range_edge(&mut eumax, v_min, v_max);
        }
    }

    if !vmininf {
        if hasiso && !vmindegen {
            evmin = b_make_edge_curve(&cvmin.expect("cvmin"), tol);
        } else {
            evmin = bat::builder_make_edge();
        }
        if vclosed {
            // OCCT L2807: the seam two-pcurve form (see #28).
            bat::builder_update_edge_pcurve(&mut evmin, lvmin.as_ref().expect("lvmin"), f, tol);
        } else {
            bat::builder_update_edge_pcurve(&mut evmin, lvmin.as_ref().expect("lvmin"), f, tol);
        }
        if !umininf {
            v00.orientation = Orientation::Forward;
            bat::builder_add_edge_vertex(&mut evmin, &v00);
        }
        if !umaxinf {
            v10.orientation = Orientation::Reversed;
            bat::builder_add_edge_vertex(&mut evmin, &v10);
        }
        bat::builder_range_edge(&mut evmin, u_min, u_max);
        if vmindegen {
            bat::builder_set_degenerated(&mut evmin, true);
        }
    }

    if !vmaxinf {
        if vclosed {
            evmax = evmin.clone();
        } else {
            if hasiso && !vmaxdegen {
                evmax = b_make_edge_curve(&cvmax.expect("cvmax"), tol);
            } else {
                evmax = bat::builder_make_edge();
            }
            bat::builder_update_edge_pcurve(&mut evmax, lvmax.as_ref().expect("lvmax"), f, tol);
            if !umininf {
                v01.orientation = Orientation::Forward;
                bat::builder_add_edge_vertex(&mut evmax, &v01);
            }
            if !umaxinf {
                v11.orientation = Orientation::Reversed;
                bat::builder_add_edge_vertex(&mut evmax, &v11);
            }
            bat::builder_range_edge(&mut evmax, u_min, u_max);
            if vmaxdegen {
                bat::builder_set_degenerated(&mut evmax, true);
            }
        }
    }

    // make the wires and add them to the face
    // OCCT L2866-2918.
    eumin.orientation = Orientation::Reversed;
    evmax.orientation = Orientation::Reversed;

    if !umininf && !umaxinf && vmininf && vmaxinf {
        // two wires in u
        let mut w = bat::builder_make_wire();
        bat::builder_add_wire_edge(&mut w, &eumin);
        bat::builder_add_face_wire(f, &w);
        let mut w = bat::builder_make_wire();
        bat::builder_add_wire_edge(&mut w, &eumax);
        bat::builder_add_face_wire(f, &w);
        bat::builder_set_closed(f, uclosed);
    } else if umininf && umaxinf && !vmininf && !vmaxinf {
        // two wires in v
        let mut w = bat::builder_make_wire();
        bat::builder_add_wire_edge(&mut w, &evmin);
        bat::builder_add_face_wire(f, &w);
        let mut w = bat::builder_make_wire();
        bat::builder_add_wire_edge(&mut w, &evmax);
        bat::builder_add_face_wire(f, &w);
        bat::builder_set_closed(f, vclosed);
    } else if !umininf || !umaxinf || !vmininf || !vmaxinf {
        // one wire
        let mut w = bat::builder_make_wire();
        if !umininf {
            bat::builder_add_wire_edge(&mut w, &eumin);
        }
        if !vmininf {
            bat::builder_add_wire_edge(&mut w, &evmin);
        }
        if !umaxinf {
            bat::builder_add_wire_edge(&mut w, &eumax);
        }
        if !vmaxinf {
            bat::builder_add_wire_edge(&mut w, &evmax);
        }
        bat::builder_add_face_wire(f, &w);
        bat::builder_set_closed(&mut w, !umininf && !umaxinf && !vmininf && !vmaxinf);
        bat::builder_set_closed(f, uclosed && vclosed);
    }
}

/// OCCT BRep_Builder::MakeEdge(E, C, Tol) — the curve-carrying edge.
fn b_make_edge_curve(c: &Curve3, tol: f64) -> Shape {
    let mut e = bat::builder_make_edge();
    if let TShape::Edge(ed) = std::sync::Arc::make_mut(&mut e.data) {
        ed.curve = Some(c.clone());
        ed.tolerance = tol;
    }
    e
}

/// OCCT Geom_Surface::UIso(U) — GAP carrier (the iso-curve construction is
/// not translated; the brep_offset_offset_b.rs #1608 precedent).
fn surface_uiso_gap(s: &Surface3, _u: f64) -> Curve3 {
    let _ = s;
    panic!("GAP: Geom_Surface::UIso (iso-curve construction not translated)");
}

/// OCCT Geom_Surface::VIso(V) — GAP carrier (idem).
fn surface_viso_gap(s: &Surface3, _v: f64) -> Curve3 {
    let _ = s;
    panic!("GAP: Geom_Surface::VIso (iso-curve construction not translated)");
}

// ---------------------------------------------------------------------------
// OCCT static EnlargeGeometry (cxx L2923-3222).
// ---------------------------------------------------------------------------

/// OCCT static EnlargeGeometry(S, U1, U2, V1, V2, IsV1degen, IsV2degen,
/// uf1, uf2, vf1, vf2, coeff, theGlobalEnlargeU, theGlobalEnlargeVfirst,
/// theGlobalEnlargeVlast, theLenBeforeUfirst, theLenAfterUlast,
/// theLenBeforeVfirst, theLenAfterVlast) (cxx L2923-3222).
#[allow(clippy::too_many_arguments)]
fn enlarge_geometry(
    s: &mut Surface3,
    u1: &mut f64,
    u2: &mut f64,
    v1: &mut f64,
    v2: &mut f64,
    is_v1degen: &mut bool,
    is_v2degen: &mut bool,
    uf1: f64,
    uf2: f64,
    vf1: f64,
    vf2: f64,
    coeff: f64,
    the_global_enlarge_u: bool,
    the_global_enlarge_vfirst: bool,
    the_global_enlarge_vlast: bool,
    the_len_before_ufirst: f64,
    the_len_after_ulast: f64,
    the_len_before_vfirst: f64,
    the_len_after_vlast: f64,
) -> bool {
    let tol_apex = 1e-5;

    // OCCT L2945-2987: the RectangularTrimmedSurface branch (the recursive
    // enlarge + the partial retrim).
    let mut surface_change = false;
    if let Surface3::Trimmed(ts) = s {
        let mut bs = (*ts.basis).clone();
        surface_change = enlarge_geometry(
            &mut bs,
            u1,
            u2,
            v1,
            v2,
            is_v1degen,
            is_v2degen,
            uf1,
            uf2,
            vf1,
            vf2,
            coeff,
            the_global_enlarge_u,
            the_global_enlarge_vfirst,
            the_global_enlarge_vlast,
            the_len_before_ufirst,
            the_len_after_ulast,
            the_len_before_vfirst,
            the_len_after_vlast,
        );
        if !the_global_enlarge_vfirst {
            *v1 = vf1;
        }
        if !the_global_enlarge_vlast {
            *v2 = vf2;
        }
        if !the_global_enlarge_vfirst || !the_global_enlarge_vlast {
            *s = Surface3::Trimmed(rcad_kernel::geom::TrimmedSurface::new(
                bs.clone(),
                *u1,
                *u2,
                *v1,
                *v2,
            ));
        } else {
            *s = bs.clone();
        }
        surface_change = true;
        return surface_change;
    }
    // OCCT L2988-3011: the OffsetSurface branch (the recursive basis
    // enlarge + SetBasisSurface).
    if let Surface3::Offset(os) = s {
        let mut surf = (*os.basis).clone();
        surface_change = enlarge_geometry(
            &mut surf,
            u1,
            u2,
            v1,
            v2,
            is_v1degen,
            is_v2degen,
            uf1,
            uf2,
            vf1,
            vf2,
            coeff,
            the_global_enlarge_u,
            the_global_enlarge_vfirst,
            the_global_enlarge_vlast,
            the_len_before_ufirst,
            the_len_after_ulast,
            the_len_before_vfirst,
            the_len_after_vlast,
        );
        os.basis = Box::new(surf);
        return surface_change;
    }
    // OCCT L3012-3108: the SurfaceOfLinearExtrusion / SurfaceOfRevolution
    // branch.
    if matches!(s, Surface3::LinearExtrusion(_) | Surface3::Revolution(_)) {
        let mut enlarge_u = the_global_enlarge_u;
        let mut enlarge_v = true;
        let mut enlarge_ufirst = enlarge_u;
        let mut enlarge_ulast = enlarge_u;
        let mut enlarge_vfirst = the_global_enlarge_vfirst;
        let mut enlarge_vlast = the_global_enlarge_vlast;
        let (mut su1, mut su2, mut sv1, mut sv2) = {
            let d = s.default_domain();
            (d[0], d[1], d[2], d[3])
        };
        if precision_is_infinite(su1) || precision_is_infinite(su2) {
            let du_first = uf2 - uf1;
            let du_last = uf2 - uf1;
            su1 = uf1 - du_first;
            su2 = uf2 + du_last;
            enlarge_u = false;
        } else if s.is_u_closed() {
            enlarge_u = false;
        } else {
            // OCCT L3035-3050: the viso Gabarit probes — GAP leaves (the
            // iso-curve construction; annotated above).
            let viso_gap = surface_viso_gap(s, vf1);
            let _du_default = gcpnts_length_gap(&viso_gap) * coeff;
            let _ = (&mut enlarge_ufirst, &mut enlarge_ulast);
        }
        if precision_is_infinite(sv1) || precision_is_infinite(sv2) {
            let dv_first = vf2 - vf1;
            let dv_last = vf2 - vf1;
            sv1 = vf1 - dv_first;
            sv2 = vf2 + dv_last;
            enlarge_v = false;
        } else if s.is_v_closed() {
            enlarge_v = false;
        } else {
            // OCCT L3064-3080: the uiso/viso Gabarit probes — GAP leaves.
            let uiso_gap = surface_uiso_gap(s, uf1);
            let _dv_default = gcpnts_length_gap(&uiso_gap) * coeff;
            let viso1_gap = surface_viso_gap(s, vf1);
            let viso2_gap = surface_viso_gap(s, vf2);
            if gabarit(&viso1_gap) <= tol_apex {
                enlarge_vfirst = false;
                *is_v1degen = true;
            }
            if gabarit(&viso2_gap) <= tol_apex {
                enlarge_vlast = false;
                *is_v2degen = true;
            }
        }
        // OCCT L3082-3107: the ExtendSurfByLength forms — GAP leaves (the
        // GeomLib batch); the OCCT Bounds write closes the branch.
        let _ = (enlarge_u, enlarge_v, enlarge_ufirst, enlarge_ulast, enlarge_vfirst, enlarge_vlast);
        let _ = (the_len_before_ufirst, the_len_after_ulast, the_len_before_vfirst, the_len_after_vlast);
        panic!("GAP: GeomLib::ExtendSurfByLength (TKTopAlgo/GeomLib not translated)");
    }
    // OCCT L3109-3210: the Bezier/BSpline branch.
    if matches!(s, Surface3::Bezier(_) | Surface3::BSpline(_)) {
        let mut enlarge_u = the_global_enlarge_u;
        let mut enlarge_v = true;
        let mut enlarge_ufirst = enlarge_u;
        let mut enlarge_ulast = enlarge_u;
        let mut enlarge_vfirst = the_global_enlarge_vfirst;
        let mut enlarge_vlast = the_global_enlarge_vlast;
        if s.is_u_closed() {
            enlarge_u = false;
        }
        if s.is_v_closed() {
            enlarge_v = false;
        }

        let duf = uf2 - uf1;
        let dvf = vf2 - vf1;
        let (su1, su2, sv1, sv2) = {
            let d = s.default_domain();
            (d[0], d[1], d[2], d[3])
        };

        // OCCT L3132-3147: the iso Gabarits — GAP leaves (the iso-curve
        // construction; annotated above).
        let uiso1_gap = surface_uiso_gap(s, su1);
        let uiso2_gap = surface_uiso_gap(s, su2);
        let viso1_gap = surface_viso_gap(s, sv1);
        let viso2_gap = surface_viso_gap(s, sv2);
        let gabarit_uiso1 = gabarit(&uiso1_gap);
        let gabarit_uiso2 = gabarit(&uiso2_gap);
        let gabarit_viso1 = gabarit(&viso1_gap);
        let gabarit_viso2 = gabarit(&viso2_gap);
        if gabarit_viso1 <= tol_apex || gabarit_viso2 <= tol_apex {
            enlarge_u = false;
        }
        if gabarit_uiso1 <= tol_apex || gabarit_uiso2 <= tol_apex {
            enlarge_v = false;
        }

        let _ = (gabarit_uiso1, gabarit_uiso2);
        if enlarge_u {
            let _du_default = gcpnts_length_gap(&viso1_gap) * coeff;
            let _ = (duf, the_len_before_ufirst, the_len_after_ulast);
            if gabarit_uiso1 <= tol_apex {
                enlarge_ufirst = false;
            }
            if gabarit_uiso2 <= tol_apex {
                enlarge_ulast = false;
            }
        }
        if enlarge_v {
            let _dv_default = gcpnts_length_gap(&uiso1_gap) * coeff;
            let _ = (dvf, the_len_before_vfirst, the_len_after_vlast);
            if gabarit_viso1 <= tol_apex {
                enlarge_vfirst = false;
                *is_v1degen = true;
            }
            if gabarit_viso2 <= tol_apex {
                enlarge_vlast = false;
                *is_v2degen = true;
            }
        }
        // OCCT L3183-3209: the ExtendSurfByLength forms — GAP leaf.
        let _ = (enlarge_u, enlarge_v, enlarge_ufirst, enlarge_ulast, enlarge_vfirst, enlarge_vlast);
        panic!("GAP: GeomLib::ExtendSurfByLength (TKTopAlgo/GeomLib not translated)");
    }
    // OCCT L3211-3220: the remaining types — the bounds clamp.
    let (uu1, uu2, vv1, vv2) = {
        let d = s.default_domain();
        (d[0], d[1], d[2], d[3])
    };
    // Pas d extension au dela des bornes de la surface.
    *u1 = uu1.max(*u1);
    *v1 = vv1.max(*v1);
    *u2 = uu2.min(*u2);
    *v2 = vv2.min(*v2);
    surface_change
}

/// OCCT GCPnts_AbscissaPoint::Length(C) — GAP carrier (architecture
/// difference #24).
fn gcpnts_length_gap(_c: &Curve3) -> f64 {
    panic!("GAP: GCPnts_AbscissaPoint::Length (TKGeomAlgo not translated)");
}

// ---------------------------------------------------------------------------
// OCCT static UpdatePCurves (cxx L3230-3263).
// ---------------------------------------------------------------------------

/// OCCT static UpdatePCurves(F, BF) (cxx L3230-3263) — the pcurve copy of F
/// onto the enlarged face BF (F and BF have to be FORWARD).
fn update_pcurves(f: &Shape, bf: &mut Shape) {
    // OCCT L3235-3238: Emap; NullPCurve.
    let mut emap = OcctIndexedShapeMap::new();
    for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        emap.add(&e);
    }

    for i in 1..=emap.extent() {
        let mut ce = emap.at_1(i).clone();
        ce.orientation = Orientation::Forward;
        // OCCT L3244: C2 = BRep_Tool::CurveOnSurface(CE, F, f, l).
        let (c2, fp, lp) = match brep_tool_curve_on_surface(&ce, f) {
            Some(t) => t,
            None => continue,
        };
        if bat::brep_tool_is_closed_on_surface(&ce, f) {
            // OCCT L3247-3253: the seam form (CE.Reverse(); C2R; the
            // NullPCurve resets; the two-pcurve rebind) — the rcad store
            // keeps one pcurve per face key (architecture difference #28).
            let mut ce_rev = ce.clone();
            ce_rev.orientation = Orientation::Reversed;
            let c2r = brep_tool_curve_on_surface(&ce_rev, f).map(|(c, _, _)| c);
            b_remove_edge_pcurve(&mut ce, f);
            let ce_tol = brep_tool_tolerance(&ce);
            bat::builder_update_edge_pcurve(&mut ce, &c2, bf, ce_tol);
            let _ = (c2r, ce_rev);
        } else {
            // OCCT L3254-3258: the simple form.
            b_remove_edge_pcurve(&mut ce, f);
            let ce_tol = brep_tool_tolerance(&ce);
            bat::builder_update_edge_pcurve(&mut ce, &c2, bf, ce_tol);
        }
        // OCCT L3260: B.Range(CE, f, l).
        bat::builder_range_edge(&mut ce, fp, lp);
    }
}

// ---------------------------------------------------------------------------
// OCCT static CompactUVBounds (cxx L3267-3305).
// ---------------------------------------------------------------------------

/// OCCT static CompactUVBounds(F, UMin, UMax, VMin, VMax) (cxx L3267-3305)
/// — Calcul serre pour que les bornes ne couvrent pas plus d une periode.
pub(crate) fn compact_uv_bounds(f: &Shape, u_min: &mut f64, u_max: &mut f64, v_min: &mut f64, v_max: &mut f64) {
    let n = 33usize;
    let mut b = BndBox2d::new();

    for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L3282-3283: BRepAdaptor_Curve2d C(E, F); BRep_Tool::Range(E, U1, U2).
        let (c, u1, u2) = match brep_tool_curve_on_surface(&e, f) {
            Some(t) => t,
            None => continue,
        };
        let mut u = u1;
        let du = (u2 - u1) / (n as f64 - 1.0);
        for _ in 1..n {
            let p = Curve2dEval::point_at(&c, u);
            u += du;
            b.update(p.x, p.y, p.x, p.y);
        }
        let p = Curve2dEval::point_at(&c, u2);
        b.update(p.x, p.y, p.x, p.y);
    }

    // OCCT L3297-3304.
    if let Some((g0, g1, g2, g3)) = b.get() {
        *u_min = g0;
        *v_min = g1;
        *u_max = g2;
        *v_max = g3;
    } else if let Some(s) = face_surface_of(f) {
        let d = s.default_domain();
        *u_min = d[0];
        *u_max = d[1];
        *v_min = d[2];
        *v_max = d[3];
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::CheckBounds (hxx L123-127; cxx L3309-3433).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::CheckBounds(F, Analyse, enlargeU, enlargeVfirst,
/// enlargeVlast) (cxx L3309-3433).
pub fn check_bounds(
    f: &Shape,
    analyse: &BRepOffsetAnalyse,
    enlarge_u: &mut bool,
    enlarge_vfirst: &mut bool,
    enlarge_vlast: &mut bool,
) {
    *enlarge_u = true;
    *enlarge_vfirst = true;
    *enlarge_vlast = true;

    let mut u_bound = 0i32;
    let mut v_bound = 0i32;
    let mut ufirst = f64::MAX;
    let mut ulast = f64::MIN;
    let mut vfirst = f64::MAX;
    let mut vlast = f64::MIN;

    // OCCT L3323-3324: CompactUVBounds(F, UF1, UF2, VF1, VF2).
    let (mut uf1, mut uf2, mut vf1, mut vf2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    compact_uv_bounds(f, &mut uf1, &mut uf2, &mut vf1, &mut vf2);

    // OCCT L3326-3330: the basis-surface collapse.
    let mut the_surf = face_surface_of(f);
    if let Some(Surface3::Trimmed(ts)) = the_surf {
        the_surf = Some((*ts.basis).clone());
    }

    // OCCT L3332-3336: the extrusion / revolution / Bezier / BSpline probe.
    let is_swept_or_freeform = matches!(
        &the_surf,
        Some(Surface3::LinearExtrusion(_))
            | Some(Surface3::Revolution(_))
            | Some(Surface3::Bezier(_))
            | Some(Surface3::BSpline(_))
    );
    if is_swept_or_freeform {
        for an_edge in explorer(f, ShapeType::Edge, ShapeType::Shape) {
            // OCCT L3340-3342: L = Analyse.Type(anEdge).
            let l = analyse.type_(&an_edge);
            if !l.is_empty() || brep_tool_degenerated(&an_edge) {
                // OCCT L3344-3345: OT = L.First().Type().
                let ot = l.first().map(|iv| iv.my_type);
                let ot_tangential = ot == Some(crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Tangential);
                if ot_tangential || brep_tool_degenerated(&an_edge) {
                    // OCCT L3347-3349: aCurve = CurveOnSurface with the
                    // trimmed basis extraction.
                    let (mut a_curve, fpar, lpar) = match brep_tool_curve_on_surface(&an_edge, f) {
                        Some(t) => t,
                        None => continue,
                    };
                    if let Curve2d::Trimmed(_) = a_curve {
                        a_curve = basis_curve2(&a_curve);
                    }

                    // OCCT L3354-3370: theLine = the Line/Bezier/BSpline
                    // ConvertToLine2d form — GAP leaf (architecture
                    // difference #24); the carrier keeps the null-line
                    // fall-through of the OCCT (annotated).
                    let the_line: Option<Line2d> = match &a_curve {
                        Curve2d::Line(l) => Some(l.clone()),
                        _ => None, // ShapeCustom_Curve2d::ConvertToLine2d GAP
                    };

                    if let Some(line) = the_line {
                        // OCCT L3374-3412: the DX2d / DY2d parallel probes.
                        let dir = line.direction;
                        let dx2d = DVec2::new(1.0, 0.0);
                        let dy2d = DVec2::new(0.0, 1.0);
                        let angular = rcad_kernel::precision::ANGULAR;
                        let is_par = |a: DVec2, b: DVec2| {
                            (a.x * b.y - a.y * b.x).abs() <= angular
                        };
                        if is_par(dir, dx2d) {
                            v_bound += 1;
                            if brep_tool_degenerated(&an_edge) {
                                if (line.origin.y - vf1).abs() <= rcad_kernel::precision::CONFUSION {
                                    *enlarge_vfirst = false;
                                } else {
                                    // theLine->Location().Y() is near VF2
                                    *enlarge_vlast = false;
                                }
                            } else {
                                if line.origin.y < vfirst {
                                    vfirst = line.origin.y;
                                }
                                if line.origin.y > vlast {
                                    vlast = line.origin.y;
                                }
                            }
                        } else if is_par(dir, dy2d) {
                            u_bound += 1;
                            if line.origin.x < ufirst {
                                ufirst = line.origin.x;
                            }
                            if line.origin.x > ulast {
                                ulast = line.origin.x;
                            }
                        }
                    }
                    let _ = (fpar, lpar);
                }
            }
        }
    }

    // OCCT L3419-3432: the bound-count decision.
    if u_bound >= 2 || v_bound >= 2 {
        if u_bound >= 2
            && (uf1 - ufirst).abs() <= rcad_kernel::precision::CONFUSION
            && (uf2 - ulast).abs() <= rcad_kernel::precision::CONFUSION
        {
            *enlarge_u = false;
        }
        if v_bound >= 2
            && (vf1 - vfirst).abs() <= rcad_kernel::precision::CONFUSION
            && (vf2 - vlast).abs() <= rcad_kernel::precision::CONFUSION
        {
            *enlarge_vfirst = false;
            *enlarge_vlast = false;
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::EnLargeFace (hxx L145-156; cxx L3437-3634).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::EnLargeFace(F, BF, CanExtentSurface, UpdatePCurve,
/// theEnlargeU, theEnlargeVfirst, theEnlargeVlast, theExtensionMode,
/// theLenBeforeUfirst, theLenAfterUlast, theLenBeforeVfirst,
/// theLenAfterVlast) (cxx L3437-3634) — returns true if the surface of NF
/// has changed.
#[allow(clippy::too_many_arguments)]
pub fn en_large_face(
    f: &Shape,
    bf: &mut Shape,
    can_extent_surface: bool,
    update_pcurve: bool,
    the_enlarge_u: bool,
    the_enlarge_vfirst: bool,
    the_enlarge_vlast: bool,
    the_extension_mode: i32,
    the_len_before_ufirst: f64,
    the_len_after_ulast: f64,
    the_len_before_vfirst: f64,
    the_len_after_vlast: f64,
) -> bool {
    //---------------------------
    // extension de la geometrie.
    //---------------------------
    // OCCT L3453-3460: S = the face surface; the uperiodic/vperiodic and
    // degeneracy flags.
    let mut s = match face_surface_of(f) {
        Some(v) => v,
        None => return false,
    };
    let mut uperiodic = false;
    let mut vperiodic = false;
    let mut is_vv1degen = false;
    let mut is_vv2degen = false;
    let mut surface_change = false;

    // OCCT L3462-3470: CompactUVBounds / BRepTools::UVBounds.
    let (uf1, uf2, vf1, vf2);
    if s.is_u_periodic() || s.is_v_periodic() {
        // Calcul serre pour que les bornes ne couvre pas plus d une periode
        let (mut a_uf1, mut a_uf2, mut a_vf1, mut a_vf2) =
            (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        compact_uv_bounds(f, &mut a_uf1, &mut a_uf2, &mut a_vf1, &mut a_vf2);
        uf1 = a_uf1;
        uf2 = a_uf2;
        vf1 = a_vf1;
        vf2 = a_vf2;
    } else {
        let b = crate::feat::loc_ope_generator_b::brep_tools_uv_bounds(f);
        uf1 = b[0];
        uf2 = b[1];
        vf1 = b[2];
        vf2 = b[3];
    }

    // OCCT L3472: S->Bounds(US1, US2, VS1, VS2).
    let s_domain = s.default_domain();
    let (us1, us2, vs1, vs2) = (s_domain[0], s_domain[1], s_domain[2], s_domain[3]);
    // OCCT L3473-3489: coeff and the UU/VV targets.
    let (mut uu1, mut uu2, mut vv1, mut vv2, coeff) = if the_extension_mode == 1 {
        (
            -THE_INFINI,
            THE_INFINI,
            -THE_INFINI,
            THE_INFINI,
            0.25f64,
        )
    } else {
        let face_du = uf2 - uf1;
        let face_dv = vf2 - vf1;
        (
            uf1 - 10.0 * face_du,
            uf2 + 10.0 * face_du,
            vf1 - 10.0 * face_dv,
            vf2 + 10.0 * face_dv,
            1.0f64,
        )
    };

    // OCCT L3491-3519: EnlargeGeometry or the bounds clamp.
    if can_extent_surface {
        surface_change = enlarge_geometry(
            &mut s,
            &mut uu1,
            &mut uu2,
            &mut vv1,
            &mut vv2,
            &mut is_vv1degen,
            &mut is_vv2degen,
            uf1,
            uf2,
            vf1,
            vf2,
            coeff,
            the_enlarge_u,
            the_enlarge_vfirst,
            the_enlarge_vlast,
            the_len_before_ufirst,
            the_len_after_ulast,
            the_len_before_vfirst,
            the_len_after_vlast,
        );
    } else {
        uu1 = us1.max(uu1);
        uu2 = uu2.min(us2);
        vv1 = vs1.max(vv1);
        vv2 = vv2.min(vs2);
    }

    // OCCT L3521-3546: the periodic bounds.
    if s.is_u_periodic() {
        uperiodic = true;
        let period = surface_u_period(&s);
        let delta = period - (uf2 - uf1);
        let alpha = 0.1;
        uu1 = uf1 - alpha * delta;
        uu2 = uf2 + alpha * delta;
        if (uu2 - uu1) > period {
            uu2 = uu1 + period;
        }
    }
    if s.is_v_periodic() {
        vperiodic = true;
        let period = surface_v_period(&s);
        let delta = period - (vf2 - vf1);
        let alpha = 0.1;
        vv1 = vf1 - alpha * delta;
        vv2 = vf2 + alpha * delta;
        if (vv2 - vv1) > period {
            vv2 = vv1 + period;
        }
    }

    // Special treatment for conical surfaces
    // OCCT L3548-3574: the cone apex clamp.
    let mut the_surf = s.clone();
    if let Surface3::Trimmed(ts) = &the_surf {
        the_surf = (*ts.basis).clone();
    }
    if let Surface3::Cone(cone) = &the_surf {
        let the_apex = cone.apex;
        let y_dir = cone.axis.cross(cone.ref_dir);
        let (_u_apex, v_apex) = elslib_cone_parameters(
            the_apex,
            cone.apex,
            cone.ref_dir,
            y_dir,
            cone.axis,
            cone.radius,
            cone.half_angle_rad,
        );
        if vv1 < v_apex && v_apex < vv2 {
            // consider that VF1 and VF2 are on the same side from apex
            let tol_apex = 1e-5;
            if v_apex - vf1 >= tol_apex || v_apex - vf2 >= tol_apex {
                vv2 = v_apex;
            } else {
                vv1 = v_apex;
            }
        }
    }

    // OCCT L3576-3588: the enlarge-flag clamps.
    if !the_enlarge_u {
        uu1 = uf1;
        uu2 = uf2;
    }
    if !the_enlarge_vfirst {
        vv1 = vf1;
    }
    if !the_enlarge_vlast {
        vv2 = vf2;
    }

    // Detect closedness in U and V directions
    // OCCT L3590-3600: BRepTools::DetectClosedness — the reduced re-host
    // probes the surface closedness (architecture difference #31).
    let mut uclosed = s.is_u_closed();
    let mut vclosed = s.is_v_closed();
    if uclosed && !uperiodic && (the_len_before_ufirst != 0.0 || the_len_after_ulast != 0.0) {
        uclosed = false;
    }
    if vclosed && !vperiodic && (the_len_before_vfirst != 0.0 && the_len_after_vlast != 0.0) {
        vclosed = false;
    }

    // OCCT L3602-3603: MakeFace(...); BF.Location(L) — the rcad storage is
    // location-baked (identity; annotated).
    make_face(
        &s,
        uu1,
        uu2,
        vv1,
        vv2,
        uclosed,
        vclosed,
        is_vv1degen,
        is_vv2degen,
        bf,
    );

    // OCCT L3623-3630: the pcurve / surface updates.
    if surface_change && update_pcurve {
        let mut f_forward = f.clone();
        f_forward.orientation = Orientation::Forward;
        update_pcurves(&f_forward, bf);
        // OCCT L3629: BB.UpdateFace(F, S, L, BRep_Tool::Tolerance(F)) — the
        // rcad caller holds a shared handle; the in-place surface update of
        // the argument face is staged with the owner batch (annotated).
        let mut f_updated = f.clone();
        b_update_face_surface(&mut f_updated, &s, brep_tool_tolerance(f));
    }

    // OCCT L3632-3633.
    bf.orientation = f.orientation;
    surface_change
}

// ---------------------------------------------------------------------------
// OCCT static TryParameter (cxx L3638-3681).
// ---------------------------------------------------------------------------

/// OCCT static TryParameter(OE, V, NE, TolConf) (cxx L3638-3681).
fn try_parameter(oe: &Shape, v: &mut Shape, ne: &Shape, tol_conf: f64) -> bool {
    // OCCT L3643-3650: the (Of, Ol, Nf, Nl, U, P, OK) locals.
    let mut ok = false;
    let (oc, of, ol) = match brep_tool_curve(oe) {
        Some(t) => t,
        None => return false,
    };
    let (nc, nf, nl) = match brep_tool_curve(ne) {
        Some(t) => t,
        None => return false,
    };
    let mut u = 0.0f64;
    let p = match bat::brep_tool_pnt(v) {
        Some(t) => t,
        None => return false,
    };

    // OCCT L3653-3660: the Of probe.
    if p.distance(CurveEval::point_at(&oc, of)) < tol_conf {
        if of > nf && of < nl && p.distance(CurveEval::point_at(&nc, of)) < tol_conf {
            ok = true;
            u = of;
        }
    }
    // OCCT L3661-3668: the Ol probe.
    if p.distance(CurveEval::point_at(&oc, ol)) < tol_conf {
        if ol > nf && ol < nl && p.distance(CurveEval::point_at(&nc, ol)) < tol_conf {
            ok = true;
            u = ol;
        }
    }
    // OCCT L3669-3679: the B.UpdateVertex form.
    if ok {
        let mut ne_fwd = ne.clone();
        ne_fwd.orientation = Orientation::Forward;
        let mut v_int = v.clone();
        v_int.orientation = Orientation::Internal;
        b_update_vertex_on_edge(&mut ne_fwd, &v_int, u, brep_tool_tolerance(ne));
        v.orientation = v.orientation;
    }
    // OCCT L3680: return OK.
    ok
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::MapVertexEdges (hxx L180-183; cxx L3685-3716).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::MapVertexEdges(S, MEV) (cxx L3685-3716) — store in
/// MVE for a vertex V in S the incident edges E in S (one entry per edge).
pub fn map_vertex_edges(
    s: &Shape,
    mev: &mut ShapeDataMap<Vec<Shape>>,
) {
    // OCCT L3690-3696: the FORWARD-ordered edge walk with the DejaVu fence.
    let mut deja_vu: OcctShapeSet = HashMap::new();
    for e in explorer(&oriented(s, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
        if set_add(&mut deja_vu, &e) {
            let (v1, v2) = top_exp_vertices(&e);
            // OCCT L3699-3704.
            if !shape_data_map::is_bound(mev, &v1) {
                shape_data_map::bind(mev, &v1, Vec::new());
            }
            shape_data_map::change_find(mev, &v1).push(e.clone());
            // OCCT L3705-3713.
            if !v1.is_same(&v2) {
                if !shape_data_map::is_bound(mev, &v2) {
                    shape_data_map::bind(mev, &v2, Vec::new());
                }
                shape_data_map::change_find(mev, &v2).push(e);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::BuildNeighbour (hxx L171-175; cxx L3720-3789).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::BuildNeighbour(W, F, NOnV1, NOnV2) (cxx
/// L3720-3789) — via the wire explorer store in NOnV1/NOnV2 the edge
/// neighbours on the extremity vertices.
pub fn build_neighbour(
    w: &Shape,
    f: &Shape,
    n_on_v1: &mut ShapeDataMap<Shape>,
    n_on_v2: &mut ShapeDataMap<Shape>,
) {
    // OCCT L3728-3734: wexp.Init(the FORWARD wire, the FORWARD face) — the
    // rcad sub-shape order (the BRepTools_WireExplorer re-host; annotated).
    let _ = f;
    let mut cur_e = Shape::null();
    let mut first_e = Shape::null();
    let mut prec_e = Shape::null();
    let mut v1 = Shape::null();
    let mut v2 = Shape::null();
    let mut vp1 = Shape::null();
    let mut vp2 = Shape::null();
    let mut fv1 = Shape::null();
    let mut fv2 = Shape::null();

    let wire_edges = bat::sub_shapes(&oriented(w, Orientation::Forward));
    let mut idx = 0usize;
    if let Some(e0) = wire_edges.first() {
        cur_e = e0.clone();
        first_e = e0.clone();
        prec_e = e0.clone();
        let (t1, t2) = top_exp_vertices(&cur_e);
        v1 = t1;
        v2 = t2;
        fv1 = vp1_or(&v1);
        fv2 = vp2_or(&v2);
        let _ = (&fv1, &fv2);
        fv1 = v1.clone();
        fv2 = v2.clone();
        vp1 = v1.clone();
        vp2 = v2.clone();
    }
    idx = 1;
    // OCCT L3740-3768: while (wexp.More()).
    while idx < wire_edges.len() {
        cur_e = wire_edges[idx].clone();
        let (t1, t2) = top_exp_vertices(&cur_e);
        v1 = t1;
        v2 = t2;
        if v1.is_same(&vp1) {
            shape_data_map::bind(n_on_v1, &prec_e, cur_e.clone());
            shape_data_map::bind(n_on_v1, &cur_e, prec_e.clone());
        }
        if v1.is_same(&vp2) {
            shape_data_map::bind(n_on_v2, &prec_e, cur_e.clone());
            shape_data_map::bind(n_on_v1, &cur_e, prec_e.clone());
        }
        if v2.is_same(&vp1) {
            shape_data_map::bind(n_on_v1, &prec_e, cur_e.clone());
            shape_data_map::bind(n_on_v2, &cur_e, prec_e.clone());
        }
        if v2.is_same(&vp2) {
            shape_data_map::bind(n_on_v2, &prec_e, cur_e.clone());
            shape_data_map::bind(n_on_v2, &cur_e, prec_e.clone());
        }
        prec_e = cur_e.clone();
        vp1 = v1.clone();
        vp2 = v2.clone();
        idx += 1;
    }
    // OCCT L3769-3788: the closure binds.
    if v1.is_same(&fv1) {
        shape_data_map::bind(n_on_v1, &first_e, cur_e.clone());
        shape_data_map::bind(n_on_v1, &cur_e, first_e.clone());
    }
    if v1.is_same(&fv2) {
        shape_data_map::bind(n_on_v2, &first_e, cur_e.clone());
        shape_data_map::bind(n_on_v1, &cur_e, first_e.clone());
    }
    if v2.is_same(&fv1) {
        shape_data_map::bind(n_on_v1, &first_e, cur_e.clone());
        shape_data_map::bind(n_on_v2, &cur_e, first_e.clone());
    }
    if v2.is_same(&fv2) {
        shape_data_map::bind(n_on_v2, &first_e, cur_e.clone());
        shape_data_map::bind(n_on_v2, &cur_e, first_e.clone());
    }
}

fn vp1_or(v: &Shape) -> Shape {
    v.clone()
}
fn vp2_or(v: &Shape) -> Shape {
    v.clone()
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::ExtentFace (hxx L158-164; cxx L3793-4333).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::ExtentFace(F, ConstShapes, ToBuild, Side, TolConf,
/// NF) (cxx L3793-4333).
pub fn extent_face(
    f: &Shape,
    const_shapes: &mut ShapeDataMap<Shape>,
    to_build: &mut ShapeDataMap<Shape>,
    side: State,
    tol_conf: f64,
    nf: &mut Shape,
) {
    // OCCT L3802-3807: Build; Extent; the edge locals; BRep_Builder B; EF.
    let mut build: ShapeDataMap<Shape> = HashMap::new();
    let mut ef = Shape::null();

    // Construction de la boite englobante de la face a etendre et des
    // bouchons pour limiter les extensions. (commented-out OCCT block)

    // OCCT L3818-3819: SurfaceChange = EnLargeFace(F, EF, true).
    let surface_change = en_large_face(f, &mut ef, true, false, true, true, true, 1, -1.0, -1.0, -1.0, -1.0);

    // OCCT L3821-3824: NF = EF.EmptyCopied(); NF.Orientation(FORWARD).
    *nf = empty_copied(&ef);
    nf.orientation = Orientation::Forward;

    if surface_change {
        //------------------------------------------------
        // Mise a jour des pcurves sur la surface de base.
        //------------------------------------------------
        // OCCT L3831-3876: the pcurve rebind walk.
        let mut f_forward = f.clone();
        f_forward.orientation = Orientation::Forward;
        let mut emap = OcctIndexedShapeMap::new();
        for e in explorer(&f_forward, ShapeType::Edge, ShapeType::Shape) {
            emap.add(&e);
        }
        for i in 1..=emap.extent() {
            let mut ce = emap.at_1(i).clone();
            ce.orientation = Orientation::Forward;
            // OCCT L3840-3847: Ecs = ConstShapes(CE) (the patch) and the
            // range read.
            let mut ecs = Shape::null();
            let (c2, fp, lp) = match brep_tool_curve_on_surface(&ce, &f_forward) {
                Some(t) => t,
                None => continue,
            };
            let (mut fp, mut lp) = (fp, lp);
            if shape_data_map::is_bound(const_shapes, &ce) {
                ecs = shape_data_map::find(const_shapes, &ce);
                let r = brep_tool_range(&ecs);
                fp = r.0;
                lp = r.1;
            }
            if bat::brep_tool_is_closed_on_surface(&ce, &f_forward) {
                // OCCT L3849-3860: the seam two-pcurve form (see #28).
                let ce_rev = oriented(&ce, Orientation::Reversed);
                let c2r = brep_tool_curve_on_surface(&ce_rev, &f_forward).map(|(c, _, _)| c);
                let ce_tol = brep_tool_tolerance(&ce);
                bat::builder_update_edge_pcurve(&mut ce, &c2, &ef, ce_tol);
                if !ecs.is_null() {
                    bat::builder_update_edge_pcurve(&mut ecs, &c2, &ef, ce_tol);
                }
                let _ = c2r;
            } else {
                // OCCT L3862-3869.
                let ce_tol = brep_tool_tolerance(&ce);
                bat::builder_update_edge_pcurve(&mut ce, &c2, &ef, ce_tol);
                if !ecs.is_null() {
                    bat::builder_update_edge_pcurve(&mut ecs, &c2, &ef, ce_tol);
                }
            }
            // OCCT L3870-3874.
            bat::builder_range_edge(&mut ce, fp, lp);
            if !ecs.is_null() {
                bat::builder_range_edge(&mut ecs, fp, lp);
            }
        }
    }

    // OCCT L3879-3968: the wire walk — Construction edges.
    for w in explorer(&oriented(f, Orientation::Forward), ShapeType::Wire, ShapeType::Shape) {
        // OCCT L3882-3888: MVE; NOnV1; NOnV2.
        let mut mve: ShapeDataMap<Vec<Shape>> = HashMap::new(); // Vertex -> Edges incidentes.
        let mut n_on_v1: ShapeDataMap<Shape> = HashMap::new();
        let mut n_on_v2: ShapeDataMap<Shape> = HashMap::new();

        map_vertex_edges(&w, &mut mve);
        build_neighbour(&w, f, &mut n_on_v1, &mut n_on_v2);

        let mut l_int1: Vec<Shape> = Vec::new();
        let mut l_int2: Vec<Shape> = Vec::new();
        let mut stop_face = Shape::null();
        //------------------------------------------------
        // Construction edges
        //------------------------------------------------
        // OCCT L3895-3918: the first edge pass (the TryProject forms).
        for e in explorer(&oriented(&w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
            if shape_data_map::is_bound(const_shapes, &e) {
                shape_data_map::un_bind(to_build, &e);
            }
            if shape_data_map::is_bound(to_build, &e) {
                let mut loe: Vec<Shape> = Vec::new();
                loe.push(e.clone());
                let tb_e = shape_data_map::find(to_build, &e);
                let ok = try_project(&tb_e, &ef, &loe, &mut l_int2, &mut l_int1, side, tol_conf);
                if ok && !l_int1.is_empty() {
                    shape_data_map::un_bind(to_build, &e);
                }
            }
        }

        // OCCT L3920-3968: the second edge pass (the Inter3D forms).
        for e in explorer(&oriented(&w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
            if shape_data_map::is_bound(const_shapes, &e) {
                shape_data_map::un_bind(to_build, &e);
            }
            if shape_data_map::is_bound(to_build, &e) {
                // OCCT L3929: EnLargeFace(ToBuild(E), StopFace, false).
                en_large_face(&shape_data_map::find(to_build, &e), &mut stop_face, false, false, true, true, true, 1, -1.0, -1.0, -1.0, -1.0);
                // OCCT L3930-3931: Inter3D(EF, StopFace, ..., E, NullFace, NullFace).
                let null_face = Shape::null();
                inter3d(&ef, &stop_face, &mut l_int1, &mut l_int2, side, &e, &null_face, &null_face);
                // No intersection, it may happen for example for a chosen
                // (non-offsetted) planar face and its neighbour offsetted
                // cylindrical face, if the offset is directed so that the
                // radius of the cylinder becomes smaller.
                if l_int1.is_empty() {
                    continue;
                }
                if l_int1.len() > 1 {
                    // l intersection est en plusieurs edges (franchissement
                    // de couture)
                    select_edge(f, &ef, &e, &mut l_int1);
                }
                let mut ne = l_int1.first().cloned().unwrap_or_else(Shape::null);
                // OCCT L3945-3946: TE->Tolerance(TE->Tolerance() * 10.).
                let ne_tol = brep_tool_tolerance(&ne);
                b_set_edge_tolerance(&mut ne, ne_tol * 10.0); //????
                if ne.orientation == e.orientation {
                    shape_data_map::bind(&mut build, &e, oriented(&ne, Orientation::Forward));
                } else {
                    shape_data_map::bind(&mut build, &e, oriented(&ne, Orientation::Reversed));
                }
                // OCCT L3955-3960.
                let e_on_v1 = shape_data_map::find(&n_on_v1, &e);
                if !shape_data_map::is_bound(to_build, &e_on_v1)
                    && !shape_data_map::is_bound(const_shapes, &e_on_v1)
                    && !shape_data_map::is_bound(&build, &e_on_v1)
                {
                    extent_edge_tool(f, &ef, &e_on_v1, &mut ne);
                    shape_data_map::bind(&mut build, &e_on_v1, oriented(&ne, Orientation::Forward));
                }
                // OCCT L3961-3966.
                let e_on_v2 = shape_data_map::find(&n_on_v2, &e);
                if !shape_data_map::is_bound(to_build, &e_on_v2)
                    && !shape_data_map::is_bound(const_shapes, &e_on_v2)
                    && !shape_data_map::is_bound(&build, &e_on_v2)
                {
                    extent_edge_tool(f, &ef, &e_on_v2, &mut ne);
                    shape_data_map::bind(&mut build, &e_on_v2, oriented(&ne, Orientation::Forward));
                }
            }
        }

        //------------------------------------------------
        // Construction Vertex.
        //------------------------------------------------
        // OCCT L3973-4137.
        let mut lv: Vec<Shape> = Vec::new();

        for e in explorer(&oriented(&w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
            let (v1, v2) = top_exp_vertices(&e);
            let (fp, lp) = brep_tool_range(&e);
            let _ = (fp, lp);
            let mut v = Shape::null();
            if shape_data_map::is_bound(&build, &e) {
                // OCCT L3986-4037: the NEOnV1 branch.
                let ne_on_v1 = shape_data_map::find(&n_on_v1, &e);
                if shape_data_map::is_bound(&build, &ne_on_v1)
                    && (shape_data_map::is_bound(to_build, &e)
                        || shape_data_map::is_bound(to_build, &ne_on_v1))
                {
                    if e.is_same(&ne_on_v1) {
                        // OCCT L3989-3992: V = FirstVertex(Build(E)).
                        let build_e = shape_data_map::find(&build, &e);
                        v = top_exp_vertices(&build_e).0;
                    } else {
                        //---------------
                        // intersection.
                        //---------------
                        if !shape_data_map::is_bound(&build, &v1) {
                            // OCCT L4000-4004.
                            inter2d(&ef, &shape_data_map::find(&build, &e), &shape_data_map::find(&build, &ne_on_v1), &mut lv, rcad_kernel::precision::CONFUSION);

                            if !lv.is_empty() {
                                // OCCT L4008-4015.
                                let build_e = shape_data_map::find(&build, &e);
                                if build_e.orientation == Orientation::Forward {
                                    v = lv.first().cloned().unwrap_or_else(Shape::null);
                                } else {
                                    v = lv.last().cloned().unwrap_or_else(Shape::null);
                                }
                            } else {
                                // OCCT L4019: return.
                                return;
                            }
                        } else {
                            // OCCT L4022-4035.
                            v = shape_data_map::find(&build, &v1);
                            if shape_data_map::value(&mve, &v1).len() > 2 {
                                v.orientation = Orientation::Forward;
                                let build_e = shape_data_map::find(&build, &e);
                                if build_e.orientation == Orientation::Reversed {
                                    v.orientation = Orientation::Reversed;
                                }

                                project_vertex_on_edge(&mut v, &shape_data_map::find(&build, &e), tol_conf);
                            }
                        }
                    }
                } else {
                    //------------
                    // projection
                    //------------
                    // OCCT L4039-4056.
                    v = v1.clone();
                    if shape_data_map::is_bound(const_shapes, &v1) {
                        v = shape_data_map::find(const_shapes, &v1);
                    }
                    v.orientation = Orientation::Forward;
                    let build_e = shape_data_map::find(&build, &e);
                    if build_e.orientation == Orientation::Reversed {
                        v.orientation = Orientation::Reversed;
                    }
                    if !try_parameter(&e, &mut v, &shape_data_map::find(&build, &e), tol_conf) {
                        project_vertex_on_edge(&mut v, &shape_data_map::find(&build, &e), tol_conf);
                    }
                }

                // OCCT L4059-4060.
                shape_data_map::bind(const_shapes, &v1, v.clone());
                shape_data_map::bind(&mut build, &v1, v.clone());
                // OCCT L4061-4113: the NEOnV2 branch.
                let ne_on_v2 = shape_data_map::find(&n_on_v2, &e);
                if shape_data_map::is_bound(&build, &ne_on_v2)
                    && (shape_data_map::is_bound(to_build, &e)
                        || shape_data_map::is_bound(to_build, &ne_on_v2))
                {
                    if e.is_same(&ne_on_v2) {
                        // OCCT L4064-4066: V = LastVertex(Build(E)).
                        let build_e = shape_data_map::find(&build, &e);
                        v = top_exp_vertices(&build_e).1;
                    } else {
                        //--------------
                        // intersection.
                        //---------------
                        if !shape_data_map::is_bound(&build, &v2) {
                            // OCCT L4076-4080.
                            inter2d(&ef, &shape_data_map::find(&build, &e), &shape_data_map::find(&build, &ne_on_v2), &mut lv, rcad_kernel::precision::CONFUSION);

                            if !lv.is_empty() {
                                // OCCT L4084-4091.
                                let build_e = shape_data_map::find(&build, &e);
                                if build_e.orientation == Orientation::Forward {
                                    v = lv.last().cloned().unwrap_or_else(Shape::null);
                                } else {
                                    v = lv.first().cloned().unwrap_or_else(Shape::null);
                                }
                            } else {
                                // OCCT L4095: return.
                                return;
                            }
                        } else {
                            // OCCT L4098-4111.
                            v = shape_data_map::find(&build, &v2);
                            if shape_data_map::value(&mve, &v2).len() > 2 {
                                v.orientation = Orientation::Reversed;
                                let build_e = shape_data_map::find(&build, &e);
                                if build_e.orientation == Orientation::Reversed {
                                    v.orientation = Orientation::Forward;
                                }

                                project_vertex_on_edge(&mut v, &shape_data_map::find(&build, &e), tol_conf);
                            }
                        }
                    }
                } else {
                    //------------
                    // projection
                    //------------
                    // OCCT L4115-4132.
                    v = v2.clone();
                    if shape_data_map::is_bound(const_shapes, &v2) {
                        v = shape_data_map::find(const_shapes, &v2);
                    }
                    v.orientation = Orientation::Reversed;
                    let build_e = shape_data_map::find(&build, &e);
                    if build_e.orientation == Orientation::Reversed {
                        v.orientation = Orientation::Forward;
                    }
                    if !try_parameter(&e, &mut v, &shape_data_map::find(&build, &e), tol_conf) {
                        project_vertex_on_edge(&mut v, &shape_data_map::find(&build, &e), tol_conf);
                    }
                }
                // OCCT L4134-4135.
                shape_data_map::bind(const_shapes, &v2, v.clone());
                shape_data_map::bind(&mut build, &v2, v.clone());
            }
        }

        //-----------------
        // Reconstruction.
        //-----------------
        // OCCT L4139-4329.
        let mut nw = bat::builder_make_wire();

        for e in explorer(&oriented(&w, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
            let (v1, v2) = top_exp_vertices(&e);
            let mut ne: Shape;
            if shape_data_map::is_bound(&build, &e) {
                // OCCT L4157-4184: the Build branch.
                ne = shape_data_map::find(&build, &e);
                let (fp, lp) = brep_tool_range(&ne);
                let or = ne.orientation;
                //-----------------------------------------------------
                // Copy pour virer les vertex deja sur la nouvelle edge.
                //-----------------------------------------------------
                let nv1 = shape_data_map::find(const_shapes, &v1);
                let nv2 = shape_data_map::find(const_shapes, &v2);

                // OCCT L4171-4176: the INTERNAL-parameter reads.
                let ne_fwd = oriented(&ne, Orientation::Forward);
                let u1 = bat::brep_tool_parameter(&nv1, &ne_fwd);
                let u2 = bat::brep_tool_parameter(&nv2, &ne_fwd);

                // OCCT L4183-4185: NE = NE.EmptyCopied(); FORWARD.
                ne = empty_copied(&ne);
                ne.orientation = Orientation::Forward;
                if nv1.is_same(&nv2) {
                    //--------------
                    // edge ferme.
                    //--------------
                    // OCCT L4188-4246.
                    let (mut u1, mut u2) = (u1, u2);
                    if or == Orientation::Forward {
                        u1 = fp;
                        u2 = lp;
                    } else {
                        u1 = lp;
                        u2 = fp;
                    }
                    if or == Orientation::Forward {
                        if u1 > u2 {
                            if (u1 - lp).abs() < rcad_kernel::precision::CONFUSION {
                                u1 = fp;
                            }
                            if (u2 - fp).abs() < rcad_kernel::precision::CONFUSION {
                                u2 = lp;
                            }
                        }
                        let mut nv1_fwd = nv1.clone();
                        nv1_fwd.orientation = Orientation::Forward;
                        bat::builder_add_edge_vertex(&mut ne, &nv1_fwd);
                        let mut nv2_rev = nv2.clone();
                        nv2_rev.orientation = Orientation::Reversed;
                        bat::builder_add_edge_vertex(&mut ne, &nv2_rev);
                        bat::builder_range_edge(&mut ne, u1, u2);
                        shape_data_map::bind(const_shapes, &e, ne.clone());
                        ne.orientation = e.orientation;
                    } else {
                        if u2 > u1 {
                            if (u2 - lp).abs() < rcad_kernel::precision::CONFUSION {
                                u2 = fp;
                            }
                            if (u1 - fp).abs() < rcad_kernel::precision::CONFUSION {
                                u1 = lp;
                            }
                        }
                        let mut nv2_fwd = nv2.clone();
                        nv2_fwd.orientation = Orientation::Forward;
                        bat::builder_add_edge_vertex(&mut ne, &nv2_fwd);
                        let mut nv1_rev = nv1.clone();
                        nv1_rev.orientation = Orientation::Reversed;
                        bat::builder_add_edge_vertex(&mut ne, &nv1_rev);
                        bat::builder_range_edge(&mut ne, u2, u1);
                        shape_data_map::bind(const_shapes, &e, oriented(&ne, Orientation::Reversed));
                        ne.orientation = bat::top_abs_reverse(e.orientation);
                    }
                } else {
                    //-------------------
                    // edge is not ferme.
                    //-------------------
                    // OCCT L4248-4305.
                    if or == Orientation::Forward {
                        if u1 > u2 {
                            let mut nv2_fwd = nv2.clone();
                            nv2_fwd.orientation = Orientation::Forward;
                            bat::builder_add_edge_vertex(&mut ne, &nv2_fwd);
                            let mut nv1_rev = nv1.clone();
                            nv1_rev.orientation = Orientation::Reversed;
                            bat::builder_add_edge_vertex(&mut ne, &nv1_rev);
                            bat::builder_range_edge(&mut ne, u2, u1);
                        } else {
                            let mut nv1_fwd = nv1.clone();
                            nv1_fwd.orientation = Orientation::Forward;
                            bat::builder_add_edge_vertex(&mut ne, &nv1_fwd);
                            let mut nv2_rev = nv2.clone();
                            nv2_rev.orientation = Orientation::Reversed;
                            bat::builder_add_edge_vertex(&mut ne, &nv2_rev);
                            bat::builder_range_edge(&mut ne, u1, u2);
                        }
                        shape_data_map::bind(const_shapes, &e, ne.clone());
                        ne.orientation = e.orientation;
                    } else if u2 > u1 {
                        let mut nv1_fwd = nv1.clone();
                        nv1_fwd.orientation = Orientation::Forward;
                        bat::builder_add_edge_vertex(&mut ne, &nv1_fwd);
                        let mut nv2_rev = nv2.clone();
                        nv2_rev.orientation = Orientation::Reversed;
                        bat::builder_add_edge_vertex(&mut ne, &nv2_rev);
                        bat::builder_range_edge(&mut ne, u1, u2);
                        shape_data_map::bind(const_shapes, &e, ne.clone());
                        ne.orientation = e.orientation;
                    } else {
                        let mut nv2_fwd = nv2.clone();
                        nv2_fwd.orientation = Orientation::Forward;
                        bat::builder_add_edge_vertex(&mut ne, &nv2_fwd);
                        let mut nv1_rev = nv1.clone();
                        nv1_rev.orientation = Orientation::Reversed;
                        bat::builder_add_edge_vertex(&mut ne, &nv1_rev);
                        bat::builder_range_edge(&mut ne, u2, u1);
                        shape_data_map::bind(const_shapes, &e, oriented(&ne, Orientation::Reversed));
                        ne.orientation = bat::top_abs_reverse(e.orientation);
                    }
                }
                // OCCT L4306: Build.UnBind(E).
                shape_data_map::un_bind(&mut build, &e);
            } else if shape_data_map::is_bound(const_shapes, &e) {
                // !Build.IsBound(E)
                // OCCT L4308-4321: the ConstShapes branch.
                ne = shape_data_map::find(const_shapes, &e);
                build_pcurves(&ne, nf);
                let or = ne.orientation;
                if or == Orientation::Reversed {
                    ne.orientation = bat::top_abs_reverse(e.orientation);
                } else {
                    ne.orientation = e.orientation;
                }
            } else {
                // OCCT L4322-4326.
                ne = e.clone();
                shape_data_map::bind(const_shapes, &e, oriented(&ne, Orientation::Forward));
            }
            // OCCT L4327: B.Add(NW, NE).
            bat::builder_add_wire_edge(&mut nw, &ne);
        }
        // OCCT L4329: B.Add(NF, NW.Oriented(W.Orientation())).
        bat::builder_add_face_wire(nf, &oriented(&nw, w.orientation));
    }
    // OCCT L4331-4332.
    nf.orientation = f.orientation;
    // BRepTools::Update(NF) — GAP no-op (architecture difference #31).
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_Tool::Deboucle3D (hxx L191-193; cxx L4337-4418).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Tool::Deboucle3D(S, Boundary) (cxx L4337-4418) — remove
/// the non valid part of an offset shape.
pub fn deboucle3d(s: &Shape, boundary: &OcctShapeSet) -> Shape {
    let mut ss = Shape::null();
    match s.shape_type() {
        ShapeType::Shell => {
            // if the shell contains free borders that do not belong to the
            // free borders of caps (Boundary) it is removed.
            // OCCT L4347-4351: MapShapesAndAncestors(S, EDGE, FACE, Map).
            let mut map: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
            crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(s, ShapeType::Edge, ShapeType::Face, &mut map);

            // OCCT L4353-4373: the JeGarde walk.
            let mut je_garde = true;
            let mut i = 1usize;
            while i <= shape_indexed_data_map::extent(&map) && je_garde {
                let a_lf = shape_indexed_data_map::value_1(&map, i);
                if a_lf.len() < 2 {
                    let an_edge = shape_indexed_data_map::find_key_1(&map, i);
                    if an_edge.orientation == Orientation::Internal {
                        let a_face = &a_lf[0];
                        if a_face.orientation != Orientation::Internal {
                            i += 1;
                            continue;
                        }
                    }
                    if !set_contains(boundary, an_edge) && !brep_tool_degenerated(an_edge) {
                        je_garde = false;
                    }
                }
                i += 1;
            }
            if je_garde {
                ss = s.clone();
            }
        }
        ShapeType::Compound | ShapeType::Solid => {
            // iterate on sub-shapes and add non-empty.
            // OCCT L4384-4409.
            let mut nb_sub = 0usize;
            if s.shape_type() == ShapeType::Compound {
                ss = bat::builder_make_compound();
            } else {
                ss = b_make_solid();
            }
            for cur_s in bat::sub_shapes(s) {
                let sub_shape = deboucle3d(&cur_s, boundary);
                if !sub_shape.is_null() {
                    if s.shape_type() == ShapeType::Compound {
                        bat::builder_add_compound_shape(&mut ss, &sub_shape);
                    } else {
                        b_add_solid_shell(&mut ss, &sub_shape);
                    }
                    nb_sub += 1;
                }
            }
            if nb_sub == 0 {
                ss = Shape::null();
            }
        }
        _ => {}
    }

    // OCCT L4417: return SS.
    ss
}
