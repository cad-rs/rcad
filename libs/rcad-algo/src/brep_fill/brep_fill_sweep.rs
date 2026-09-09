//! OCCT BRepFill_Sweep (TKBool/BRepFill) — 1:1 translation of
//! BRepFill_Sweep.hxx (L43-189) + BRepFill_Sweep.cxx, part A (this file):
//! the class shell, the simple setters/accessors, SetBounds, BuildWire,
//! MergeVertex / UpdateVertex / EvalExtrapol and the self-contained file
//! statics (Translate / Box / BuildVertex / NullEdge / HasPCurves /
//! ReverseEdgeInFirstOrLastWire / ReverseModifiedEdges).
//!
//! Part B (re-dispatch, listed in the delivery report) — the remaining
//! translation units of BRepFill_Sweep.cxx: NumberOfPoles, Couture,
//! CheckSameParameter / CheckSameParameterExact / SameParameter /
//! CorrectSameParameter, Oriente, UpdateEdgeOnPlane, BuildFace, BuildEdge
//! (both forms), Filling, Substitute, KeepEdge, UpdateEdge, IsDegen,
//! CorrectApproxParameters, BuildShell, Build, PerformCorner,
//! RebuildTopOrBottomEdge.
//!
//! Architecture differences:
//! - `handle(BRepFill_SectionLaw) mySec` /
//!   `handle(BRepFill_LocationLaw) myLoc` map to
//!   `Rc<RefCell<dyn ...Ops>>` (the trait slots of the brep_fill law files).
//! - `NCollection_HArray2<TopoDS_Shape>` maps to
//!   [`crate::brep_fill::brep_fill_pipe_shell_b::ShapeHArray2`]
//!   (`Vec<Vec<Shape>>`, row-major, 1-based OCCT indexing -> `[r-1][c-1]`).
//! - `NCollection_Map / NCollection_DataMap` keyed by
//!   TopTools_ShapeMapHasher map to HashSet/HashMap keyed by `ShapeKey`.
//! - `BRep_Builder` maps to `rcad_kernel::topo::topods::BRepBuilder`; the
//!   owning pool is the leading `brep` argument.
//! - `BRep_Tool::Tolerance / Pnt` map to the stored TVertexData fields.
//! - `GeomAdaptor/BRepAdaptor` surface bounds map to
//!   `surface_adaptor_basis_and_bounds`; the UIso/VIso dispatch is the local
//!   ElSLib / BSpline re-host ([`surface_uiso`] / [`surface_viso`]).

use std::cell::RefCell;
use std::collections::HashSet;
use std::rc::Rc;

use glam::DVec3;

use rcad_kernel::geom::{CurveEval, BSplineSurface, Circle3, Curve3, Ellipse3, Line3, Surface3, TrimmedCurve3};
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::topo::topods::{surface_adaptor_basis_and_bounds, BRep, BRepBuilder, Shape};

use crate::brep_fill::brep_fill_location_law::BRepFillLocationLawOps;
use crate::brep_fill::brep_fill_pipe_shell_b::{
    ShapeHArray2, ShapeSet,
};
use crate::brep_fill::generator::{shape_reversed, ShapeKey};
use crate::geomalgo::geomfill::section_law::SectionLaw;
use crate::geomalgo::geomfill::sweep::GeomFillApproxStyle;

/// OCCT BRepFill_TransitionStyle (BRepFill_TransitionStyle.hxx).
pub use crate::brep_fill::brep_fill_pipe_shell_b::BRepFillTransitionStyle;

// ---------------------------------------------------------------------------
// UIso / VIso re-hosts (OCCT Geom_Surface::UIso / VIso dispatch)
// ---------------------------------------------------------------------------

/// OCCT ElSLib + Geom_BSplineSurface::UIso over the rcad Surface3 — the
/// u-varying iso curve at `u`.  The Offset / remaining variants keep the
/// OCCT failure path (kernel GAP).
pub(super) fn surface_uiso(surf: &Surface3, u: f64) -> Curve3 {    match surf {
        // ElSLib::PlaneUIso (ElSLib.cxx): line through P(u, 0) along the V
        // direction.
        Surface3::Plane(pl) => Curve3::Line(Line3 {
            origin: pl.origin + pl.u_dir * u,
            direction: pl.normal.cross(pl.u_dir).normalize_or_zero(),
        }),
        // ElSLib::CylinderUIso: circle (Loc, axis, R).
        Surface3::Cylinder(cy) => Curve3::Circle(Circle3 {
            center: cy.origin,
            normal: cy.axis,
            x_dir: cy.ref_dir,
            y_dir: cy.axis.cross(cy.ref_dir).normalize_or_zero(),
            radius: cy.radius,
        }),
        // ElSLib::SphereUIso: the meridian great circle of longitude u.
        Surface3::Sphere(sp) => {
            let x = sp.ref_dir;
            let y = sp.axis.cross(sp.ref_dir).normalize_or_zero();
            let normal = (x * u.cos() + y * u.sin()).normalize_or_zero();
            Curve3::Circle(Circle3 {
                center: sp.center,
                normal,
                x_dir: sp.axis,
                y_dir: normal.cross(sp.axis).normalize_or_zero(),
                radius: sp.radius,
            })
        }
        // Geom_ConicalSurface::UIso: the ruling line at longitude u (the
        // fixed-u iso varies v along the ruling; P(u0, v) = apex_true
        // + v*(cos(a) * axis + sin(a) * radial(u0))).
        Surface3::Cone(co) => {
            let a = co.half_angle_rad;
            let apex = co.apex - (co.radius / a.tan()) * co.axis;
            let x = co.ref_dir;
            let y = co.axis.cross(co.ref_dir).normalize_or_zero();
            let radial = (x * u.cos() + y * u.sin()).normalize_or_zero();
            Curve3::Line(Line3 {
                origin: apex,
                direction: (co.axis * a.cos() + radial * a.sin()).normalize_or_zero(),
            })
        }
        // ElSLib::TorusUIso (the geomfill/sweep.rs re-host): the minor
        // circle at longitude u.
        Surface3::Torus(to) => {
            let x = to.ref_dir;
            let y = to.axis.cross(to.ref_dir).normalize_or_zero();
            let radial = (x * u.cos() + y * u.sin()).normalize_or_zero();
            Curve3::Circle(Circle3 {
                center: to.center + radial * to.major_radius,
                normal: radial,
                x_dir: to.axis,
                y_dir: radial.cross(to.axis).normalize_or_zero(),
                radius: to.minor_radius,
            })
        }
        // Geom_BSplineSurface::UIso: the poles of the iso are the
        // homogeneous De Boor evaluation of the u basis per V column.
        Surface3::BSpline(bs) => Curve3::BSpline(bspline_surface_uiso(bs, u)),
        // OCCT Geom_SurfaceOfLinearExtrusion::UIso
        // (Geom_SurfaceOfLinearExtrusion.cxx L275-281): the ruling line
        // (basisCurve->Value(U), direction).
        Surface3::LinearExtrusion(le) => Curve3::Line(Line3 {
            origin: le.profile.point_at(u),
            direction: le.direction,
        }),
        // OCCT Geom_SurfaceOfRevolution::UIso (Geom_SurfaceOfRevolution.cxx
        // L372-381): the basis curve rotated by U about the axis
        // (C->Rotate(Ax1(loc, direction), U)).
        Surface3::Revolution(rev) => {
            curve3_rotated_about_axis(&rev.profile, rev.axis_origin, rev.axis_dir, u)
        }
        // OCCT Geom_RectangularTrimmedSurface::UIso
        // (Geom_RectangularTrimmedSurface.cxx L444-458): the basis UIso,
        // restricted to the v-trim range when isvtrimmed (the rcad flag maps
        // to "the stored v-trim bounds differ from the basis natural domain").
        Surface3::Trimmed(ts) => {
            let c = surface_uiso(&ts.basis, u);
            if basis_v_bounds_of(&ts.basis) != (ts.trim[2], ts.trim[3]) {
                Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(c),
                    first: ts.trim[2],
                    last: ts.trim[3],
                })
            } else {
                c
            }
        }
        // OCCT Geom_OffsetSurface::UIso (Geom_OffsetSurface.cxx L601-655)
        // approximates the offset iso through AdvApprox_ApproxAFunction and
        // the Geom_OffsetSurface_UIsoEvaluator — the AdvApprox evaluator
        // wiring is a kernel GAP (plan §0.6).
        Surface3::Offset(_) => panic!(
            "GAP: Geom_OffsetSurface::UIso (TKG3d/Geom, AdvApprox evaluator \
             wiring) — BRepFill_Sweep::BuildWire"
        ),
        _ => panic!(
            "GAP: Geom_*Surface::UIso (TKMath/TKG3d kernel re-host) is not \
             translated for this surface type — BRepFill_Sweep::BuildWire"
        ),
    }
}

/// The natural v-bounds of a basis surface (the OCCT `isvtrimmed` check
/// compares the trim bounds against the basis domain).
fn basis_v_bounds_of(surf: &Surface3) -> (f64, f64) {
    use rcad_kernel::geom::SurfaceEval;
    let d = surf.default_domain();
    (d[2], d[3])
}

/// OCCT ElSLib + Geom_BSplineSurface::VIso over the rcad Surface3 — the
/// v-varying iso curve at `v`.  Same GAP note as [`surface_uiso`].
pub(super) fn surface_viso(surf: &Surface3, v: f64) -> Curve3 {
    match surf {
        // ElSLib::PlaneVIso: line through P(0, v) along the U direction.
        Surface3::Plane(pl) => Curve3::Line(Line3 {
            origin: pl.origin + pl.normal.cross(pl.u_dir) * v,
            direction: pl.u_dir,
        }),
        // ElSLib::CylinderVIso: the ruling line at longitude 0, height v.
        Surface3::Cylinder(cy) => Curve3::Line(Line3 {
            origin: cy.origin + cy.axis * v,
            direction: cy.axis,
        }),
        // ElSLib::SphereVIso: the parallel circle at latitude v.
        Surface3::Sphere(sp) => Curve3::Circle(Circle3 {
            center: sp.center + sp.axis * (sp.radius * v.cos()),
            normal: sp.axis,
            x_dir: sp.ref_dir,
            y_dir: sp.axis.cross(sp.ref_dir).normalize_or_zero(),
            radius: sp.radius * v.sin(),
        }),
        // ElSLib::ConeVIso: the parallel circle at parameter v.
        Surface3::Cone(co) => {
            let a = co.half_angle_rad;
            let apex = co.apex - (co.radius / a.tan()) * co.axis;
            let center = apex + v * a.cos() * co.axis;
            Curve3::Circle(Circle3 {
                center,
                normal: co.axis,
                x_dir: co.ref_dir,
                y_dir: co.axis.cross(co.ref_dir).normalize_or_zero(),
                radius: v * a.sin(),
            })
        }
        // ElSLib::TorusVIso: the major circle at parameter v.
        Surface3::Torus(to) => Curve3::Circle(Circle3 {
            center: to.center + to.axis * (to.minor_radius * v.sin()),
            normal: to.axis,
            x_dir: to.ref_dir,
            y_dir: to.axis.cross(to.ref_dir).normalize_or_zero(),
            radius: to.major_radius + to.minor_radius * v.cos(),
        }),
        // Geom_BSplineSurface::VIso.
        Surface3::BSpline(bs) => Curve3::BSpline(bspline_surface_viso_full(bs, v)),
        // OCCT Geom_SurfaceOfLinearExtrusion::VIso
        // (Geom_SurfaceOfLinearExtrusion.cxx L285-291): the basis curve
        // translated by V*direction (Vdir.Multiply(V); basis->Translated).
        Surface3::LinearExtrusion(le) => {
            curve3_translated(&le.profile, le.direction * v)
        }
        // OCCT Geom_SurfaceOfRevolution::VIso (Geom_SurfaceOfRevolution.cxx
        // L383-411): the parallel circle (Loc, axis, Rad) built at the basis
        // point of parameter V.
        Surface3::Revolution(rev) => {
            // Pnt Pc = basisCurve->Value(V);
            let pc = rev.profile.point_at(v);
            // gp_Lin L1(loc, direction); Rad = L1.Distance(Pc).
            let axis = rev.axis_dir;
            let rad = (pc - rev.axis_origin - (pc - rev.axis_origin).dot(axis) * axis).length();
            // Ax2 Rep — the circle frame: center C on the axis, XDir D the
            // radial direction (P = Pc - C normalized when Rad > Resolution).
            let c = rev.axis_origin + (pc - rev.axis_origin).dot(axis) * axis;
            let x_dir = if rad > GP_RESOLUTION {
                let d = pc - c;
                if d.length() > GP_RESOLUTION {
                    d.normalize_or_zero()
                } else {
                    // gp_Ax2(C, direction) default X — any orthogonal
                    // direction (the OCCT default picks a solver-dependent
                    // one; the radial projection is degenerate here).
                    DVec3::Z.cross(axis).normalize_or_zero()
                }
            } else {
                DVec3::Z.cross(axis).normalize_or_zero()
            };
            Curve3::Circle(Circle3 {
                center: c,
                normal: axis,
                x_dir,
                y_dir: axis.cross(x_dir).normalize_or_zero(),
                radius: rad,
            })
        }
        // OCCT Geom_RectangularTrimmedSurface::VIso
        // (Geom_RectangularTrimmedSurface.cxx L463-477): the basis VIso,
        // restricted to the u-trim range when isutrimmed.
        Surface3::Trimmed(ts) => {
            let c = surface_viso(&ts.basis, v);
            if basis_v_bounds_of(&ts.basis) != (ts.trim[0], ts.trim[1]) {
                Curve3::Trimmed(TrimmedCurve3 {
                    curve: Box::new(c),
                    first: ts.trim[0],
                    last: ts.trim[1],
                })
            } else {
                c
            }
        }
        // OCCT Geom_OffsetSurface::VIso (Geom_OffsetSurface.cxx L657-706) —
        // the same AdvApprox GAP as UIso.
        Surface3::Offset(_) => panic!(
            "GAP: Geom_OffsetSurface::VIso (TKG3d/Geom, AdvApprox evaluator \
             wiring) — BRepFill_Sweep::BuildWire"
        ),
        _ => panic!(
            "GAP: Geom_*Surface::VIso (TKMath/TKG3d kernel re-host) is not \
             translated for this surface type — BRepFill_Sweep::BuildWire"
        ),
    }
}

/// OCCT gp_GTrsf-based `Geom_Curve::Translated(T)` over the rcad Curve3 —
/// supported for the analytic / poles-based variants (the OCCT operation
/// translates the geometry representation in place); the remaining curve
/// kinds keep the OCCT failure path.
fn curve3_translated(c: &Curve3, t: DVec3) -> Curve3 {
    match c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: l.origin + t,
            direction: l.direction,
        }),
        Curve3::Circle(ci) => Curve3::Circle(Circle3 {
            center: ci.center + t,
            ..ci.clone()
        }),
        Curve3::Ellipse(el) => Curve3::Ellipse(Ellipse3 {
            center: el.center + t,
            ..el.clone()
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(rcad_kernel::geom::BSplineCurve3 {
            control_points: bs.control_points.iter().map(|p| p + t).collect(),
            ..bs.clone()
        }),
        Curve3::Bezier(bz) => Curve3::Bezier(rcad_kernel::geom::BezierCurve3 {
            control_points: bz.control_points.iter().map(|p| p + t).collect(),
            weights: bz.weights.clone(),
        }),
        Curve3::Trimmed(tr) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(curve3_translated(&tr.curve, t)),
            first: tr.first,
            last: tr.last,
        }),
        _ => panic!(
            "GAP: Geom_Curve::Translated (TKG3d/Geom kernel re-host) for this \
             curve type — Geom_SurfaceOfLinearExtrusion::VIso"
        ),
    }
}

/// OCCT `Geom_Curve::Rotate(Ax1, Angle)` over the rcad Curve3 — the
/// rotation about the (origin, direction) axis applied to the geometry
/// representation; supported for the analytic / poles-based variants.
fn curve3_rotated_about_axis(c: &Curve3, origin: DVec3, direction: DVec3, angle: f64) -> Curve3 {
    use glam::DAffine3;
    let axis = DAffine3::from_rotation_translation(
        glam::DQuat::from_axis_angle(direction, angle),
        DVec3::ZERO,
    );
    let rot = |p: DVec3| origin + axis.transform_point3(p - origin);    match c {
        Curve3::Line(l) => Curve3::Line(Line3 {
            origin: rot(l.origin),
            direction: axis.transform_vector3(l.direction),
        }),
        Curve3::Circle(ci) => Curve3::Circle(Circle3 {
            center: rot(ci.center),
            normal: axis.transform_vector3(ci.normal),
            x_dir: axis.transform_vector3(ci.x_dir),
            y_dir: axis.transform_vector3(ci.y_dir),
            radius: ci.radius,
        }),
        Curve3::Ellipse(el) => Curve3::Ellipse(Ellipse3 {
            center: rot(el.center),
            normal: axis.transform_vector3(el.normal),
            major_dir: axis.transform_vector3(el.major_dir),
            major_radius: el.major_radius,
            minor_radius: el.minor_radius,
        }),
        Curve3::BSpline(bs) => Curve3::BSpline(rcad_kernel::geom::BSplineCurve3 {
            control_points: bs.control_points.iter().map(|p| rot(*p)).collect(),
            ..bs.clone()
        }),
        Curve3::Bezier(bz) => Curve3::Bezier(rcad_kernel::geom::BezierCurve3 {
            control_points: bz.control_points.iter().map(|p| rot(*p)).collect(),
            weights: bz.weights.clone(),
        }),
        Curve3::Trimmed(tr) => Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(curve3_rotated_about_axis(&tr.curve, origin, direction, angle)),
            first: tr.first,
            last: tr.last,
        }),
        _ => panic!(
            "GAP: Geom_Curve::Rotate (TKG3d/Geom kernel re-host) for this \
             curve type — Geom_SurfaceOfRevolution::UIso"
        ),
    }
}

/// OCCT Geom_BSplineSurface::UIso — the iso poles are the homogeneous De
/// Boor evaluation of the u basis per V column (the nsections.rs
/// bspline_surface_viso precedent, transposed).
fn bspline_surface_uiso(surf: &BSplineSurface, u: f64) -> rcad_kernel::geom::BSplineCurve3 {
    use rcad_kernel::math::bspl::de_boor_homo;
    let rational = surf.weights.iter().any(|row| row.iter().any(|&w| w != 1.0));
    let nb_u = surf.control_points.len();
    let nb_v = surf.control_points.first().map(|r| r.len()).unwrap_or(0);
    let mut poles = Vec::with_capacity(nb_v);
    let mut weights = Vec::with_capacity(nb_v);
    for jj in 0..nb_v {
        let row: Vec<DVec3> = (0..nb_u).map(|ii| surf.control_points[ii][jj]).collect();
        let row_w: Vec<f64> = (0..nb_u)
            .map(|ii| if rational { surf.weights[ii][jj] } else { 1.0 })
            .collect();
        let h = de_boor_homo(surf.degree_u, &surf.knots_u, &row, &row_w, u);
        weights.push(h[3]);
        poles.push(DVec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]));
    }
    rcad_kernel::geom::BSplineCurve3 {
        degree: surf.degree_v,
        knots: surf.knots_v.clone(),
        control_points: poles,
        weights,
        is_periodic: false,
    }
}

/// OCCT Geom_BSplineSurface::VIso (same form as the BRepFill_NSections
/// iso re-host).
fn bspline_surface_viso_full(surf: &BSplineSurface, v: f64) -> rcad_kernel::geom::BSplineCurve3 {
    use rcad_kernel::math::bspl::de_boor_homo;
    let rational = surf.weights.iter().any(|row| row.iter().any(|&w| w != 1.0));
    let nb_u = surf.control_points.len();
    let nb_v = surf.control_points.first().map(|r| r.len()).unwrap_or(0);
    let mut poles = Vec::with_capacity(nb_u);
    let mut weights = Vec::with_capacity(nb_u);
    for ii in 0..nb_u {
        let column: Vec<DVec3> = (0..nb_v).map(|jj| surf.control_points[ii][jj]).collect();
        let column_w: Vec<f64> = (0..nb_v)
            .map(|jj| if rational { surf.weights[ii][jj] } else { 1.0 })
            .collect();
        let h = de_boor_homo(surf.degree_v, &surf.knots_v, &column, &column_w, v);
        weights.push(h[3]);
        poles.push(DVec3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]));
    }
    rcad_kernel::geom::BSplineCurve3 {
        degree: surf.degree_u,
        knots: surf.knots_u.clone(),
        control_points: poles,
        weights,
        is_periodic: false,
    }
}

// ---------------------------------------------------------------------------
// File statics (part A)
// ---------------------------------------------------------------------------

/// OCCT static Translate (cxx L178-189) — Copy a column from one table to
/// another (1-based OCCT indexing -> `[i-1]`).
pub(super) fn translate(
    array_in: &ShapeHArray2,
    in_col: i32,
    array_out: &mut ShapeHArray2,
    out_col: i32,
) {
    let nb = array_out[0].len();
    for ii in 1..=nb {
        array_out[ii - 1][out_col as usize - 1] = array_in[ii - 1][in_col as usize - 1].clone();
    }
}

/// OCCT static Box (cxx L193-206) — the section poles bounding box (the
/// Bnd_Box SetVoid / Add / Get triple reduced to the accumulated min/max,
/// the pure-math Bnd_Box re-host).
pub(super) fn box_of_section(sec: &Rc<RefCell<dyn SectionLaw>>, u: f64) -> ([f64; 3], [f64; 3]) {
    let mut nb_poles = 0usize;
    let mut bid = 0usize;
    {
        let mut bid2 = 0usize;
        sec.borrow().section_shape(&mut nb_poles, &mut bid, &mut bid2);
    }
    let mut poles = vec![DVec3::ZERO; nb_poles];
    let mut w = vec![0.0f64; nb_poles];
    sec.borrow().d0(u, &mut poles, &mut w);
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    for p in &poles {
        for k in 0..3 {
            min[k] = min[k].min(p[k]);
            max[k] = max[k].max(p[k]);
        }
    }
    (min, max)
}

/// OCCT static BuildVertex (cxx L1440-1458).
pub(super) fn build_vertex(
    brep: &mut BRep,
    iso: &Curve3,
    isfirst: bool,
    first: f64,
    last: f64,
    vertex: &mut Shape,
) {
    let mut b = BRepBuilder::new();
    let val = if isfirst { first } else { last };
    // B.MakeVertex(TopoDS::Vertex(Vertex), Iso->Value(val), Confusion());
    *vertex = b.add_vertex(brep, iso.point_at(val), rcad_kernel::core::precision::CONFUSION);
}

/// OCCT static NullEdge (cxx L1462-1472).
pub(super) fn null_edge(brep: &mut BRep, vertex: &mut Shape) -> Shape {
    use rcad_kernel::topo::topods::Orientation;
    let mut b = BRepBuilder::new();
    let mut e = b.add_edge(brep, None, vertex.clone(), shape_reversed(vertex), [0.0, 0.0]);
    // Vertex.Orientation(TopAbs_FORWARD); B.Add(E, Vertex); B.Add(E,
    // Vertex.Reversed()); B.Degenerated(E, true);
    vertex.orientation = Orientation::Forward;
    b.set_edge_degenerated(brep, e.clone(), true);
    e
}

/// OCCT static HasPCurves (cxx L156-172) — the edge carries a
/// curve-on-surface representation.
pub(super) fn has_pcurves(brep: &BRep, e: &Shape) -> bool {
    !brep.edge(e.clone()).pcurves.is_empty()
}

/// OCCT static ReverseEdgeInFirstOrLastWire (cxx L1877-1900).
pub(super) fn reverse_edge_in_first_or_last_wire(brep: &mut BRep, the_wire: &mut Shape, the_edge: &Shape) {
    // TopoDS_Iterator itw(theWire) — the wire-ordered children.
    let mut edge_to_reverse: Option<Shape> = None;
    for an_edge in crate::brep_fill::compatible_wires::wire_edges(brep, the_wire) {
        if an_edge.is_same(the_edge) {
            edge_to_reverse = Some(an_edge);
            break;
        }
    }

    if let Some(mut edge) = edge_to_reverse {
        let mut bb = BRepBuilder::new();
        edge.orientation = match edge.orientation {
            rcad_kernel::topo::topods::Orientation::Forward => {
                rcad_kernel::topo::topods::Orientation::Reversed
            }
            rcad_kernel::topo::topods::Orientation::Reversed => {
                rcad_kernel::topo::topods::Orientation::Forward
            }
            o => o,
        };
        bb.remove_from_wire(brep, the_wire.clone(), edge.clone());
        bb.add_to_wire(brep, the_wire.clone(), edge);
    }
}

/// OCCT static ReverseModifiedEdges (cxx L1902-1936).
pub(super) fn reverse_modified_edges(brep: &mut BRep, the_wire: &mut Shape, the_emap: &ShapeSet) {
    if the_emap.is_empty() {
        return;
    }

    let mut bb = BRepBuilder::new();
    let ledges = crate::brep_fill::compatible_wires::wire_edges(brep, the_wire);
    for e in &ledges {
        bb.remove_from_wire(brep, the_wire.clone(), e.clone());
    }

    for e in ledges {
        let mut an_edge = e;
        if the_emap.contains(&ShapeKey(an_edge.ptr_id())) {
            an_edge = shape_reversed(&an_edge);
        }
        bb.add_to_wire(brep, the_wire.clone(), an_edge);
    }
}

// ---------------------------------------------------------------------------
// BRepFill_Sweep
// ---------------------------------------------------------------------------

/// OCCT BRepFill_Sweep (hxx L43-189): Topological Sweep Algorithm.
pub struct BRepFillSweep {
    /// OCCT bool isDone.
    pub is_done: bool,
    /// OCCT bool KPart.
    pub kpart: bool,
    /// OCCT double myTol3d.
    pub my_tol3d: f64,
    /// OCCT double myBoundTol.
    pub my_bound_tol: f64,
    /// OCCT double myTol2d.
    pub my_tol2d: f64,
    /// OCCT double myTolAngular.
    pub my_tol_angular: f64,
    /// OCCT double myAngMin.
    pub my_ang_min: f64,
    /// OCCT double myAngMax.
    pub my_ang_max: f64,
    /// OCCT GeomFill_ApproxStyle myApproxStyle.
    pub my_approx_style: GeomFillApproxStyle,
    /// OCCT GeomAbs_Shape myContinuity.
    pub my_continuity: GeomAbsShape,
    /// OCCT int myDegmax.
    pub my_degmax: i32,
    /// OCCT int mySegmax.
    pub my_segmax: i32,
    /// OCCT bool myForceApproxC1.
    pub my_force_approx_c1: bool,
    /// OCCT TopoDS_Shape myShape.
    pub my_shape: Shape,
    /// OCCT handle(BRepFill_LocationLaw) myLoc.
    pub my_loc: Rc<RefCell<dyn BRepFillLocationLawOps>>,
    /// OCCT handle(BRepFill_SectionLaw) mySec.
    pub my_sec: Rc<RefCell<dyn crate::brep_fill::brep_fill_section_law::BRepFillSectionLawOps>>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myUEdges.
    pub my_u_edges: Option<ShapeHArray2>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myVEdges.
    pub my_v_edges: Option<ShapeHArray2>,
    /// OCCT NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> myVEdgesModified.
    pub my_v_edges_modified: std::collections::HashMap<ShapeKey, Shape>,
    /// OCCT handle(NCollection_HArray2<TopoDS_Shape>) myFaces.
    pub my_faces: Option<ShapeHArray2>,
    /// OCCT NCollection_List<TopoDS_Shape> myAuxShape.
    pub my_aux_shape: Vec<Shape>,
    /// OCCT handle(NCollection_HArray1<TopoDS_Shape>) myTapes.
    pub my_tapes: Option<Vec<Shape>>,
    /// OCCT double Error.
    pub error: f64,
    /// OCCT TopoDS_Wire FirstShape.
    pub first_shape: Shape,
    /// OCCT TopoDS_Wire LastShape.
    pub last_shape: Shape,
    /// The BuildWire v-close flag (the wire of cxx L2094-2096 — the
    /// wire under construction, carried for the pool-based builder).
    pub(crate) brep_wire: Option<Shape>,
}

impl BRepFillSweep {
    /// OCCT BRepFill_Sweep(Section, Location, WithKPart) (cxx L1940-1958).
    pub fn new(
        brep: &mut BRep,
        section: Rc<RefCell<dyn crate::brep_fill::brep_fill_section_law::BRepFillSectionLawOps>>,
        location: Rc<RefCell<dyn BRepFillLocationLawOps>>,
        with_kpart: bool,
    ) -> Self {
        let mut this = BRepFillSweep {
            is_done: false,
            kpart: with_kpart,
            my_tol3d: 0.0,
            my_bound_tol: 0.0,
            my_tol2d: 0.0,
            my_tol_angular: 0.0,
            my_ang_min: 0.0,
            my_ang_max: 0.0,
            my_approx_style: GeomFillApproxStyle::GeomFill_Location,
            my_continuity: GeomAbsShape::C2,
            my_degmax: 11,
            my_segmax: 30,
            my_force_approx_c1: false,
            my_shape: Shape::null(),
            my_loc: location,
            my_sec: section,
            my_u_edges: None,
            my_v_edges: None,
            my_v_edges_modified: std::collections::HashMap::new(),
            my_faces: None,
            my_aux_shape: Vec::new(),
            my_tapes: None,
            error: 0.0,
            first_shape: Shape::null(),
            last_shape: Shape::null(),
            brep_wire: None,
        };

        this.set_tolerance(brep, 1.0e-4, 1.0, 1.0e-5, 1.0e-2);
        this.set_angular_control(0.01, 6.0);
        // myAuxShape.Clear();
        this.my_aux_shape.clear();

        this.my_approx_style = GeomFillApproxStyle::GeomFill_Location;
        this.my_continuity = GeomAbsShape::C2;
        this.my_degmax = 11;
        this.my_segmax = 30;
        this.my_force_approx_c1 = false;
        this
    }

    /// OCCT SetBounds (cxx L1962-2007) — It is necessary to check the
    /// SameRange on its (PRO13551).  `BRepLib::CheckSameRange(E)` maps to
    /// the stored `TEdgeData::same_range` flag.
    pub fn set_bounds(&mut self, brep: &mut BRep, first: &Shape, last: &Shape) {
        self.first_shape = first.clone();
        self.last_shape = last.clone();

        let mut b = BRepBuilder::new();
        if !self.first_shape.is_null() {
            for e in crate::brep_fill::compatible_wires::wire_edges(brep, &self.first_shape) {
                let same_range = brep.edge(e.clone()).same_range;
                if !same_range {
                    b.set_edge_same_range(brep, e.clone(), false);
                    b.set_edge_same_parameter(brep, e.clone(), false);
                }
            }
        }

        if !self.last_shape.is_null() {
            for e in crate::brep_fill::compatible_wires::wire_edges(brep, &self.last_shape) {
                let same_range = brep.edge(e.clone()).same_range;
                if !same_range {
                    b.set_edge_same_range(brep, e.clone(), false);
                    b.set_edge_same_parameter(brep, e.clone(), false);
                }
            }
        }
    }

    /// OCCT SetTolerance (cxx L2011-2020).
    pub fn set_tolerance(
        &mut self,
        _brep: &mut BRep,
        tol3d: f64,
        bound_tol: f64,
        tol2d: f64,
        tol_angular: f64,
    ) {
        self.my_tol3d = tol3d;
        self.my_bound_tol = bound_tol;
        self.my_tol2d = tol2d;
        self.my_tol_angular = tol_angular;
    }

    /// OCCT SetAngularControl (cxx L2024-2028).
    pub fn set_angular_control(&mut self, min_angle: f64, max_angle: f64) {
        self.my_ang_min = min_angle.max(rcad_kernel::core::precision::ANGULAR);
        self.my_ang_max = max_angle.min(6.28);
    }

    /// OCCT SetForceApproxC1 (cxx L2036-2039).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_force_approx_c1 = force_approx_c1;
    }

    /// OCCT BuildWire (cxx L2078-2209) — the vertex-section sweep: one
    /// GeomFill_Sweep per path edge, the iso edge built over each surface.
    pub fn build_wire(&mut self, brep: &mut BRep, _transition: BRepFillTransitionStyle) -> bool {
        let mut isec = 1i32;
        let mut p1;

        let mut b = BRepBuilder::new();
        let nb_path = self.my_loc.borrow().base().nb_law();
        // vclose = (myLoc->IsClosed() && (myLoc->IsG1(0, myTol3d) >= 0));
        let vclose = self.my_loc.borrow().base().is_closed(brep)
            && self.my_loc.borrow().base().is_g1(brep, 0, self.my_tol3d, 1.0e-4) >= 0;
        self.error = 0.0;

        // TopoDS_Wire wire; B.MakeWire(wire);
        let wire = b.make_wire(brep);

        // (1.1) Construction of Tables
        // myFaces  = new HArray2(1, 1, 1, NbPath);
        // myUEdges = new HArray2(1, 2, 1, NbPath);
        // myVEdges = new HArray2(1, 1, 1, NbPath + 1);
        self.my_faces = Some(vec![vec![Shape::null(); nb_path as usize]]);
        self.my_u_edges = Some(vec![vec![Shape::null(); nb_path as usize]; 2]);
        self.my_v_edges = Some(vec![vec![Shape::null(); nb_path as usize + 1]]);

        // (1.2) Calculate curves / vertex / edge
        for ipath in 1..=nb_path {
            // GeomFill_Sweep Sweep(myLoc->Law(ipath), KPart);
            let loc_law = self.my_loc.borrow().base().law(ipath);
            let mut sweep =
                crate::geomalgo::geomfill::sweep::Sweep::new(loc_law, self.kpart);
            sweep.set_tolerance_sweep(
                self.my_tol3d,
                self.my_bound_tol,
                self.my_tol2d,
                self.my_tol_angular,
            );
            sweep.set_force_approx_c1(self.my_force_approx_c1);
            let sec_law = self.my_sec.borrow().base().law(isec);
            sweep.build(
                sec_law,
                self.my_approx_style,
                self.my_continuity,
                self.my_degmax,
                self.my_segmax,
            );
            if !sweep.is_done() {
                return false;
            }
            let surf = sweep.surface().expect("null surface").clone();
            // S->Bounds(...) and the iso parameter (cxx L2119-2141).
            let (_basis, bounds) = surface_adaptor_basis_and_bounds(&surf);
            let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
            let (mut first, mut last, val);
            let iso: Curve3;
            if sweep.exchange_uv() {
                if sweep.u_reversed() {
                    first = umin;
                    last = umax;
                    val = vmax;
                } else {
                    first = umin;
                    last = umax;
                    val = vmin;
                }
                iso = surface_viso(&surf, val);
            } else {
                if sweep.u_reversed() {
                    val = umax;
                    first = vmin;
                    last = vmax;
                } else {
                    val = umin;
                    first = vmin;
                    last = vmax;
                }
                iso = surface_uiso(&surf, val);
            }

            // Vertex by position (cxx L2143-2189)
            let my_v_edges = self.my_v_edges.as_mut().expect("myVEdges");
            if ipath < nb_path {
                let mut v = my_v_edges[0][ipath as usize].clone();
                build_vertex(brep, &iso, false, first, last, &mut v);
                my_v_edges[0][ipath as usize] = v;
            } else if vclose {
                let mut v = my_v_edges[0][0].clone();
                my_v_edges[0][ipath as usize] = v.clone();
                // Iso->D0(Last, P1); Tol = P1.Distance(BRep_Tool::Pnt(V));
                p1 = iso.point_at(last);
                let tol = p1.distance(brep.vertex(v.clone()).point);
                b.update_vertex_tolerance(brep, v, tol);
            } else if !self.first_shape.is_null() {
                // myVEdges->SetValue(1, NbPath, FirstShape);
                my_v_edges[0][nb_path as usize - 1] = self.first_shape.clone();
            } else {
                let mut v = my_v_edges[0][nb_path as usize].clone();
                build_vertex(brep, &iso, false, first, last, &mut v);
                my_v_edges[0][nb_path as usize] = v;
            }

            if ipath > 1 {
                p1 = iso.point_at(first);
                let v = my_v_edges[0][(ipath - 1) as usize].clone();
                let tol = p1.distance(brep.vertex(v.clone()).point);
                b.update_vertex_tolerance(brep, v, tol);
            }
            if ipath == 1 {
                if !self.first_shape.is_null() {
                    my_v_edges[0][0] = self.first_shape.clone();
                } else {
                    let mut v = my_v_edges[0][0].clone();
                    build_vertex(brep, &iso, true, first, last, &mut v);
                    my_v_edges[0][0] = v;
                }
            }

            // Construction of the edge (cxx L2192-2205)
            // BRepLib_MakeEdge MkE; MkE.Init(Iso, V(1, ipath), V(1, ipath+1),
            //                                Iso->FirstParameter(),
            //                                Iso->LastParameter());
            let v_edges = self.my_v_edges.as_ref().expect("myVEdges");
            let v1 = v_edges[0][(ipath - 1) as usize].clone();
            let v2 = v_edges[0][ipath as usize].clone();
            let e = b.add_edge(
                brep,
                Some(iso.clone()),
                v1,
                v2,
                [
                    curve_first_parameter_of(&iso),
                    curve_last_parameter_of(&iso),
                ],
            );
            // if (!MkE.IsDone()) return false; — the rcad add_edge is
            // infallible at this layer (the OCCT error arm is empty too).
            b.update_edge_tolerance(brep, e.clone(), sweep.error_on_surface());
            b.add_to_wire(brep, wire.clone(), e.clone());
            self.my_faces.as_mut().expect("myFaces")[0][(ipath - 1) as usize] = e;
        }
        self.my_shape = wire;
        true
    }

    /// OCCT MergeVertex (cxx L3984-4004) — Make V2 = V1 if V2 is too close
    /// to V1.
    pub fn merge_vertex(&self, brep: &BRep, v1: &Shape, v2: &mut Shape) -> bool {
        let mut tol = brep.vertex(v1.clone()).tolerance.max(brep.vertex(v2.clone()).tolerance);
        if tol < self.my_tol3d {
            tol = self.my_tol3d;
        }
        let p1 = brep.vertex(v1.clone()).point;
        let p2 = brep.vertex(v2.clone()).point;
        if p1.distance(p2) <= tol {
            *v2 = v1.clone();
            return true;
        }
        false
    }

    /// OCCT UpdateVertex (cxx L4010-4040) — Update the Tolerance of
    /// Vertices depending on Laws.
    pub fn update_vertex(
        &self,
        brep: &mut BRep,
        ipath: i32,
        isec: i32,
        err_app: f64,
        param: f64,
        v: &mut Shape,
    ) {
        let mut b = BRepBuilder::new();
        let the_v = v.clone();
        let mut vv = Shape::null();
        let (sec_vertex, sec_tol) = {
            let sec = self.my_sec.borrow();
            let v = sec.vertex(brep, isec, param);
            let t = sec.vertex_tol(brep, isec - 1, param);
            (v, t)
        };
        self.my_loc
            .borrow()
            .base()
            .perform_vertex(brep, ipath, &sec_vertex, err_app + sec_tol, &mut vv, 0);
        let p1 = brep.vertex(vv.clone()).point;
        let p2 = brep.vertex(the_v.clone()).point;

        let mut tol = brep.vertex(vv.clone()).tolerance;
        tol += p1.distance(p2);

        if tol > brep.vertex(the_v.clone()).tolerance {
            b.update_vertex_tolerance(brep, the_v, tol);
        }
    }

    /// OCCT EvalExtrapol (cxx L3910-3978).  (the &mut receiver carries the
    /// CurvilinearBounds myLength mutation through the OCCT array handle —
    /// the const-mutation pattern documented in brep_fill_location_law.rs)
    pub fn eval_extrapol(
        &mut self,
        brep: &mut BRep,
        index: i32,
        transition: BRepFillTransitionStyle,
    ) -> f64 {
        let mut extrap = 0.0f64;
        if transition == BRepFillTransitionStyle::RightCorner {
            let i1;
            let i2;
            if (index == 1) || (index == self.my_loc.borrow().base().nb_law() + 1) {
                if !self.my_loc.borrow().base().is_closed(brep)
                    || !self.my_sec.borrow().base().is_vclosed()
                {
                    return extrap;
                }
                i1 = self.my_loc.borrow().base().nb_law();
                i2 = 1;
            } else {
                i1 = index - 1;
                i2 = index;
            }

            let mut v1 = DVec3::ZERO;
            let mut v2 = DVec3::ZERO;
            let mut m1 = crate::geomalgo::geomfill::gp_mat::GpMat::identity();
            let mut m2 = crate::geomalgo::geomfill::gp_mat::GpMat::identity();

            let mut f = 0.0;
            let mut l = 0.0;
            self.my_loc.borrow().base().law(i1).borrow().get_domain(&mut f, &mut l);
            self.my_loc.borrow().base().law(i1).borrow().d0(l, &mut m1, &mut v1);
            let t1 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m1, 3);
            self.my_loc.borrow().base().law(i2).borrow().get_domain(&mut f, &mut l);
            self.my_loc.borrow().base().law(i2).borrow().d0(f, &mut m2, &mut v2);
            let t2 = crate::brep_fill::brep_fill_location_law::gp_mat_column(&m2, 3);

            let alpha = gv_angle(t1, t2);
            if (alpha > self.my_ang_max) || (alpha < self.my_ang_min) {
                // Angle too great => No "straight" connection
                // Angle too small => No connection
                return extrap; // = 0.0
            }

            // Sec = mySec->ConcatenedLaw();
            let sec = self.my_sec.borrow().concatened_law(brep);
            let Some(sec) = sec else { return extrap };

            // Calculating parameter U (cxx L3954-3960)
            let (mut lf, mut length) = (0.0f64, 0.0f64);
            self.my_loc.borrow_mut().base_mut().curvilinear_bounds(
                self.my_loc.borrow().base().nb_law(),
                &mut lf,
                &mut length,
            );
            let (mut sec_first, mut sec_len) = (0.0f64, 0.0f64);
            self.my_sec
                .borrow()
                .base()
                .law(1)
                .borrow()
                .get_domain(&mut sec_first, &mut sec_len);
            sec_len -= sec_first;
            let (mut lf2, mut ll) = (0.0f64, 0.0f64);
            self.my_loc
                .borrow_mut()
                .base_mut()
                .curvilinear_bounds(i1, &mut lf2, &mut ll);
            let u = sec_first + (ll / length) * sec_len;

            // Box(Sec, U, box); box.Get(Xmin, ..., Zmax);
            let (min, max) = box_of_section(&sec, u);

            let mut r = min[0]
                .abs()
                .max(max[0].abs())
                .max(min[1].abs().max(max[1].abs()));
            // double coef = 2.;
            let coef = 2.0;
            r *= coef;
            extrap = min[2].abs().max(max[2].abs()) + 100.0 * self.my_tol3d;
            extrap += r * (alpha / 2.0).tan();
        }
        extrap
    }

    // -------------------------------------------------------------------
    // Accessors (cxx L3534-3583)
    // -------------------------------------------------------------------

    /// OCCT IsDone (cxx L3534-3537).
    pub fn is_done(&self) -> bool {
        self.is_done
    }

    /// OCCT Shape (cxx L3541-3544).
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT ErrorOnSurface (cxx L3548-3551).
    pub fn error_on_surface(&self) -> f64 {
        self.error
    }

    /// OCCT SubShape (cxx L3555-3558) — myFaces.
    pub fn sub_shape(&self) -> Option<ShapeHArray2> {
        self.my_faces.clone()
    }

    /// OCCT InterFaces (cxx L3562-3565) — myUEdges.
    pub fn inter_faces(&self) -> Option<ShapeHArray2> {
        self.my_u_edges.clone()
    }

    /// OCCT Sections (cxx L3571-3574) — myVEdges.
    pub fn sections(&self) -> Option<ShapeHArray2> {
        self.my_v_edges.clone()
    }

    /// OCCT Tape (cxx L3580-3583).
    pub fn tape(&self, index: i32) -> Shape {
        self.my_tapes
            .as_ref()
            .map(|t| t[(index - 1) as usize].clone())
            .unwrap_or_else(Shape::null)
    }
}

/// OCCT gp_Vec::Angle (gp_XYZ::Angle).
fn gv_angle(v1: DVec3, v2: DVec3) -> f64 {
    let an_norm = v1.length();
    let a_no_norm = v2.length();
    let mut value = v1.dot(v2) / (an_norm * a_no_norm);
    if value > 1.0 {
        value = 1.0;
    } else if value < -1.0 {
        value = -1.0;
    }
    value.acos()
}

/// OCCT Geom_Curve::FirstParameter over the rcad Curve3.
fn curve_first_parameter_of(c: &Curve3) -> f64 {
    crate::geomalgo::geomfill::trihedron_law::curve_first_parameter(c)
}

/// OCCT Geom_Curve::LastParameter over the rcad Curve3.
fn curve_last_parameter_of(c: &Curve3) -> f64 {
    crate::geomalgo::geomfill::trihedron_law::curve_last_parameter(c)
}

// The gp::Resolution threshold consumed by the re-hosts.
const _: f64 = GP_RESOLUTION;

// The unused-import guards for the part-B symbols consumed from this module.
#[allow(unused_imports)]
use crate::brep_fill::brep_fill_pipe_shell_b::ShapeToArray2Map as _PartBShapeToArray2Map;
#[allow(dead_code)]
fn _part_b_type_guards(_t: Option<TrimmedCurve3>, _s: Option<HashSet<ShapeKey>>) {}
