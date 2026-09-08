//! OCCT Draft (Draft.hxx L28-41 + Draft.cxx L34-100) — the static draft
//! angle query.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/Draft.hxx
//!         $OCCT_SRC/src/ModelingAlgorithms/TKOffset/Draft/Draft.cxx
//!
//! Architecture differences:
//! 1. Geom_Surface::Transformed(Lo.Transformation()) (cxx L53) — rcad carries
//!    the location as an id and the draft flows are identity-location; the
//!    re-application of the location transformation is an identity (the
//!    loc_ope_split_drafts.rs #12 precedent).
//! 2. ElSLib::D1(U, V, Co, P, D1u, D1v) — the rcad stand-in evaluates the
//!    cone through SurfaceEval::derivatives (the GeomAdaptor vehicle,
//!    geom/eval.rs).

use glam::DVec3;
use rcad_kernel::geom::{Surface3, SurfaceEval};
use rcad_kernel::precision::ANGULAR;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::Orientation;

use super::draft_modification::{brep_tool_surface, brep_tools_uv_bounds};

/// OCCT Draft::Angle(const TopoDS_Face& F, const gp_Dir& Direction)
/// (Draft.hxx L33-40; Draft.cxx L34-100).
pub fn angle(the_f: &Shape, d: DVec3) -> f64 {
    // OCCT L37-38: TopLoc_Location Lo; S = BRep_Tool::Surface(F, Lo).
    let mut s = brep_tool_surface(the_f).expect("BRep_Tool::Surface");

    // OCCT L39-44: TypeS = S->DynamicType();
    // if (TypeS == STANDARD_TYPE(Geom_RectangularTrimmedSurface))
    //   { S = BasisSurface(); TypeS = S->DynamicType(); }
    if let Surface3::Trimmed(t) = s {
        s = t.basis.as_ref().clone();
    }

    // OCCT L46-50: only Plane / Conical / Cylindrical are valid.
    match &s {
        Surface3::Plane(_) | Surface3::Cone(_) | Surface3::Cylinder(_) => {}
        _ => {
            // OCCT: throw Standard_DomainError();
            panic!("Standard_DomainError");
        }
    }

    // OCCT L53: S = down_cast<Geom_Surface>(S->Transformed(Lo.Transformation()))
    // — identity (architecture difference #1).
    let s = s;

    let angle;
    match s {
        // OCCT L54-67: the plane case.
        Surface3::Plane(pl) => {
            // OCCT L56-57: gp_Ax3 ax3(Pln().Position()); gp_Vec normale(ax3.Direction());
            let mut normale = pl.normal;
            // OCCT L58-61: if (!ax3.Direct()) normale.Reverse(); — the direct
            // sense of the rcad plane frame is (U ^ V) . N > 0.
            let direct = pl.u_dir.cross(pl.v_dir).dot(pl.normal) > 0.0;
            if !direct {
                normale = -normale;
            }
            // OCCT L62-65: if (F.Orientation() == TopAbs_REVERSED) normale.Reverse();
            if the_f.orientation == Orientation::Reversed {
                normale = -normale;
            }
            // OCCT L66: Angle = std::asin(normale.Dot(D));
            angle = normale.dot(d).asin();
        }
        // OCCT L68-77: the cylinder case.
        Surface3::Cylinder(cy) => {
            // OCCT L71: testdir = D.Dot(Cy.Axis().Direction());
            let testdir = d.dot(cy.axis);
            if testdir.abs() <= 1.0 - ANGULAR {
                // OCCT L74: throw Standard_DomainError();
                panic!("Standard_DomainError");
            }
            // OCCT L76: Angle = 0.;
            angle = 0.0;
        }
        // OCCT L78-98: the cone case.
        Surface3::Cone(co) => {
            // OCCT L81-85: the axis must be colinear with D.
            let testdir = d.dot(co.axis);
            if testdir.abs() <= 1.0 - ANGULAR {
                // OCCT L84: throw Standard_DomainError();
                panic!("Standard_DomainError");
            }
            // OCCT L86-87: BRepTools::UVBounds(F, umin, umax, vmin, vmax);
            let bounds = brep_tools_uv_bounds(the_f);
            let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);
            // OCCT L90: ElSLib::D1(umin + umax / 2., vmin + vmax / 2., Co, ptbid, d1u, d1v);
            let (_, d1u_raw, d1v) = Surface3::Cone(co).derivatives(umin + umax / 2.0, vmin + vmax / 2.0);
            // OCCT L91-92: d1u.Cross(d1v); d1u.Normalize();
            let mut d1u = d1u_raw.cross(d1v).normalize_or_zero();
            // OCCT L93-96: if (F.Orientation() == TopAbs_REVERSED) d1u.Reverse();
            if the_f.orientation == Orientation::Reversed {
                d1u = -d1u;
            }
            // OCCT L97: Angle = std::asin(d1u.Dot(D));
            angle = d1u.dot(d).asin();
        }
        // The OCCT L46-50 domain guard already excluded the other types.
        _ => unreachable!(),
    }
    // OCCT L99: return Angle;
    angle
}
