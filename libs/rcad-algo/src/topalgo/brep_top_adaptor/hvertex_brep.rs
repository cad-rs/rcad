//! BRepTopAdaptor_HVertex (TKTopAlgo) — the topological vertex of the
//! BRepTopAdaptor_TopolTool, with the parametric-resolution refinement.
//!
//! 1:1 translation of OCCT `BRepTopAdaptor_HVertex.hxx` (L17-57) + `.cxx`
//! (L24-192).

use glam::DVec2;
use rcad_kernel::geom::{CurveEval, SurfaceEval};
use glam::DVec3;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, Shape};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::topalgo::adaptor3d::HVertex;
use crate::topalgo::brep_adaptor::curve2d::BRepCurve2d;

/// OCCT BRepTopAdaptor_HVertex — a TopoDS_Vertex bound to its pcurve
/// adaptor.
pub struct BRepHVertex<'c, 'a> {
    /// OCCT Adaptor3d_HVertex base subobject (never initialized by the BRep
    /// subclass ctor; Value/Resolution/IsSame are all overridden).
    pub(crate) base: HVertex,
    /// OCCT myVtx.
    my_vtx: Shape,
    /// OCCT myCurve.
    my_curve: &'c BRepCurve2d<'a>,
    /// The owning BRep (BRep_Tool accessors).
    brep: &'c BRep,
}

impl<'c, 'a> BRepHVertex<'c, 'a> {
    /// OCCT BRepTopAdaptor_HVertex(V, C) (cxx L26-30).
    pub fn new(v: &Shape, c: &'c BRepCurve2d<'a>) -> Self {
        BRepHVertex {
            base: HVertex::new(),
            my_vtx: v.clone(),
            my_curve: c,
            brep: c.brep(),
        }
    }

    /// OCCT Vertex() (hxx L34).
    pub fn vertex(&self) -> &Shape {
        &self.my_vtx
    }

    /// OCCT Value() (cxx L32-36) — "do nothing": always (RealFirst,
    /// RealFirst).
    pub fn value(&self) -> DVec2 {
        DVec2::new(-f64::MAX, -f64::MAX)
    }

    /// OCCT Parameter(C) (cxx L38-43) — BRep_Tool::Parameter(myVtx,
    /// brhc->Edge(), brhc->Face()).
    pub fn parameter(&self, _c: &BRepCurve2d<'a>) -> f64 {
        self.brep
            .parameter_on_edge(&self.my_vtx, self.my_curve.edge(), self.my_curve.face())
            .unwrap_or(0.0)
    }

    /// OCCT Resolution(C) (cxx L45-170) — the parametric resolution of the
    /// vertex on its pcurve, refined against the face surface.
    pub fn resolution(&self, c: &BRepCurve2d<'a>) -> f64 {
        let f = c.face();
        // OCCT BRepAdaptor_Surface S(F, false) — the face surface; the
        // pcurve and the surface share the face-local frame.
        let surf = self
            .brep
            .face_surface(f)
            .expect("BRepAdaptor_Surface: face has no surface");

        // double tv = BRep_Tool::Tolerance(myVtx);
        let tv = self.brep.vertex_tolerance(&self.my_vtx);
        // double p = BRep_Tool::Parameter(myVtx, brhc->Edge(), brhc->Face());
        let p = self.parameter(c);
        let or = self.orientation();
        // gp_Pnt P, P1; gp_Vec DU, DV, DC;
        // C->D1(p, p2d, v2d);
        let (p2d, v2d) = Curve2dAdaptor::d1(c, p);
        // S.D1(p2d.X(), p2d.Y(), P, DU, DV);
        let (big_p, du, dv) = surf.derivatives(p2d.x, p2d.y);
        // DC.SetLinearForm(v2d.X(), DU, v2d.Y(), DV);
        let dc = v2d.x * du + v2d.y * dv;
        let mut res_uv;
        let mag = dc.length();

        // double URes = S.UResolution(tv); double VRes = S.VResolution(tv);
        let u_res = self.brep.u_resolution(f, tv);
        let v_res = self.brep.v_resolution(f, tv);
        // double tURes = C->Resolution(URes); double tVRes = C->Resolution(VRes);
        let t_u_res = Curve2dAdaptor::resolution(c, u_res);
        let t_v_res = Curve2dAdaptor::resolution(c, v_res);
        let res_uv1 = t_u_res.max(t_v_res);

        if mag < 1e-12 {
            return res_uv1;
        }

        // For lack of better options limit the parametric solution to
        // 10 million * tolerance of the point.
        if tv > 1.0e7 * mag {
            res_uv = 1.0e7;
        } else {
            res_uv = tv / mag;
        }

        let u_min = Curve2dAdaptor::first_parameter(c);
        let u_max = Curve2dAdaptor::last_parameter(c);
        let point_at = |pp: f64| -> DVec3 {
            let p2 = Curve2dAdaptor::value(c, pp);
            surf.point_at(p2.x, p2.y)
        };

        // Control
        let mut pp = if or == Orientation::Reversed {
            p + res_uv
        } else {
            p - res_uv
        };
        pp = pp.clamp(u_min, u_max);

        // C->D0(pp, p2d); S.D0(p2d.X(), p2d.Y(), P1);
        let p1 = point_at(pp);

        let mut dist = big_p.distance(p1);
        if (dist > 1e-12) && ((dist > 1.1 * tv) || (dist < 0.8 * tv)) {
            // Refine if possible
            let pp = if or == Orientation::Reversed {
                p + tv / dist
            } else {
                p - tv / dist
            };
            let pp = pp.clamp(u_min, u_max);

            let (p2d, v2d) = Curve2dAdaptor::d1(c, pp);
            let (p1, du, dv) = surf.derivatives(p2d.x, p2d.y);
            let dc = v2d.x * du + v2d.y * dv;
            let mut dist1 = big_p.distance(p1);
            if (dist1 - tv).abs() < (dist - tv).abs() {
                // Take the result of interpolation
                res_uv = tv / dist;
                dist = dist1;
            }

            let mut mag = dc.length();
            if tv > 1.0e7 * mag {
                mag = tv * 1.0e-7;
            }
            let pp = if or == Orientation::Reversed {
                p + tv / mag
            } else {
                p - tv / mag
            };
            let pp = pp.clamp(u_min, u_max);

            let p1 = point_at(pp);
            dist1 = big_p.distance(p1);
            if (dist1 - tv).abs() < (dist - tv).abs() {
                // Take the new estimation
                res_uv = tv / mag;
                let _ = dist1;
            }
        }

        res_uv.min(res_uv1)
    }

    /// OCCT Orientation() (cxx L172-175).
    pub fn orientation(&self) -> Orientation {
        self.my_vtx.orientation
    }

    /// OCCT IsSame(Other) (cxx L177-182) — downcast + TopoDS IsSame (same
    /// TShape and location).
    pub fn is_same(&self, other: &BRepHVertex<'c, 'a>) -> bool {
        self.my_vtx.is_same(&other.my_vtx)
    }
}
