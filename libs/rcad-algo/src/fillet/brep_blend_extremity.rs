//! OCCT BRepBlend_Extremity (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_Extremity.hxx (L36-155) + BRepBlend_Extremity.cxx (L19-142) +
//! BRepBlend_Extremity.lxx (L8-87).
//!
//! Architecture mapping: `occ::handle<Adaptor3d_HVertex>` maps to an owned
//! rcad [`HVertex`] value; `occ::handle<Adaptor2d_Curve2d>` maps to a cloned
//! rcad `Curve2d` value (rcad has no shared-handle adaptor layer).

use glam::DVec3;

use rcad_kernel::geom::Curve2d;

use crate::fillet::brep_blend_point_on_rst::BRepBlendPointOnRst;
use crate::geomalgo::int_patch::transitions::Transition as IntSurfTransition;
use crate::topalgo::adaptor3d::hvertex::HVertex;

/// OCCT BRepBlend_Extremity — definition of an extremity of a 3d curve on a
/// surface or on a restriction (BRepBlend_Extremity.hxx L30).
#[derive(Debug, Clone)]
pub struct BRepBlendExtremity {
    vtx: Option<HVertex>,
    seqpt: Vec<BRepBlendPointOnRst>,
    pt: DVec3,
    tang: DVec3,
    param: f64,
    u: f64,
    v: f64,
    tol: f64,
    isvtx: bool,
    hastang: bool,
}

impl BRepBlendExtremity {
    /// OCCT BRepBlend_Extremity() (BRepBlend_Extremity.cxx L19-28) —
    /// empty constructor.
    pub fn new() -> Self {
        BRepBlendExtremity {
            vtx: None,
            seqpt: Vec::new(),
            pt: DVec3::ZERO,
            tang: DVec3::ZERO,
            param: 0.0,
            u: 0.0,
            v: 0.0,
            tol: 0.0,
            isvtx: false,
            hastang: false,
        }
    }

    /// OCCT BRepBlend_Extremity(P, U, V, Param, Tol) (L30-40).
    pub fn new_on_surface(p: DVec3, u: f64, v: f64, param: f64, tol: f64) -> Self {
        BRepBlendExtremity {
            vtx: None,
            seqpt: Vec::new(),
            pt: p,
            tang: DVec3::ZERO,
            param,
            u,
            v,
            tol,
            isvtx: false,
            hastang: false,
        }
    }

    /// OCCT BRepBlend_Extremity(P, U, V, Param, Tol, Vtx) (L42-53).
    pub fn new_on_surface_with_vertex(
        p: DVec3,
        u: f64,
        v: f64,
        param: f64,
        tol: f64,
        vtx: &HVertex,
    ) -> Self {
        BRepBlendExtremity {
            vtx: Some(vtx.clone()),
            seqpt: Vec::new(),
            pt: p,
            tang: DVec3::ZERO,
            param,
            u,
            v,
            tol,
            isvtx: true,
            hastang: false,
        }
    }

    /// OCCT BRepBlend_Extremity(P, W, Param, Tol) (L55-64).
    pub fn new_on_curve(p: DVec3, w: f64, param: f64, tol: f64) -> Self {
        BRepBlendExtremity {
            vtx: None,
            seqpt: Vec::new(),
            pt: p,
            tang: DVec3::ZERO,
            param,
            u: w,
            v: 0.0,
            tol,
            isvtx: false,
            hastang: false,
        }
    }

    /// OCCT SetValue(P, U, V, Param, Tol) (BRepBlend_Extremity.cxx L66-77).
    pub fn set_value(&mut self, p: DVec3, u: f64, v: f64, param: f64, tol: f64) {
        self.pt = p;
        self.u = u;
        self.v = v;
        self.param = param;
        self.tol = tol;
        self.isvtx = false;
        self.seqpt.clear();
    }

    /// OCCT SetValue(P, U, V, Param, Tol, Vtx) (L79-91).
    pub fn set_value_with_vertex(
        &mut self,
        p: DVec3,
        u: f64,
        v: f64,
        param: f64,
        tol: f64,
        vtx: &HVertex,
    ) {
        self.pt = p;
        self.u = u;
        self.v = v;
        self.param = param;
        self.tol = tol;
        self.isvtx = true;
        self.vtx = Some(vtx.clone());
        self.seqpt.clear();
    }

    /// OCCT SetValue(P, W, Param, Tol) (L93-103).
    pub fn set_value_on_curve(&mut self, p: DVec3, w: f64, param: f64, tol: f64) {
        self.pt = p;
        self.u = w;
        self.param = param;
        self.tol = tol;
        self.isvtx = false;
        self.seqpt.clear();
    }

    /// OCCT SetVertex(V) (L105-109).
    pub fn set_vertex(&mut self, v: &HVertex) {
        self.isvtx = true;
        self.vtx = Some(v.clone());
    }

    /// OCCT AddArc(A, Param, TLine, TArc) (L111-117).
    pub fn add_arc(
        &mut self,
        a: &Curve2d,
        param: f64,
        t_line: IntSurfTransition,
        t_arc: IntSurfTransition,
    ) {
        self.seqpt
            .push(BRepBlendPointOnRst::new_with_arc(a, param, t_line, t_arc));
    }

    /// OCCT SetTangent(Tangent) (BRepBlend_Extremity.lxx L9-13).
    pub fn set_tangent(&mut self, tangent: DVec3) {
        self.hastang = true;
        self.tang = tangent;
    }

    /// OCCT Value() (lxx accessors return pt).
    pub fn value(&self) -> DVec3 {
        self.pt
    }

    /// OCCT HasTangent() (lxx L15-18).
    pub fn has_tangent(&self) -> bool {
        self.hastang
    }

    /// OCCT Tangent() (lxx L20-27) — throws Standard_DomainError when the
    /// tangent was not set.
    pub fn tangent(&self) -> DVec3 {
        if !self.hastang {
            panic!("Standard_DomainError: BRepBlend_Extremity::Tangent");
        }
        self.tang
    }

    /// OCCT Parameters(U, V) (lxx L29-33).
    pub fn parameters(&self) -> (f64, f64) {
        (self.u, self.v)
    }

    /// OCCT Tolerance() (lxx L35-38).
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// OCCT IsVertex() (lxx L40-43).
    pub fn is_vertex(&self) -> bool {
        self.isvtx
    }

    /// OCCT Vertex() (lxx L45-52) — throws Standard_DomainError when not a
    /// vertex.
    pub fn vertex(&self) -> &HVertex {
        if !self.isvtx {
            panic!("Standard_DomainError: BRepBlend_Extremity::Vertex");
        }
        self.vtx.as_ref().expect("isvtx without a vertex")
    }

    /// OCCT NbPointOnRst() (lxx L54-57).
    pub fn nb_point_on_rst(&self) -> i32 {
        self.seqpt.len() as i32
    }

    /// OCCT PointOnRst(Index) (lxx L59-62) — 1-based Index.
    pub fn point_on_rst(&self, index: i32) -> &BRepBlendPointOnRst {
        &self.seqpt[(index - 1) as usize]
    }

    /// OCCT Parameter() (lxx L64-69).
    pub fn parameter(&self) -> f64 {
        self.u
    }

    /// OCCT ParameterOnGuide() (lxx L71-76).
    pub fn parameter_on_guide(&self) -> f64 {
        self.param
    }
}
