//! OCCT BRepAdaptor_Surface / BRepAdaptor_Curve (TKBRep/BRepAdaptor) — the
//! 1:1 translation of `BRepAdaptor_Surface.hxx` (L48-83) + `.cxx`
//! (L113-180), `BRepAdaptor_Curve.hxx` (L37-98) + `.cxx` (L132-276), and the
//! BRepTools::UVBounds support (BRepTools.cxx L64-158, L170-362).
//!
//! OCCT inheritance `GeomAdaptor_TransformedSurface <- BRepAdaptor_Surface`
//! becomes composition: the rcad struct owns the
//! [`TransformedSurfaceAdaptor`] base (the OCCT inherited member set) plus
//! `myFace` / `myEdge`, and every inherited member forwards to the base.
//!
//! Architecture differences (rcad encodings):
//!   - the OCCT TopoDS_Shape embeds the TShape handle AND the location
//!     datum; rcad splits them — `Shape.location` indexes `BRep.locations` —
//!     so the adaptor members that read `BRep_Tool::Surface(F, L)` /
//!     `TopLoc_Location::Transformation()` take the owning `&BRep` next to
//!     the shape (the same convention as every rcad BRep_Tool read);
//!   - the `TopLoc_Location::Transformation()` gp_Trsf read is
//!     [`location_transformation`].

use std::sync::Arc;

use crate::math::bnd::BndBox2d;
use crate::math::gp::Trsf;
use crate::math::GeomAbsShape;
use crate::topo::brep_tool::brep_tool_tolerance;
use crate::topo::topods::{BRep, BRepTool as _, Orientation, Shape, TShape};

use super::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, Adaptor3dSurface, Curve2dHandle, CurveOnSurface,
    Geom2dCurveAdaptor, SurfaceHandle,
};
use super::geom_adaptor_surface::GeomSurfaceAdaptor;
use super::geom_adaptor_transformed_curve::TransformedCurveAdaptor;
use super::geom_adaptor_transformed_surface::TransformedSurfaceAdaptor;
use super::proj_lib_projected_curve::{Adaptor3dCurveGeom, Adaptor3dSurfaceGeom};
use crate::geom::{Curve2d, Surface3, SurfaceEval};

// =========================================================================
// OCCT TopLoc_Location::Transformation — the location Trsf read
// =========================================================================

/// OCCT `TopLoc_Location::Transformation()` — the gp_Trsf of the accumulated
/// location.  The rcad location slot is a `glam::DAffine3`; the encoding
/// splits it back into the gp_Trsf member layout (matrix / loc / scale,
/// gp_Trsf.hxx L379-385).
pub fn location_transformation(brep: &BRep, s: &Shape) -> Trsf {
    let loc = brep.get_location(s.location);
    let mut t = Trsf::identity();
    if loc == glam::DAffine3::IDENTITY {
        return t;
    }
    let scale = loc.x_axis.length();
    let cols = [loc.x_axis, loc.y_axis, loc.z_axis];
    for (j, col) in cols.iter().enumerate() {
        t.matrix[0][j] = col.x / scale;
        t.matrix[1][j] = col.y / scale;
        t.matrix[2][j] = col.z / scale;
    }
    t.loc = loc.w_axis;
    t.scale = scale;
    t
}

// OCCT `BRep_Tool::Tolerance(S)` now has a single kernel-local definition,
// [`crate::topo::brep_tool::brep_tool_tolerance`] (re-exported from the
// BRep_Tool home `crate::topo::topods`); the duplicate that used to live here
// was deleted.

// =========================================================================
// OCCT BRepTools::UVBounds (BRepTools.cxx L64-158 + the per-edge body
// L170-362)
// =========================================================================

/// OCCT `BRepTools::UVBounds(F, UMin, UMax, VMin, VMax)` (BRepTools.cxx
/// L64-77): Bnd_Box2d + AddUVBounds + the void answer.
pub fn brep_tools_uv_bounds(brep: &BRep, f: &Shape) -> (f64, f64, f64, f64) {
    // OCCT L66: Bnd_Box2d B; AddUVBounds(F, B).
    let mut b = BndBox2d::new();
    brep_tools_add_uv_bounds(brep, f, &mut b);
    // OCCT L67-73: if (!B.IsVoid()) B.Get(...) else all zero.
    match b.get() {
        Some((u_min, v_min, u_max, v_max)) => (u_min, u_max, v_min, v_max),
        None => (0.0, 0.0, 0.0, 0.0),
    }
}

/// OCCT `BRepTools::AddUVBounds(const TopoDS_Face& FF, Bnd_Box2d& B)`
/// (BRepTools.cxx L123-158).
pub fn brep_tools_add_uv_bounds(brep: &BRep, ff: &Shape, b: &mut BndBox2d) {
    // OCCT L125-126: F = FF; F.Orientation(TopAbs_FORWARD).
    let mut f = ff.clone();
    f.orientation = Orientation::Forward;
    // OCCT L128: TopExp_Explorer ex(F, TopAbs_EDGE) — the edge walk.
    // OCCT L132-136: fill box for the given face.
    let mut a_box = BndBox2d::new();
    for e in face_edge_shapes(&f) {
        brep_tools_add_uv_bounds_edge(brep, &f, &e, &mut a_box);
    }

    // OCCT L139-154: if the box is empty (face without edges or without
    // pcurves), get natural bounds.
    if a_box.is_void() {
        // OCCT L142-146: aSurf = BRep_Tool::Surface(F, L); null -> return.
        let Some(a_surf) = brep.face_surface(&f) else {
            return;
        };
        // OCCT L148: aSurf->Bounds(UMin, UMax, VMin, VMax).
        let [u_min, u_max, v_min, v_max] = SurfaceEval::default_domain(a_surf);
        // OCCT L150: aBox.Update(UMin, VMin, UMax, VMax).
        a_box.update(u_min, v_min, u_max, v_max);
    }

    // OCCT L156: B.Add(aBox).
    b.add_box(&a_box);
}

/// OCCT `BRepTools::AddUVBounds(const TopoDS_Face& aF, const TopoDS_Edge&
/// aE, Bnd_Box2d& aB)` (BRepTools.cxx L170-362).
fn brep_tools_add_uv_bounds_edge(brep: &BRep, a_f: &Shape, a_e: &Shape, a_b: &mut BndBox2d) {
    //
    // OCCT L179: aC2D = BRep_Tool::CurveOnSurface(aE, aF, aT1, aT2).
    let Some((c2d, a_t1, a_t2)) = brep.curve_on_surface(a_e, a_f) else {
        // OCCT L180-182: null pcurve -> return.
        return;
    };
    //
    // OCCT L191: aS = BRep_Tool::Surface(aF, aLoc).
    let Some(a_s) = brep.face_surface(a_f) else {
        return;
    };
    // OCCT L185-360 (the rest of the body) is data-source independent.
    let a_box_s = brep_tools_add_uv_bounds_curve_box(&c2d, a_t1, a_t2, a_s);
    // OCCT L360: aB.Add(aBoxS).
    a_b.add_box(&a_box_s);
}

/// The data-source independent remainder of `BRepTools::AddUVBounds(aF, aE,
/// aB)` (BRepTools.cxx L185-360): the 2D box of the edge pcurve, clamped to
/// the face's UV domain by the U (L203-299) and V (L302-355) periodicity
/// rules.  OCCT obtains `c2d` from `BRep_Tool::CurveOnSurface` and `a_s` from
/// `BRep_Tool::Surface`; those two lookups are the only steps that depend on
/// which shape pool the face lives in, so callers supply their results.
pub fn brep_tools_add_uv_bounds_curve_box(
    c2d: &Curve2d,
    a_t1: f64,
    a_t2: f64,
    a_s: &Surface3,
) -> BndBox2d {
    // OCCT L173-176: the scalars and the scratch boxes.
    let mut a_box_s = BndBox2d::new();
    //
    // OCCT L185: BndLib_Add2dCurve::Add(aC2D, aT1, aT2, 0., aBoxC) followed
    // by L186-188 `if (!aBoxC.IsVoid()) aBoxC.Get(...)` — the kernel
    // curve2d_bounding_box answers the box directly.
    let a_box_c = crate::curve2d_bounding_box(c2d, a_t1, a_t2, 0.0);
    let (mut a_x_min, mut a_y_min, mut a_x_max, mut a_y_max) =
        (a_box_c[0], a_box_c[2], a_box_c[1], a_box_c[3]);
    //
    // OCCT L192: aS->Bounds(aUmin, aUmax, aVmin, aVmax).
    let [a_umin, a_umax, a_vmin, a_vmax] = SurfaceEval::default_domain(a_s);

    // OCCT L194-200: unwrap one Geom_RectangularTrimmedSurface level.
    let a_s: &Surface3 = match a_s {
        Surface3::Trimmed(t) => t.basis.as_ref(),
        s => s,
    };

    // OCCT L203: if (!aS->IsUPeriodic()).
    if !a_s.is_u_periodic() {
        let mut is_u_periodic = false;

        // OCCT L210-219: the additional U-periodicity verification for the
        // B-spline surfaces (two- then three-or-six-point sampling).
        if matches!(a_s, Surface3::BSpline(_)) && (a_x_min < a_umin || a_x_max > a_umax) {
            let a_tol2 =
                100.0 * crate::core::precision::CONFUSION * crate::core::precision::CONFUSION;
            is_u_periodic = true;
            // OCCT L222-235: 1. verify U-closedness by sampling.
            if !a_s.is_u_closed() {
                let a_v_step = a_vmax - a_vmin;
                let mut a_v = a_vmin;
                while a_v <= a_vmax {
                    let p1 = a_s.point_at(a_umin, a_v);
                    let p2 = a_s.point_at(a_umax, a_v);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_u_periodic = false;
                        break;
                    }
                    a_v += a_v_step;
                }
            }
            // OCCT L236-285: 2. verify the periodicity inside the edge UV
            // box (3 or 6 point pairs).
            if is_u_periodic {
                let a_v = (a_vmin + a_vmax) * 0.5;
                let mut a_u = [0.0f64; 6];
                let mut a_upp = [0.0f64; 6];
                let mut a_nb_pnt = 0usize;
                if a_x_min < a_umin {
                    a_u[0] = a_x_min;
                    a_u[1] = (a_x_min + a_umin) * 0.5;
                    a_u[2] = a_umin;
                    a_upp[0] = a_u[0] + a_umax - a_umin;
                    a_upp[1] = a_u[1] + a_umax - a_umin;
                    a_upp[2] = a_u[2] + a_umax - a_umin;
                    a_nb_pnt += 3;
                }
                if a_x_max > a_umax {
                    a_u[a_nb_pnt] = a_umax;
                    a_u[a_nb_pnt + 1] = (a_x_max + a_umax) * 0.5;
                    a_u[a_nb_pnt + 2] = a_x_max;
                    a_upp[a_nb_pnt] = a_u[a_nb_pnt] - a_umax + a_umin;
                    a_upp[a_nb_pnt + 1] = a_u[a_nb_pnt + 1] - a_umax + a_umin;
                    a_upp[a_nb_pnt + 2] = a_u[a_nb_pnt + 2] - a_umax + a_umin;
                    a_nb_pnt += 3;
                }
                for an_ind in 0..a_nb_pnt {
                    let p1 = a_s.point_at(a_u[an_ind], a_v);
                    let p2 = a_s.point_at(a_upp[an_ind], a_v);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_u_periodic = false;
                        break;
                    }
                }
            }
        }

        // OCCT L288-299.
        if !is_u_periodic {
            if (a_x_min < a_umin) && (a_umin < a_x_max) {
                a_x_min = a_umin;
            }
            if (a_x_min < a_umax) && (a_umax < a_x_max) {
                a_x_max = a_umax;
            }
        }
    }

    // OCCT L302: the V arm — the mirror of the U arm.
    if !a_s.is_v_periodic() {
        let mut is_v_periodic = false;

        if matches!(a_s, Surface3::BSpline(_)) && (a_y_min < a_vmin || a_y_max > a_vmax) {
            let a_tol2 =
                100.0 * crate::core::precision::CONFUSION * crate::core::precision::CONFUSION;
            is_v_periodic = true;
            // OCCT L314-327: 1. verify V-closedness by sampling.
            if !a_s.is_v_closed() {
                let a_u_step = a_umax - a_umin;
                let mut a_u = a_umin;
                while a_u <= a_umax {
                    let p1 = a_s.point_at(a_u, a_vmin);
                    let p2 = a_s.point_at(a_u, a_vmax);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_v_periodic = false;
                        break;
                    }
                    a_u += a_u_step;
                }
            }
            // OCCT L328-341: 2. verify the periodicity inside the edge UV
            // box.
            if is_v_periodic {
                let a_u = (a_umin + a_umax) * 0.5;
                let mut a_v = [0.0f64; 6];
                let mut a_vpp = [0.0f64; 6];
                let mut a_nb_pnt = 0usize;
                if a_y_min < a_vmin {
                    a_v[0] = a_y_min;
                    a_v[1] = (a_y_min + a_vmin) * 0.5;
                    a_v[2] = a_vmin;
                    a_vpp[0] = a_v[0] + a_vmax - a_vmin;
                    a_vpp[1] = a_v[1] + a_vmax - a_vmin;
                    a_vpp[2] = a_v[2] + a_vmax - a_vmin;
                    a_nb_pnt += 3;
                }
                if a_y_max > a_vmax {
                    a_v[a_nb_pnt] = a_vmax;
                    a_v[a_nb_pnt + 1] = (a_y_max + a_vmax) * 0.5;
                    a_v[a_nb_pnt + 2] = a_y_max;
                    a_vpp[a_nb_pnt] = a_v[a_nb_pnt] - a_vmax + a_vmin;
                    a_vpp[a_nb_pnt + 1] = a_v[a_nb_pnt + 1] - a_vmax + a_vmin;
                    a_vpp[a_nb_pnt + 2] = a_v[a_nb_pnt + 2] - a_vmax + a_vmin;
                    a_nb_pnt += 3;
                }
                for an_ind in 0..a_nb_pnt {
                    let p1 = a_s.point_at(a_u, a_v[an_ind]);
                    let p2 = a_s.point_at(a_u, a_vpp[an_ind]);
                    if p1.distance_squared(p2) > a_tol2 {
                        is_v_periodic = false;
                        break;
                    }
                }
            }
        }

        // OCCT L344-355.
        if !is_v_periodic {
            if (a_y_min < a_vmin) && (a_vmin < a_y_max) {
                a_y_min = a_vmin;
            }
            if (a_y_min < a_vmax) && (a_vmax < a_y_max) {
                a_y_max = a_vmax;
            }
        }
    }

    // OCCT L358: aBoxS.Update(aXmin, aYmin, aXmax, aYmax).
    a_box_s.update(a_x_min, a_y_min, a_x_max, a_y_max);
    a_box_s
}

/// The OCCT `TopExp_Explorer(F, TopAbs_EDGE)` walk over the face wires.
pub fn face_edge_shapes(f: &Shape) -> Vec<Shape> {
    let mut out = Vec::new();
    if let TShape::Face(fd) = f.data.as_ref() {
        for wire in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
            if let TShape::Wire(wd) = wire.data.as_ref() {
                out.extend(wd.edges.iter().cloned());
            }
        }
    }
    out
}

// =========================================================================
// OCCT BRepAdaptor_Surface (BRepAdaptor_Surface.hxx L48-83, .cxx L113-180)
// =========================================================================

/// OCCT BRepAdaptor_Surface — a Face of the BRep topology used like a 3D
/// surface.
#[derive(Clone)]
pub struct BRepAdaptorSurface {
    /// The OCCT base class GeomAdaptor_TransformedSurface (hxx L48).
    pub transformed: TransformedSurfaceAdaptor,
    /// OCCT: TopoDS_Face myFace (hxx L82).
    pub my_face: Shape,
}

impl BRepAdaptorSurface {
    /// OCCT BRepAdaptor_Surface() (cxx L113) — the undefined adaptor with no
    /// face loaded (the rcad null shape).
    pub fn new() -> Self {
        BRepAdaptorSurface {
            transformed: TransformedSurfaceAdaptor::new(),
            my_face: Shape::null(),
        }
    }

    /// OCCT BRepAdaptor_Surface(F, R) (cxx L117-120) -> Initialize(F, R).
    /// The owning brep rides the call (the rcad location-pool read — see the
    /// module header).
    pub fn with_face(brep: &BRep, f: &Shape, r: bool) -> Self {
        let mut a = BRepAdaptorSurface::new();
        a.initialize(brep, f, r);
        a
    }

    /// OCCT ShallowCopy (cxx L124-137) — the member copy.
    pub fn shallow_copy_of(&self) -> Self {
        self.clone()
    }

    /// OCCT Initialize(F, Restriction) (cxx L141-166).
    pub fn initialize(&mut self, brep: &BRep, f: &Shape, restriction: bool) {
        // OCCT L143-146: null face -> return.
        if f.is_null() {
            return;
        }

        // OCCT L148: myFace = F.
        self.my_face = f.clone();
        // OCCT L149-150: TopLoc_Location L; aSurface = BRep_Tool::Surface(F,
        // L) — the LOCAL surface (the location rides the shape slot).
        let Some(a_surface) = brep.face_surface(f).cloned() else {
            // OCCT L151-154: null surface -> return.
            return;
        };
        let the_trsf = location_transformation(brep, f);

        if restriction {
            // OCCT L156-161: BRepTools::UVBounds(F, umin, umax, vmin, vmax);
            // Load(aSurface, umin, umax, vmin, vmax, L.Transformation()).
            let (umin, umax, vmin, vmax) = brep_tools_uv_bounds(brep, f);
            self.transformed
                .load_with_window(a_surface, umin, umax, vmin, vmax, the_trsf, 0.0, 0.0);
        } else {
            // OCCT L162-165: the natural-domain branch —
            // Load(aSurface, L.Transformation()).
            self.transformed.load(a_surface, the_trsf);
        }
    }

    /// OCCT Face() (cxx L170-173).
    pub fn face(&self) -> &Shape {
        &self.my_face
    }

    /// OCCT BRepAdaptor_Surface::Tolerance (cxx L92-95):
    /// `return BRep_Tool::Tolerance(myFace);` — i.e. the plain BRep_Tool
    /// floor, delegated to the kernel canonical reader.
    pub fn tolerance(&self) -> f64 {
        brep_tool_tolerance(self.my_face.data.as_ref())
    }
}

impl Default for BRepAdaptorSurface {
    fn default() -> Self {
        Self::new()
    }
}

// =========================================================================
// OCCT BRepAdaptor_Curve (BRepAdaptor_Curve.hxx L37-98, .cxx L132-276)
// =========================================================================

/// OCCT BRepAdaptor_Curve — an Edge of the BRep topology used like a 3D
/// curve.
#[derive(Clone)]
pub struct BRepAdaptorCurve {
    /// The OCCT base class GeomAdaptor_TransformedCurve (hxx L37).
    pub transformed: TransformedCurveAdaptor,
    /// OCCT: TopoDS_Edge myEdge (hxx L97).
    pub my_edge: Shape,
}

impl BRepAdaptorCurve {
    /// OCCT BRepAdaptor_Curve() (cxx L132).
    pub fn new() -> Self {
        BRepAdaptorCurve {
            transformed: TransformedCurveAdaptor::new(),
            my_edge: Shape::null(),
        }
    }

    /// OCCT BRepAdaptor_Curve(E) (cxx L136-139) -> Initialize(E).
    pub fn with_edge(brep: &BRep, e: &Shape) -> Self {
        let mut a = BRepAdaptorCurve::new();
        a.initialize(brep, e);
        a
    }

    /// OCCT BRepAdaptor_Curve(E, F) (cxx L143-146) -> Initialize(E, F).
    pub fn with_edge_face(brep: &BRep, e: &Shape, f: &Shape) -> Self {
        let mut a = BRepAdaptorCurve::new();
        a.initialize_on_face(brep, e, f);
        a
    }

    /// OCCT ShallowCopy (cxx L150-166) — the member copy.
    pub fn shallow_copy_of(&self) -> Self {
        self.clone()
    }

    /// OCCT Reset (cxx L170-176).
    pub fn reset(&mut self) {
        self.transformed.my_curve.reset();
        self.transformed.my_con_surf = None;
        self.my_edge = Shape::null();
        self.transformed.my_trsf = Trsf::identity();
    }

    /// OCCT Initialize(E) (cxx L180-213).
    pub fn initialize(&mut self, brep: &BRep, e: &Shape) {
        // OCCT L182: myConSurf.Nullify().
        self.transformed.my_con_surf = None;
        // OCCT L183: myEdge = E.
        self.my_edge = e.clone();
        // OCCT L184-187: double pf, pl; TopLoc_Location L; C =
        // BRep_Tool::Curve(E, L, pf, pl) — the local curve + the range.
        let c = brep.edge_curve_data(e);
        let range = brep.edge_range(e);

        if let Some(c) = c {
            // OCCT L189-192: myCurve.Load(C, pf, pl).
            self.transformed
                .my_curve
                .load_with_range(c, range[0], range[1]);
        } else {
            // OCCT L194-210: the curve-on-surface fallback — PC, S, L, pf,
            // pl = BRep_Tool::CurveOnSurface(E, PC, S, L, pf, pl).
            let Some((pc, pf, pl, s)) = first_curve_on_surface(brep, e) else {
                // OCCT L209: throw Standard_NullObject("No geometry").
                panic!("Standard_NullObject: BRepAdaptor_Curve::No geometry");
            };
            // OCCT L200-201: HS = new GeomAdaptor_Surface(); HS->Load(S).
            let mut hs = GeomSurfaceAdaptor::empty();
            hs.load(s);
            // OCCT L202-203: HC = new Geom2dAdaptor_Curve();
            // HC->Load(PC, pf, pl).
            let hc = Geom2dCurveAdaptor::with_range(pc, pf, pl);
            // OCCT L204-205: myConSurf = new Adaptor3d_CurveOnSurface();
            // myConSurf->Load(HC, HS).
            self.transformed.my_con_surf = Some(Arc::new(CurveOnSurface::new(
                Arc::new(hc) as Curve2dHandle,
                Arc::new(hs) as SurfaceHandle,
            )));
        }
        // OCCT L212: myTrsf = L.Transformation().
        self.transformed.my_trsf = location_transformation(brep, e);
    }

    /// OCCT Initialize(E, F) (cxx L217-235).
    pub fn initialize_on_face(&mut self, brep: &BRep, e: &Shape, f: &Shape) {
        // OCCT L219: myConSurf.Nullify().
        self.transformed.my_con_surf = None;
        // OCCT L221: myEdge = E.
        self.my_edge = e.clone();
        // OCCT L222-225: L; pf, pl; S = BRep_Tool::Surface(F, L);
        // PC = BRep_Tool::CurveOnSurface(E, F, pf, pl).
        let s = brep.face_surface(f).cloned();
        let pc = brep.curve_on_surface(e, f);
        if let (Some(s), Some((pc, pf, pl))) = (s, pc) {
            // OCCT L227-228: HS = new GeomAdaptor_Surface(); HS->Load(S).
            let mut hs = GeomSurfaceAdaptor::empty();
            hs.load(s);
            // OCCT L229-230: HC = new Geom2dAdaptor_Curve();
            // HC->Load(PC, pf, pl).
            let hc = Geom2dCurveAdaptor::with_range(pc, pf, pl);
            // OCCT L231-232: myConSurf = new Adaptor3d_CurveOnSurface();
            // myConSurf->Load(HC, HS).
            self.transformed.my_con_surf = Some(Arc::new(CurveOnSurface::new(
                Arc::new(hc) as Curve2dHandle,
                Arc::new(hs) as SurfaceHandle,
            )));
        }
        // OCCT L234: myTrsf = L.Transformation().
        self.transformed.my_trsf = location_transformation(brep, f);
    }

    /// OCCT Edge() (cxx L239-242).
    pub fn edge(&self) -> &Shape {
        &self.my_edge
    }

    /// OCCT BRepAdaptor_Curve::Tolerance (cxx L146-149):
    /// `return BRep_Tool::Tolerance(myEdge);` — i.e. the plain BRep_Tool
    /// floor, delegated to the kernel canonical reader.
    pub fn tolerance(&self) -> f64 {
        brep_tool_tolerance(self.my_edge.data.as_ref())
    }

    /// OCCT Trim(First, Last, Tol) (cxx L253-276): the copy-and-restore
    /// dance keeps the transformation in the copy.  The rcad encoding clones
    /// `self` and rebinds the clone's base state.
    pub fn trim_of(&self, first: f64, last: f64, tol: f64) -> Self {
        let mut res = self.clone();
        if self.transformed.my_con_surf.is_none() {
            // OCCT L261-265: the 3D curve branch — Load(C, First, Last) in
            // the copy (C = myCurve.Curve()).
            let c = res.transformed.my_curve.curve.clone();
            res.transformed.my_curve.load_with_range(c, first, last);
        } else {
            // OCCT L269-273: the curve-on-surface branch — the trimmed COS
            // in the copy (the Adaptor3d_CurveOnSurface::Trim body,
            // Adaptor3d_CurveOnSurface.cxx L1133-1141).
            let con = self.transformed.my_con_surf.as_ref().expect("con surf");
            res.transformed.my_con_surf = Some(Arc::new(CurveOnSurface::new(
                con.my2d_curve.trim(first, last, tol),
                con.my_surface.clone(),
            )));
        }
        res
    }
}

impl Default for BRepAdaptorCurve {
    fn default() -> Self {
        Self::new()
    }
}

/// The OCCT `BRep_Tool::CurveOnSurface(E, PC, S, L, pf, pl)` read (the first
/// BRep_CurveOnSurface representation, BRep_Tool.cxx L327-372): the pcurve,
/// its range and the owner-face LOCAL surface.  The rcad pcurve map keys on
/// the owner-face TShape pointer identity.
fn first_curve_on_surface(brep: &BRep, e: &Shape) -> Option<(Curve2d, f64, f64, Surface3)> {
    let ed = brep.edge_data(e)?;
    for ((face_key, _loc), (pc, pf, pl)) in ed.pcurves.iter() {
        let face_pos = brep
            .tshapes
            .iter()
            .position(|t| Arc::as_ptr(t) as u64 == *face_key);
        let Some(face_pos) = face_pos else {
            continue;
        };
        let Some(TShape::Face(fd)) = brep.tshapes.get(face_pos).map(|t| t.as_ref()) else {
            continue;
        };
        let Some(s) = fd.surface.as_ref() else {
            continue;
        };
        return Some((pc.clone(), *pf, *pl, s.clone()));
    }
    None
}

// =========================================================================
// The trait implementations — the inherited-member forwarding
// =========================================================================

macro_rules! forward_surface_trait {
    ($imp:ty) => {
        impl Adaptor3dSurface for $imp {
            fn first_u_parameter(&self) -> f64 {
                self.transformed.first_u_parameter()
            }
            fn last_u_parameter(&self) -> f64 {
                self.transformed.last_u_parameter()
            }
            fn first_v_parameter(&self) -> f64 {
                self.transformed.first_v_parameter()
            }
            fn last_v_parameter(&self) -> f64 {
                self.transformed.last_v_parameter()
            }
            fn value(&self, u: f64, v: f64) -> glam::DVec3 {
                self.transformed.value(u, v)
            }
            fn d1(&self, u: f64, v: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3) {
                self.transformed.d1(u, v)
            }
            fn d2(
                &self,
                u: f64,
                v: f64,
            ) -> (
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
            ) {
                self.transformed.d2(u, v)
            }
            #[allow(clippy::type_complexity)]
            fn d3(
                &self,
                u: f64,
                v: f64,
            ) -> (
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
                glam::DVec3,
            ) {
                self.transformed.d3(u, v)
            }
            fn dn(&self, u: f64, v: f64, nu: i32, nv: i32) -> glam::DVec3 {
                self.transformed.dn(u, v, nu, nv)
            }
            fn u_resolution(&self, r3d: f64) -> f64 {
                self.transformed.u_resolution(r3d)
            }
            fn v_resolution(&self, r3d: f64) -> f64 {
                self.transformed.v_resolution(r3d)
            }
            fn get_type(&self) -> super::adaptor::GeomAbsSurfaceType {
                self.transformed.get_type()
            }
            fn is_u_periodic(&self) -> bool {
                self.transformed.is_u_periodic()
            }
            fn u_period(&self) -> f64 {
                self.transformed.u_period()
            }
            fn is_v_periodic(&self) -> bool {
                self.transformed.is_v_periodic()
            }
            fn v_period(&self) -> f64 {
                self.transformed.v_period()
            }
            fn u_continuity(&self) -> GeomAbsShape {
                self.transformed.u_continuity()
            }
            fn v_continuity(&self) -> GeomAbsShape {
                self.transformed.v_continuity()
            }
            fn nb_u_intervals(&self, s: GeomAbsShape) -> usize {
                self.transformed.nb_u_intervals(s)
            }
            fn nb_v_intervals(&self, s: GeomAbsShape) -> usize {
                self.transformed.nb_v_intervals(s)
            }
            fn u_intervals(&self, s: GeomAbsShape) -> Vec<f64> {
                self.transformed.u_intervals(s)
            }
            fn v_intervals(&self, s: GeomAbsShape) -> Vec<f64> {
                self.transformed.v_intervals(s)
            }
            fn u_trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dSurface> {
                self.transformed.u_trim(first, last, tol)
            }
            fn v_trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dSurface> {
                self.transformed.v_trim(first, last, tol)
            }
            fn is_u_closed(&self) -> bool {
                self.transformed.is_u_closed()
            }
            fn is_v_closed(&self) -> bool {
                self.transformed.is_v_closed()
            }
            fn u_degree(&self) -> usize {
                self.transformed.u_degree()
            }
            fn nb_u_poles(&self) -> usize {
                self.transformed.nb_u_poles()
            }
            fn v_degree(&self) -> usize {
                self.transformed.v_degree()
            }
            fn nb_v_poles(&self) -> usize {
                self.transformed.nb_v_poles()
            }
            fn nb_u_knots(&self) -> usize {
                self.transformed.nb_u_knots()
            }
            fn nb_v_knots(&self) -> usize {
                self.transformed.nb_v_knots()
            }
            fn is_u_rational(&self) -> bool {
                self.transformed.is_u_rational()
            }
            fn is_v_rational(&self) -> bool {
                self.transformed.is_v_rational()
            }
            fn bezier(&self) -> crate::geom::BezierSurface {
                self.transformed.bezier()
            }
            fn bspline(&self) -> crate::geom::BSplineSurface {
                self.transformed.bspline()
            }
            fn direction(&self) -> glam::DVec3 {
                self.transformed.direction()
            }
            fn basis_curve(&self) -> Arc<dyn Adaptor3dCurve> {
                self.transformed.basis_curve()
            }
            fn basis_surface(&self) -> Arc<dyn Adaptor3dSurface> {
                self.transformed.basis_surface()
            }
            fn offset_value(&self) -> f64 {
                self.transformed.offset_value()
            }
            fn shallow_copy(&self) -> Arc<dyn Adaptor3dSurface> {
                Arc::new(self.shallow_copy_of())
            }
            fn kernel_surface(&self) -> Option<&Surface3> {
                self.transformed.kernel_surface()
            }
        }

        impl Adaptor3dSurfaceGeom for $imp {
            fn plane(&self) -> crate::geom::Plane {
                self.transformed.plane()
            }
            fn cylinder(&self) -> crate::geom::CylindricalSurface {
                self.transformed.cylinder()
            }
            fn cone(&self) -> crate::geom::ConicalSurface {
                self.transformed.cone()
            }
            fn sphere(&self) -> crate::geom::SphericalSurface {
                self.transformed.sphere()
            }
            fn torus(&self) -> crate::geom::ToroidalSurface {
                self.transformed.torus()
            }
            fn axe_of_revolution(&self) -> (glam::DVec3, glam::DVec3) {
                self.transformed.axe_of_revolution()
            }
        }
    };
}

forward_surface_trait!(BRepAdaptorSurface);

macro_rules! forward_curve_trait {
    ($imp:ty) => {
        impl Adaptor3dCurve for $imp {
            fn first_parameter(&self) -> f64 {
                self.transformed.first_parameter()
            }
            fn last_parameter(&self) -> f64 {
                self.transformed.last_parameter()
            }
            fn value(&self, u: f64) -> glam::DVec3 {
                self.transformed.value(u)
            }
            fn d1(&self, u: f64) -> (glam::DVec3, glam::DVec3) {
                self.transformed.d1(u)
            }
            fn d2(&self, u: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3) {
                self.transformed.d2(u)
            }
            fn d3(&self, u: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3, glam::DVec3) {
                self.transformed.d3(u)
            }
            fn dn(&self, u: f64, n: i32) -> glam::DVec3 {
                self.transformed.dn(u, n)
            }
            fn continuity(&self) -> GeomAbsShape {
                self.transformed.continuity()
            }
            fn nb_intervals(&self, s: GeomAbsShape) -> usize {
                self.transformed.nb_intervals(s)
            }
            fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
                self.transformed.intervals(s)
            }
            fn trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurve> {
                self.transformed.trim(first, last, tol)
            }
            fn is_closed(&self) -> bool {
                self.transformed.is_closed()
            }
            fn is_periodic(&self) -> bool {
                self.transformed.is_periodic()
            }
            fn period(&self) -> f64 {
                self.transformed.period()
            }
            fn resolution(&self, r3d: f64) -> f64 {
                self.transformed.resolution(r3d)
            }
            fn degree(&self) -> usize {
                self.transformed.degree()
            }
            fn is_rational(&self) -> bool {
                self.transformed.is_rational()
            }
            fn nb_poles(&self) -> usize {
                self.transformed.nb_poles()
            }
            fn nb_knots(&self) -> usize {
                self.transformed.nb_knots()
            }
            fn bezier(&self) -> crate::geom::BezierCurve3 {
                self.transformed.bezier()
            }
            fn bspline(&self) -> crate::geom::BSplineCurve3 {
                self.transformed.bspline()
            }
            fn offset_curve(&self) -> crate::geom::OffsetCurve3 {
                self.transformed.offset_curve()
            }
            fn shallow_copy(&self) -> Arc<dyn Adaptor3dCurve> {
                Arc::new(self.shallow_copy_of())
            }
        }

        impl Adaptor3dCurveGeom for $imp {
            fn get_type(&self) -> super::CurveType {
                self.transformed.get_type()
            }
            fn line(&self) -> crate::geom::Line3 {
                self.transformed.line()
            }
            fn circle(&self) -> crate::geom::Circle3 {
                self.transformed.circle()
            }
            fn ellipse(&self) -> crate::geom::Ellipse3 {
                self.transformed.ellipse()
            }
            fn parabola(&self) -> crate::geom::Parabola3 {
                self.transformed.parabola()
            }
            fn hyperbola(&self) -> crate::geom::Hyperbola3 {
                self.transformed.hyperbola()
            }
            fn trim_geom(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurveGeom> {
                Arc::new(self.trim_of(first, last, tol))
            }
        }
    };
}

forward_curve_trait!(BRepAdaptorCurve);
