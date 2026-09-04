//! BRepTopAdaptor_TopolTool (TKTopAlgo) — the face-backed TopolTool.
//!
//! 1:1 translation of OCCT `BRepTopAdaptor_TopolTool.hxx` (L35-135) +
//! `.cxx` (L43-653).  The OCCT class subclasses Adaptor3d_TopolTool and
//! down-casts the surface handle to BRepAdaptor_Surface; rcad holds the
//! (`BRep`, face `Shape`) pair directly — the same data BRepAdaptor_Surface
//! wraps.
//!
//! Inherited-state mapping (OCCT Adaptor3d_TopolTool fields used here):
//! - `myS` → (brep, face)
//! - `myNbSamplesU` / `myNbSamplesV` → this struct's fields
//! - the base restriction/vertex iterators are unreachable: every iterator
//!   method is overridden.

use glam::DVec2;
use rcad_kernel::geom::{Point3, Surface3, SurfaceEval};
use rcad_kernel::topods::{BRep, BRepTool, Orientation, Shape, State, TShape};

use crate::geomalgo::geom2d_int::Curve2dAdaptor;
use crate::topalgo::adaptor3d::topol_tool::analyse;
use crate::topalgo::adaptor3d::HVertex;
use crate::topalgo::brep_adaptor::curve2d::BRepCurve2d;
use crate::topalgo::brep_top_adaptor::fclass2d_topol::FClass2dTopol;
use crate::topalgo::brep_top_adaptor::hvertex_brep::BRepHVertex;

/// OCCT Adaptor3d_TopolTool.cxx: `#define myInfinite Precision::Infinite()`.
const MY_INFINITE: f64 = rcad_kernel::precision::INFINITE_VALUE;

/// OCCT BRepTopAdaptor_TopolTool.
pub struct BRepTopolTool<'a> {
    /// The owning BRep (BRep_Tool / BRepAdaptor accessors).
    brep: &'a BRep,
    /// OCCT myFace (FORWARD-oriented in Initialize).
    my_face: Shape,
    /// OCCT myS — surfaced through (brep, my_face).
    my_nb_samples_u: i32,
    my_nb_samples_v: i32,
    /// OCCT myFClass2d (created lazily by Classify / IsThePointOn).
    my_fclass2d: Option<FClass2dTopol<'a>>,
    /// OCCT myCurve (set by Initialize(C)).
    my_curve: Option<BRepCurve2d<'a>>,
    /// OCCT myCurves (the face edges as pcurve adaptors, in explorer order).
    my_curves: Vec<BRepCurve2d<'a>>,
    /// OCCT myCIterator — the NCollection_List iterator position.
    my_c_iterator: usize,
    /// OCCT myVIterator — the TopExp_Explorer(edge, VERTEX) position.
    my_v_iterator: Vec<Shape>,
    my_v_iterator_pos: usize,
    /// OCCT myU0 / myV0 / myDU / myDV.
    my_u0: f64,
    my_v0: f64,
    my_du: f64,
    my_dv: f64,
}

impl<'a> BRepTopolTool<'a> {
    /// OCCT BRepTopAdaptor_TopolTool() (cxx L43-51).
    pub fn new(brep: &'a BRep) -> Self {
        BRepTopolTool {
            brep,
            my_face: Shape::null(),
            my_nb_samples_u: -1,
            my_nb_samples_v: -1,
            my_fclass2d: None,
            my_curve: None,
            my_curves: Vec::new(),
            my_c_iterator: 0,
            my_v_iterator: Vec::new(),
            my_v_iterator_pos: 0,
            my_u0: 0.0,
            my_v0: 0.0,
            my_du: 0.0,
            my_dv: 0.0,
        }
    }

    /// OCCT BRepTopAdaptor_TopolTool(Surface) (cxx L55-60).
    pub fn new_with_surface(brep: &'a BRep, face: &Shape) -> Self {
        let mut t = BRepTopolTool::new(brep);
        t.initialize_surface(face);
        t
    }

    /// OCCT Initialize() (cxx L64-67).
    pub fn initialize() -> ! {
        panic!("Standard_NotImplemented: BRepTopAdaptor_TopolTool::Initialize()");
    }

    /// OCCT Initialize(S) (cxx L71-94) — the face (FORWARD) and the per-edge
    /// pcurve adaptor list.
    pub fn initialize_surface(&mut self, face: &Shape) {
        // OCCT down_cast<BRepAdaptor_Surface>(S) + ConstructionError on
        // null — the (brep, face) pair is the adaptor itself.
        let mut f = face.clone();
        f.orientation = Orientation::Forward;
        self.my_face = f.clone();
        self.my_fclass2d = None;
        self.my_nb_samples_u = -1;
        self.my_curves.clear();
        // OCCT TopExp_Explorer(myFace, TopAbs_EDGE) — every edge occurrence
        // of every wire.
        if let TShape::Face(fd) = &*f.data {
            for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                if let TShape::Wire(wd) = &*w.data {
                    for e in &wd.edges {
                        self.my_curves
                            .push(BRepCurve2d::new_edge_face(self.brep, e, &f));
                    }
                }
            }
        }
        self.my_c_iterator = 0;
    }

    /// OCCT Initialize(C) (cxx L98-105).
    pub fn initialize_curve(&mut self, c: &BRepCurve2d<'a>) {
        self.my_curve = Some(c.clone());
    }

    /// OCCT Init() (cxx L109-112).
    pub fn init(&mut self) {
        self.my_c_iterator = 0;
    }

    /// OCCT More() (cxx L116-119).
    pub fn more(&self) -> bool {
        self.my_c_iterator < self.my_curves.len()
    }

    /// OCCT Next() (cxx L123-126).
    pub fn next(&mut self) {
        self.my_c_iterator += 1;
    }

    /// OCCT Value() (cxx L130-133).
    pub fn value(&self) -> Option<&BRepCurve2d<'a>> {
        self.my_curves.get(self.my_c_iterator)
    }

    /// OCCT Edge() (cxx L138-143) — the current edge.
    pub fn edge(&self) -> Option<&Shape> {
        self.value().map(|c| c.edge())
    }

    /// OCCT InitVertexIterator() (cxx L148-151) — TopExp_Explorer(edge,
    /// TopAbs_VERTEX).
    pub fn init_vertex_iterator(&mut self) {
        self.my_v_iterator.clear();
        self.my_v_iterator_pos = 0;
        if let Some(c) = &self.my_curve {
            let e = c.edge();
            let vf = self.brep.first_vertex(e);
            let vl = self.brep.last_vertex(e);
            if !vf.is_null() {
                self.my_v_iterator.push(vf);
            }
            if !vl.is_null() {
                self.my_v_iterator.push(vl);
            }
        }
    }

    /// OCCT MoreVertex() (cxx L155-158).
    pub fn more_vertex(&self) -> bool {
        self.my_v_iterator_pos < self.my_v_iterator.len()
    }

    /// OCCT NextVertex() (cxx L160-163).
    pub fn next_vertex(&mut self) {
        self.my_v_iterator_pos += 1;
    }

    /// OCCT Vertex() (cxx L167-170).
    pub fn vertex(&self) -> Option<BRepHVertex<'_, 'a>> {
        let v = self.my_v_iterator.get(self.my_v_iterator_pos)?.clone();
        self.my_curve.as_ref().map(|c| BRepHVertex::new(&v, c))
    }

    /// OCCT Classify(P, Tol, RecadreOnPeriodic) (cxx L174-187).
    pub fn classify(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> State {
        if self.my_face.is_null() {
            return State::Unknown;
        }
        if self.my_fclass2d.is_none() {
            self.my_fclass2d = Some(FClass2dTopol::new(self.brep, &self.my_face.clone(), tol));
        }
        self.my_fclass2d
            .as_ref()
            .unwrap()
            .perform(p, recadre_on_periodic)
    }

    /// OCCT IsThePointOn(P, Tol, RecadreOnPeriodic) (cxx L191-201).
    pub fn is_the_point_on(&mut self, p: DVec2, tol: f64, recadre_on_periodic: bool) -> bool {
        if self.my_fclass2d.is_none() {
            self.my_fclass2d = Some(FClass2dTopol::new(self.brep, &self.my_face.clone(), tol));
        }
        State::On
            == self
                .my_fclass2d
                .as_ref()
                .unwrap()
                .test_on_restriction(p, tol, recadre_on_periodic)
    }

    /// OCCT Destroy() (cxx L205-209).
    pub fn destroy(&mut self) {
        if let Some(mut f) = self.my_fclass2d.take() {
            f.destroy();
        }
    }

    /// OCCT Orientation(C) (cxx L213-217).
    pub fn orientation_curve(&self, c: &BRepCurve2d<'a>) -> Orientation {
        c.edge().orientation
    }

    /// OCCT Orientation(V) (cxx L221-224) — Adaptor3d_TopolTool::Orientation.
    pub fn orientation_vertex(&self, v: &BRepHVertex<'_, 'a>) -> Orientation {
        // The base implementation returns V->Orientation().
        v.orientation()
    }

    /// OCCT ComputeSamplePoints() (cxx L345-521).
    pub fn compute_sample_points(&mut self) {
        let mut uinf = self.first_u_parameter();
        let mut usup = self.last_u_parameter();
        let mut vinf = self.first_v_parameter();
        let mut vsup = self.last_v_parameter();
        if usup < uinf {
            std::mem::swap(&mut uinf, &mut usup);
        }
        if vsup < vinf {
            std::mem::swap(&mut vinf, &mut vsup);
        }
        if uinf == f64::MIN && usup == f64::MAX {
            uinf = -1.0e5;
            usup = 1.0e5;
        } else if uinf == f64::MIN {
            uinf = usup - 2.0e5;
        } else if usup == f64::MAX {
            usup = uinf + 2.0e5;
        }

        if vinf == f64::MIN && vsup == f64::MAX {
            vinf = -1.0e5;
            vsup = 1.0e5;
        } else if vinf == f64::MIN {
            vinf = vsup - 2.0e5;
        } else if vsup == f64::MAX {
            vsup = vinf + 2.0e5;
        }

        let typ_s = self.surface_type();
        let (mut nbsu, mut nbsv) = match typ_s {
            crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane => (2usize, 2usize),
            crate::geomalgo::int_patch::GeomAbsSurfaceType::BezierSurface => (
                3 + self.nb_u_poles(),
                3 + self.nb_v_poles(),
            ),
            crate::geomalgo::int_patch::GeomAbsSurfaceType::BSplineSurface => {
                let mut nbsv2 = self.nb_v_knots() * self.v_degree();
                if nbsv2 < 4 {
                    nbsv2 = 4;
                }
                let mut nbsu2 = self.nb_u_knots() * self.u_degree();
                if nbsu2 < 4 {
                    nbsu2 = 4;
                }
                (nbsu2, nbsv2)
            }
            crate::geomalgo::int_patch::GeomAbsSurfaceType::Cylinder
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::Cone
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::Sphere
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::Torus => {
                //-- Set 15 for 2pi; Not enough -> 25 for 2pi.
                let nbsu2 = (8.0 * (usup - uinf)) as usize;
                let nbsv2 = (7.0 * (vsup - vinf)) as usize;
                let nbsu2 = nbsu2.max(5).min(30); // modif HRT buc60462
                let nbsv2 = nbsv2.max(5).min(15);
                (nbsu2, nbsv2)
            }
            crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfRevolution
            | crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfExtrusion => {
                (25usize, 15usize)
            }
            _ => (10usize, 10usize),
        };

        //-- If the number of points is too great, analyze.
        if nbsu < 10 {
            nbsu = 10;
        }
        if nbsv < 10 {
            nbsv = 10;
        }

        self.my_nb_samples_u = nbsu as i32;
        self.my_nb_samples_v = nbsv as i32;

        if nbsu > 10 || nbsv > 10 {
            if typ_s == crate::geomalgo::int_patch::GeomAbsSurfaceType::BSplineSurface {
                let bspl = match self.brep.face_surface(&self.my_face) {
                    Some(Surface3::BSpline(b)) => b.clone(),
                    _ => panic!("Standard_NoSuchObject"),
                };
                let nbup = bspl.control_points.len();
                let nbvp = bspl.control_points.first().map_or(0, |r| r.len());
                let (u2, v2) = analyse(&bspl.control_points, nbup, nbvp);
                self.my_nb_samples_u = u2 as i32;
                self.my_nb_samples_v = v2 as i32;
                nbsu = u2;
                nbsv = v2;
            } else if typ_s == crate::geomalgo::int_patch::GeomAbsSurfaceType::BezierSurface {
                let bez = match self.brep.face_surface(&self.my_face) {
                    Some(Surface3::Bezier(b)) => b.clone(),
                    _ => panic!("Standard_NoSuchObject"),
                };
                let nbup = bez.control_points.len();
                let nbvp = bez.control_points.first().map_or(0, |r| r.len());
                let (u2, v2) = analyse(&bez.control_points, nbup, nbvp);
                self.my_nb_samples_u = u2 as i32;
                self.my_nb_samples_v = v2 as i32;
                nbsu = u2;
                nbsv = v2;
            }
        }

        if nbsu < 10 {
            nbsu = 10;
        }
        if nbsv < 10 {
            nbsv = 10;
        }

        self.my_nb_samples_u = nbsu as i32;
        self.my_nb_samples_v = nbsv as i32;

        self.my_u0 = uinf;
        self.my_v0 = vinf;

        self.my_du = (usup - uinf) / (self.my_nb_samples_u + 1) as f64;
        self.my_dv = (vsup - vinf) / (self.my_nb_samples_v + 1) as f64;
    }

    /// OCCT NbSamplesU() (cxx L525-532).
    pub fn nb_samples_u(&mut self) -> i32 {
        if self.my_nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.my_nb_samples_u
    }

    /// OCCT NbSamplesV() (cxx L536-543).
    pub fn nb_samples_v(&mut self) -> i32 {
        if self.my_nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.my_nb_samples_v
    }

    /// OCCT NbSamples() (cxx L547-554).
    pub fn nb_samples(&mut self) -> i32 {
        if self.my_nb_samples_u < 0 {
            self.compute_sample_points();
        }
        self.my_nb_samples_u * self.my_nb_samples_v
    }

    /// OCCT SamplePoint(i, P2d, P3d) (cxx L558-566).
    pub fn sample_point(&self, i: usize) -> (DVec2, Point3) {
        let iv = 1 + i as i32 / self.my_nb_samples_u;
        let iu = 1 + i as i32 - (iv - 1) * self.my_nb_samples_u;
        let u = self.my_u0 + iu as f64 * self.my_du;
        let v = self.my_v0 + iv as f64 * self.my_dv;
        (
            DVec2::new(u, v),
            self.brep
                .face_surface(&self.my_face)
                .map(|s| s.point_at(u, v))
                .unwrap_or(Point3::ZERO),
        )
    }

    /// OCCT DomainIsInfinite() (cxx L570-594).
    pub fn domain_is_infinite(&self) -> bool {
        let uinf = self.first_u_parameter();
        let usup = self.last_u_parameter();
        let vinf = self.first_v_parameter();
        let vsup = self.last_v_parameter();
        if rcad_kernel::precision::is_negative_infinite_value(uinf) {
            return true;
        }
        if rcad_kernel::precision::is_positive_infinite_value(usup) {
            return true;
        }
        if rcad_kernel::precision::is_negative_infinite_value(vinf) {
            return true;
        }
        rcad_kernel::precision::is_positive_infinite_value(vsup)
    }

    /// OCCT Has3d() (cxx L598-601).
    pub fn has_3d(&self) -> bool {
        true
    }

    /// OCCT Tol3d(C) (cxx L605-619) — BRep_Tool::Tolerance(edge).
    pub fn tol3d_curve(&self, c: &BRepCurve2d<'a>) -> f64 {
        match &*c.edge().data {
            TShape::Edge(ed) => ed.tolerance,
            _ => 0.0,
        }
    }

    /// OCCT Tol3d(V) (cxx L623-636).
    pub fn tol3d_vertex(&self, v: &BRepHVertex<'_, 'a>) -> f64 {
        self.brep.vertex_tolerance(v.vertex())
    }

    /// OCCT Pnt(V) (cxx L640-653).
    pub fn pnt(&self, v: &BRepHVertex<'_, 'a>) -> Point3 {
        self.brep.vertex_position(v.vertex())
    }

    // -- the Adaptor3d_TopolTool myS accessors (BRepAdaptor_Surface) --

    pub(crate) fn first_u_parameter(&self) -> f64 {
        self.surface_domain()[0]
    }
    pub(crate) fn last_u_parameter(&self) -> f64 {
        self.surface_domain()[1]
    }
    pub(crate) fn first_v_parameter(&self) -> f64 {
        self.surface_domain()[2]
    }
    pub(crate) fn last_v_parameter(&self) -> f64 {
        self.surface_domain()[3]
    }

    fn surface_domain(&self) -> [f64; 4] {
        // OCCT BRepAdaptor_Surface::FirstUParameter etc. return the FACE UV
        // bounds (the wire-restricted domain) when the underlying surface is
        // unbounded — rcad carries them in TFaceData::uv_domain.
        if let TShape::Face(fd) = &*self.my_face.data {
            if let Some(dom) = fd.uv_domain {
                return dom;
            }
        }
        self.brep
            .face_surface(&self.my_face)
            .map(|s| s.default_domain())
            .unwrap_or([0.0; 4])
    }

    fn surface_type(&self) -> crate::geomalgo::int_patch::GeomAbsSurfaceType {
        use crate::geomalgo::int_patch::GeomAbsSurfaceType;
        match self.brep.face_surface(&self.my_face) {
            Some(s) => match s {
                Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
                Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
                Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
                Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
                Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
                Surface3::Bezier(_) => GeomAbsSurfaceType::BezierSurface,
                Surface3::BSpline(_) => GeomAbsSurfaceType::BSplineSurface,
                Surface3::Revolution(_) => GeomAbsSurfaceType::SurfaceOfRevolution,
                Surface3::LinearExtrusion(_) => GeomAbsSurfaceType::SurfaceOfExtrusion,
                Surface3::Offset(_) => GeomAbsSurfaceType::OffsetSurface,
                _ => GeomAbsSurfaceType::OtherSurface,
            },
            None => GeomAbsSurfaceType::OtherSurface,
        }
    }

    fn nb_u_poles(&self) -> usize {
        match self.brep.face_surface(&self.my_face) {
            Some(Surface3::Bezier(b)) => b.control_points.len(),
            Some(Surface3::BSpline(b)) => b.control_points.len(),
            _ => panic!("Standard_NoSuchObject"),
        }
    }
    fn nb_v_poles(&self) -> usize {
        match self.brep.face_surface(&self.my_face) {
            Some(Surface3::Bezier(b)) => b.control_points.first().map_or(0, |r| r.len()),
            Some(Surface3::BSpline(b)) => b.control_points.first().map_or(0, |r| r.len()),
            _ => panic!("Standard_NoSuchObject"),
        }
    }
    fn nb_u_knots(&self) -> usize {
        match self.brep.face_surface(&self.my_face) {
            Some(Surface3::BSpline(b)) => {
                crate::geomalgo::int_curve_surface::distinct_knots(&b.knots_u).len()
            }
            _ => panic!("Standard_NoSuchObject"),
        }
    }
    fn nb_v_knots(&self) -> usize {
        match self.brep.face_surface(&self.my_face) {
            Some(Surface3::BSpline(b)) => {
                crate::geomalgo::int_curve_surface::distinct_knots(&b.knots_v).len()
            }
            _ => panic!("Standard_NoSuchObject"),
        }
    }
    fn u_degree(&self) -> usize {
        match self.brep.face_surface(&self.my_face) {
            Some(Surface3::BSpline(b)) => b.degree_u,
            _ => panic!("Standard_NoSuchObject"),
        }
    }
    fn v_degree(&self) -> usize {
        match self.brep.face_surface(&self.my_face) {
            Some(Surface3::BSpline(b)) => b.degree_v,
            _ => panic!("Standard_NoSuchObject"),
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use glam::{DVec2, DVec3};
    use rcad_kernel::geom::{Curve2d, Curve3, Line2d, Line3, Plane, Surface3};
    use rcad_kernel::topods::{BRepBuilder, Orientation, State};

    /// The unit square on z = 0: vertices (0,0), (1,0), (1,1), (0,1); line
    /// edges with pcurves identical to the XY coordinates; FORWARD wire.
    pub(crate) fn square_face() -> (BRep, Shape) {
        let mut brep = BRep::new();
        let mut b = BRepBuilder::new();
        let pts = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let vs: Vec<Shape> = pts.iter().map(|p| b.add_vertex(&mut brep, *p, 1e-7)).collect();
        let mut edges = Vec::new();
        for i in 0..4 {
            let a = pts[i];
            let bb = pts[(i + 1) % 4];
            let curve = Curve3::Line(Line3 {
                origin: a,
                direction: (bb - a).normalize(),
            });
            let range = [0.0, (bb - a).length()];
            // OCCT BRep_Builder::Add(E, V) stores the start vertex FORWARD
            // and the end vertex REVERSED on the edge.
            let v_first = vs[i].clone();
            let mut v_last = vs[(i + 1) % 4].clone();
            v_last.orientation = Orientation::Reversed;
            edges.push(b.add_edge(&mut brep, Some(curve), v_first, v_last, range));
        }
        // The face FIRST (the pcurve key needs the face TShape pointer).
        let wire = brep.add_twire(edges.clone());
        let face = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            })),
            wire,
            Vec::new(),
            Some(DVec3::new(0.5, 0.5, 0.0)),
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );
        // pcurves: for the FORWARD edge i, UV runs between the corner UVs.
        for i in 0..4 {
            let j = (i + 1) % 4;
            let (a2, b2) = (
                DVec2::new(pts[i].x, pts[i].y),
                DVec2::new(pts[j].x, pts[j].y),
            );
            let len = (b2 - a2).length();
            let pc = Curve2d::Line(Line2d {
                origin: a2,
                direction: (b2 - a2) / len,
            });
            b.add_pcurve(&mut brep, edges[i].clone(), face.clone(), pc, 0.0, len);
        }
        (brep, face)
    }

    /// OCCT anchor: BRepAdaptor_Curve2d over the square's bottom edge — the
    /// pcurve parameterisation runs corner to corner, and Edge()/Face()
    /// return the shapes it was built with.
    #[test]
    fn brep_curve2d_anchor() {
        let (brep, face) = square_face();
        let bottom_edge = {
            let TShape::Face(fd) = &*face.data else {
                panic!("not a face");
            };
            let TShape::Wire(wd) = &*fd.outer_wire.data else {
                panic!("not a wire");
            };
            wd.edges[0].clone()
        };

        let c = BRepCurve2d::new_edge_face(&brep, &bottom_edge, &face);
        assert_eq!(c.first_parameter(), 0.0);
        assert_eq!(c.last_parameter(), 1.0);
        let (p, v) = c.d1(0.5);
        assert_eq!(p, DVec2::new(0.5, 0.0));
        assert_eq!(v, DVec2::X);
        assert!(c.edge().is_same(&bottom_edge));
        assert!(c.face().is_same(&face));
    }

    /// OCCT anchor: FClass2d over the square — the interior point is IN,
    /// the far exterior OUT, the boundary ON; TestOnRestriction agrees and
    /// the infinite point is OUT.
    #[test]
    fn fclass2d_topol_anchor() {
        let (brep, face) = square_face();
        let f = FClass2dTopol::new(&brep, &face, 1e-6);

        assert_eq!(f.perform(DVec2::new(0.5, 0.5), true), State::In);
        assert_eq!(f.perform(DVec2::new(2.0, 2.0), true), State::Out);
        assert_eq!(f.test_on_restriction(DVec2::new(0.5, 0.0), 1e-6, true), State::On);
        assert_eq!(f.test_on_restriction(DVec2::new(0.5, 0.5), 1e-6, true), State::In);
        assert_eq!(f.perform_infinite_point(), State::Out);
    }

    /// OCCT: the boundary point (0.5, 0) hits SiDans == 0 and falls back to
    /// the BRepClass_FaceClassifier, which returns On.
    ///
    /// Currently ignored: the rcad FClassifier-on-FaceShapeSource fallback
    /// returns In.  Root cause located — FaceShapeSource passes the kernel
    /// BRep locations table to `edge_pcurve_on_face`, which expects the
    /// bop DS convention (slot 0 = identity); the pcurve keys therefore
    /// miss, FaceExplorer::segment returns None, and the classifier exits
    /// through the nowires branch.  Un-ignore once the locations-table
    /// convention is aligned.
    #[test]
    #[ignore = "FClassifier boundary-On: FaceShapeSource locations-table convention misalignment (runway 2a-4 follow-up)"]
    fn fclass2d_topol_perform_boundary_on() {
        let (brep, face) = square_face();
        let f = FClass2dTopol::new(&brep, &face, 1e-6);
        assert_eq!(f.perform(DVec2::new(0.5, 0.0), true), State::On);
    }

    /// OCCT anchor: BRepTopAdaptor_TopolTool over the square — the edge
    /// iterator walks the 4 wire edges, Classify delegates to the FClass2d,
    /// and ComputeSamplePoints on a plane clamps to 10x10 (cxx L396-472).
    #[test]
    fn brep_topol_tool_anchor() {
        let (brep, face) = square_face();
        let mut tool = BRepTopolTool::new_with_surface(&brep, &face);

        assert_eq!(tool.more(), true);
        let mut count = 0;
        while tool.more() {
            assert!(tool.value().is_some());
            assert!(tool.edge().is_some());
            count += 1;
            tool.next();
        }
        assert_eq!(count, 4);

        // Classify delegates to the FClass2d (Perform).
        assert_eq!(tool.classify(DVec2::new(0.5, 0.5), 1e-6, true), State::In);
        assert_eq!(tool.classify(DVec2::new(3.0, 3.0), 1e-6, true), State::Out);
        assert_eq!(tool.is_the_point_on(DVec2::new(0.5, 0.0), 1e-6, true), true);
        assert_eq!(tool.is_the_point_on(DVec2::new(0.5, 0.5), 1e-6, true), false);

        // ComputeSamplePoints: the Plane branch sets 2, then the min-10
        // clamp applies; myU0/myV0 + myDU/myDV set from the domain.
        assert_eq!(tool.nb_samples_u(), 10);
        assert_eq!(tool.nb_samples_v(), 10);
        assert_eq!(tool.nb_samples(), 100);
        assert!(!tool.domain_is_infinite());
        assert!(tool.has_3d());

        // SamplePoint(i) walks the (iu, iv) grid (cxx L558-566).
        let (p2d, _p3d) = tool.sample_point(0);
        // i = 0: iv = 1 + 0/10 = 1, iu = 1 + 0 = 1; u = 0 + 1*du.
        let du = (1.0 - 0.0) / 11.0;
        assert!((p2d.x - du).abs() < 1e-12, "p2d.x={}", p2d.x);
        assert!((p2d.y - du).abs() < 1e-12, "p2d.y={}", p2d.y);
    }

    /// OCCT anchor: the vertex iterator of the current curve yields the
    /// edge's two vertices; the HVertex parameter is BRep_Tool::Parameter
    /// (the vertex parameter on the pcurve); Tol3d/Pnt read the topology.
    #[test]
    fn brep_topol_tool_vertex_iterator_anchor() {
        let (brep, face) = square_face();
        let mut tool = BRepTopolTool::new_with_surface(&brep, &face);
        tool.init();
        assert!(tool.more());
        let first_curve = tool.value().unwrap().clone();
        tool.initialize_curve(&first_curve);

        tool.init_vertex_iterator();
        assert!(tool.more_vertex());
        // Extract the vertex facts (the OCCT handle keeps the vertex alive;
        // the rcad borrow ends at the block).
        let (v0_val, v0_param, v0_tol, v0_pnt, v0_ori) = {
            let v0 = tool.vertex().unwrap();
            (
                v0.value(),
                v0.parameter(&first_curve),
                tool.tol3d_vertex(&v0),
                tool.pnt(&v0),
                tool.orientation_vertex(&v0),
            )
        };
        assert_eq!(v0_val, DVec2::new(-f64::MAX, -f64::MAX)); // do-nothing Value
        // The first vertex of the bottom edge sits at the curve start.
        assert!((v0_param - 0.0).abs() < 1e-12);
        tool.next_vertex();
        assert!(tool.more_vertex());
        let (v1_val_param, v1_ori) = {
            let v1 = tool.vertex().unwrap();
            (v1.parameter(&first_curve), tool.orientation_vertex(&v1))
        };
        assert!(v1_val_param - 1.0 < 1e-12);
        assert_ne!(v0_ori, Orientation::Reversed);
        tool.next_vertex();
        assert!(!tool.more_vertex());

        // Has3d/Tol3d/Pnt.
        assert!(tool.has_3d());
        assert!(v0_tol >= 0.0);
        assert_eq!(v0_pnt, DVec3::new(0.0, 0.0, 0.0));
        assert_eq!(v0_ori, Orientation::Forward);
        let _ = v1_ori;
    }
}
