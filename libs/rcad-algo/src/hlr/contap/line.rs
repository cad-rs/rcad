// OCCT Contap_Line (TKHLR) — a contour line: a gp_Lin, a gp_Circ, a
// walking polyline (IntSurf_LineOn2S) or a restriction arc.
//
// Contap_Line.hxx L38-99 + Contap_Line.cxx L25-120 + Contap_Line.lxx
// L28-87.
//
// Ownership notes: OCCT shares the IntSurf_LineOn2S handle between the
// walking line and the Contap_Line (SetLineOn2S).  In the Contap_Contour
// flow the IWLine is never read after SetLineOn2S and every mutation goes
// through the Contap_Line, so an owned copy is semantically equivalent;
// the same holds for the vertex sequence (ChangeSequence mutation through
// the sequence maps to &mut access in Rust).

use rcad_kernel::geom::{Circle3, Line3, Point3, Vec3};

use crate::geomalgo::int_patch::transitions::TypeTrans;
use crate::geomalgo::int_surf::{LineOn2S, PntOn2S};

use super::i_type::IType;
use super::point::{Arc, Point};

/// OCCT Contap_Line.
#[derive(Clone)]
pub struct Line {
    trans: TypeTrans,
    curv: Option<LineOn2S>,
    svtx: Vec<Point>,
    thearc: Option<Arc>,
    typ_l: IType,
    pt: Point3,
    dir1: Vec3,
    dir2: Vec3,
    rad: f64,
}

impl Line {
    /// OCCT Contap_Line() (cxx L25-29).
    pub fn new() -> Self {
        Line {
            trans: TypeTrans::Undecided,
            curv: None,
            svtx: Vec::new(),
            thearc: None,
            typ_l: IType::Walking,
            pt: Point3::ZERO,
            dir1: Vec3::ZERO,
            dir2: Vec3::ZERO,
            rad: 0.0,
        }
    }

    /// OCCT LineOn2S (lxx L28-31).
    pub fn line_on_2s(&self) -> Option<&LineOn2S> {
        self.curv.as_ref()
    }

    /// OCCT LineOn2S — mutable view (the handle content is mutated in
    /// ComputeInternalPoints).
    pub fn line_on_2s_mut(&mut self) -> &mut LineOn2S {
        self.curv.as_mut().expect("Contap_Line::LineOn2S")
    }

    /// OCCT ResetSeqOfVertex (cxx L31-34).
    pub fn reset_seq_of_vertex(&mut self) {
        self.svtx = Vec::new();
    }

    /// OCCT Add(const IntSurf_PntOn2S&) (lxx L33-36).
    pub fn add_pnt_on_2s(&mut self, p: &PntOn2S) {
        self.curv.as_mut().expect("Contap_Line::Add").add(p);
    }

    /// OCCT Add(const Contap_Point&) (cxx L36-63) — vertices are kept
    /// sorted by ParameterOnLine.
    pub fn add(&mut self, p: Point) {
        let n = self.svtx.len();
        if n == 0 {
            self.svtx.push(p);
        } else {
            let prm = p.parameter_on_line();
            if prm > self.svtx[n - 1].parameter_on_line() {
                self.svtx.push(p);
            } else {
                for i in (1..n).rev() {
                    if prm > self.svtx[i - 1].parameter_on_line() {
                        self.svtx.insert(i, p);
                        return;
                    }
                }
                self.svtx.insert(0, p);
            }
        }
    }

    /// OCCT Clear (cxx L65-73).
    pub fn clear(&mut self) {
        if let Some(c) = &mut self.curv {
            c.clear();
        }
        self.svtx = Vec::new();
        self.typ_l = IType::Walking;
    }

    /// OCCT SetValue(const gp_Lin&) (cxx L75-80).
    pub fn set_value_line(&mut self, l: &Line3) {
        self.pt = l.origin;
        self.dir1 = l.direction;
        self.typ_l = IType::Lin;
    }

    /// OCCT SetValue(const gp_Circ&) (cxx L82-89).
    pub fn set_value_circle(&mut self, c: &Circle3) {
        self.pt = c.center;
        self.dir1 = c.normal;
        self.dir2 = c.x_dir;
        self.rad = c.radius;
        self.typ_l = IType::Circle;
    }

    /// OCCT SetValue(const handle<Adaptor2d_Curve2d>&) (cxx L91-95).
    pub fn set_value_arc(&mut self, a: Arc) {
        self.thearc = Some(a);
        self.typ_l = IType::Restriction;
    }

    /// OCCT SetLineOn2S (cxx L97-101).
    pub fn set_line_on_2s(&mut self, l: &LineOn2S) {
        self.curv = Some(l.clone());
        self.typ_l = IType::Walking;
    }

    /// OCCT SetTransitionOnS (cxx L103-106).
    pub fn set_transition_on_s(&mut self, t: TypeTrans) {
        self.trans = t;
    }

    /// OCCT TransitionOnS (cxx L108-111).
    pub fn transition_on_s(&self) -> TypeTrans {
        self.trans
    }

    /// OCCT Arc (cxx L113-120) — raises Standard_DomainError when the line
    /// is not a restriction.
    pub fn arc(&self) -> Arc {
        if self.typ_l != IType::Restriction {
            panic!("Standard_DomainError: Contap_Line::Arc");
        }
        self.thearc.clone().unwrap()
    }

    /// OCCT NbVertex (lxx L38-41).
    pub fn nb_vertex(&self) -> usize {
        self.svtx.len()
    }

    /// OCCT Vertex(Index) (lxx L43-46) — the 1-based mutable
    /// ChangeSequence view.
    pub fn vertex_mut(&mut self, index: usize) -> &mut Point {
        &mut self.svtx[index - 1]
    }

    /// OCCT Vertex(Index) — read-only view.
    pub fn vertex(&self, index: usize) -> &Point {
        &self.svtx[index - 1]
    }

    /// OCCT TypeContour (lxx L48-51).
    pub fn type_contour(&self) -> IType {
        self.typ_l
    }

    /// OCCT NbPnts (lxx L53-60).
    pub fn nb_pnts(&self) -> usize {
        if self.typ_l != IType::Walking {
            panic!("Standard_DomainError: Contap_Line::NbPnts");
        }
        self.curv.as_ref().expect("Contap_Line::NbPnts").nb_points()
    }

    /// OCCT Point(Index) (lxx L62-69) — 1-based (myLine->Value(Index) with
    /// the 1-based IntSurf_LineOn2S::Value).
    pub fn point(&self, index: usize) -> &PntOn2S {
        if self.typ_l != IType::Walking {
            panic!("Standard_DomainError: Contap_Line::Point");
        }
        self.curv.as_ref().unwrap().value(index - 1)
    }

    /// OCCT Line (lxx L71-78).
    pub fn line(&self) -> Line3 {
        if self.typ_l != IType::Lin {
            panic!("Standard_DomainError: Contap_Line::Line");
        }
        Line3::new(self.pt, self.dir1)
    }

    /// OCCT Circle (lxx L80-87) — gp_Circ(gp_Ax2(pt, dir1, dir2), rad):
    /// the Y direction of the axis is dir1 ^ dir2.
    pub fn circle(&self) -> Circle3 {
        if self.typ_l != IType::Circle {
            panic!("Standard_DomainError: Contap_Line::Circle");
        }
        Circle3 {
            center: self.pt,
            normal: self.dir1,
            x_dir: self.dir2,
            y_dir: self.dir1.cross(self.dir2),
            radius: self.rad,
        }
    }
}

impl Default for Line {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::Circle3;

    /// OCCT anchor: Add keeps the vertex sequence sorted by the parameter
    /// on the line (cxx L36-63); NbVertex/Vertex round-trip.
    #[test]
    fn contap_line_vertices_sorted() {
        let mut l = Line::new();
        let mut a = Point::with_uv(glam::DVec3::ONE, 0.0, 0.0);
        a.set_parameter(3.0);
        let mut b = Point::with_uv(glam::DVec3::ZERO, 0.0, 0.0);
        b.set_parameter(1.0);
        let mut c = Point::with_uv(-glam::DVec3::ONE, 0.0, 0.0);
        c.set_parameter(2.0);
        l.add(a);
        l.add(b);
        l.add(c);
        assert_eq!(l.nb_vertex(), 3);
        assert!((l.vertex(1).parameter_on_line() - 1.0).abs() < 1e-15);
        assert!((l.vertex(2).parameter_on_line() - 2.0).abs() < 1e-15);
        assert!((l.vertex(3).parameter_on_line() - 3.0).abs() < 1e-15);
    }

    /// OCCT anchor: SetValue(gp_Circ) stores location/axis/xdir/radius and
    /// Circle() reconstructs the same circle (lxx L80-87).
    #[test]
    fn contap_line_circle_roundtrip() {
        let circ = Circle3 {
            center: glam::DVec3::new(1.0, 0.0, 0.0),
            normal: glam::DVec3::new(0.0, 0.0, 1.0),
            x_dir: glam::DVec3::new(1.0, 0.0, 0.0),
            y_dir: glam::DVec3::new(0.0, 1.0, 0.0),
            radius: 2.5,
        };
        let mut l = Line::new();
        l.set_value_circle(&circ);
        assert_eq!(l.type_contour(), IType::Circle);
        let got = l.circle();
        assert!(got.center.distance(circ.center) < 1e-12);
        assert!((got.radius - 2.5).abs() < 1e-12);
        assert!(got.y_dir.distance(circ.y_dir) < 1e-12);
    }
}
