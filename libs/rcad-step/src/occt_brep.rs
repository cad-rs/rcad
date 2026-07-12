//! ASCII OCCT native BREP (`.brep`) import.
//!
//! Implements the text serialization described in *Open CASCADE Technology —
//! BRep Format* (header `DBRep_DrawableShape`, `CASCADE Topology V1/V2`,
//! geometry pools, `TShapes` section). Binary `.brep` is rejected.
//!
//! Supported: analytic surfaces (plane, cylinder, cone, sphere, torus), linear
//! extrusion, revolution, B-spline surface, rectangular trim, offset surface;
//! common 3D/2D curves; full solid topology with stored face triangulations.
//! Location indices on shape references are accepted but not yet applied (composition
//! transforms from the `Locations` section are ignored).
//! Trimmed 3D curves (record type 8) are unsupported.

use glam::{DVec2, DVec3};
use rcad_kernel::{topods, any_perpendicular, BRep};
use rcad_kernel::geom::{
    BSplineCurve2, BSplineCurve3, BSplineSurface, BezierCurve2, BezierCurve3, Circle2d, Circle3,
    ConicalSurface, Curve2d, Curve3, CylindricalSurface, Ellipse2d, Ellipse3, Hyperbola3, Line2d,
    Line3, OffsetCurve3, OffsetSurface, Parabola3, Plane, SphericalSurface, Surface3,
    ToroidalSurface, TrimmedSurface,
};
use std::borrow::Cow;
use std::path::Path;

#[derive(Debug, Clone)]
pub enum OcctBrepError {
    Io(String),
    InvalidHeader(Cow<'static, str>),
    UnexpectedEof,
    ParseFloat(String),
    ParseInt(String),
    BadToken {
        expected: Cow<'static, str>,
        got: Option<String>,
    },
    Unsupported(Cow<'static, str>),
    EmptyResult(Cow<'static, str>),
}

impl std::fmt::Display for OcctBrepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(s) => write!(f, "I/O error: {s}"),
            Self::InvalidHeader(s) => write!(f, "invalid OCCT BREP header: {s}"),
            Self::UnexpectedEof => write!(f, "unexpected end of BREP stream"),
            Self::ParseFloat(s) => write!(f, "invalid float: {s}"),
            Self::ParseInt(s) => write!(f, "invalid integer: {s}"),
            Self::BadToken { expected, got } => write!(
                f,
                "expected {expected}, got {}",
                got.as_deref().unwrap_or("<eof>")
            ),
            Self::Unsupported(s) => write!(f, "unsupported OCCT BREP feature: {s}"),
            Self::EmptyResult(s) => write!(f, "BREP import produced empty result: {s}"),
        }
    }
}

impl std::error::Error for OcctBrepError {}

fn tokenize(input: &str) -> Vec<String> {
    input.split_whitespace().map(|s| s.to_string()).collect()
}

struct Cursor<'a> {
    tok: &'a [String],
    i: usize,
}

impl<'a> Cursor<'a> {
    fn new(tok: &'a [String]) -> Self {
        Self { tok, i: 0 }
    }

    fn peek(&self) -> Option<&str> {
        self.tok.get(self.i).map(|s| s.as_str())
    }

    fn next(&mut self) -> Result<&'a str, OcctBrepError> {
        let t = self
            .tok
            .get(self.i)
            .map(|s| s.as_str())
            .ok_or(OcctBrepError::UnexpectedEof)?;
        self.i += 1;
        Ok(t)
    }

    fn expect(&mut self, lit: &'static str) -> Result<(), OcctBrepError> {
        let t = self.next()?;
        if t != lit {
            return Err(OcctBrepError::BadToken {
                expected: Cow::Borrowed(lit),
                got: Some(t.to_string()),
            });
        }
        Ok(())
    }

    fn parse_i32(&mut self) -> Result<i32, OcctBrepError> {
        let s = self.next()?.to_string();
        s.parse().map_err(|_| OcctBrepError::ParseInt(s))
    }

    fn parse_usize(&mut self) -> Result<usize, OcctBrepError> {
        let v = self.parse_i32()?;
        if v < 0 {
            return Err(OcctBrepError::ParseInt(v.to_string()));
        }
        Ok(v as usize)
    }

    fn parse_f64(&mut self) -> Result<f64, OcctBrepError> {
        let s = self.next()?.to_string();
        s.parse().map_err(|_| OcctBrepError::ParseFloat(s))
    }

    fn parse_point3(&mut self) -> Result<DVec3, OcctBrepError> {
        Ok(DVec3::new(
            self.parse_f64()?,
            self.parse_f64()?,
            self.parse_f64()?,
        ))
    }

    fn parse_dir3(&mut self) -> Result<DVec3, OcctBrepError> {
        self.parse_point3().map(|v| v.normalize_or_zero())
    }

    fn parse_point2(&mut self) -> Result<DVec2, OcctBrepError> {
        Ok(DVec2::new(self.parse_f64()?, self.parse_f64()?))
    }

    fn parse_dir2(&mut self) -> Result<DVec2, OcctBrepError> {
        self.parse_point2().map(|v| v.normalize_or_zero())
    }
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
enum LocEntry {
    Matrix([[f64; 4]; 3]),
    Compound(Vec<(i32, i32)>),
}

#[derive(Clone, Debug)]
struct Triangulation {
    nodes: Vec<DVec3>,
    triangles: Vec<[usize; 3]>,
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
enum TShapeKind {
    Vertex {
        tol: f64,
        p: DVec3,
    },
    Edge {
        tol: f64,
        curve3d: Option<(usize, f64, f64)>,
        v0: (char, i32, i32),
        v1: (char, i32, i32),
    },
    Wire {
        edges: Vec<(char, i32, i32)>,
    },
    Face {
        natural: u8,
        tol: f64,
        surface: usize,
        loc: i32,
        triangulation: Option<usize>,
    },
    Shell {
        faces: Vec<(char, i32, i32)>,
    },
    Solid {
        shells: Vec<(char, i32, i32)>,
    },
    CompSolid {
        children: Vec<(char, i32, i32)>,
    },
    Compound {
        children: Vec<(char, i32, i32)>,
    },
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct TShape {
    kind: TShapeKind,
    flags: String,
    children: Vec<(char, i32, i32)>,
}

fn parse_locations(c: &mut Cursor<'_>) -> Result<Vec<LocEntry>, OcctBrepError> {
    c.expect("Locations")?;
    let n = c.parse_usize()?;
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        match c.next()? {
            "1" => {
                let mut m = [[0f64; 4]; 3];
                for row in &mut m {
                    for cell in row.iter_mut() {
                        *cell = c.parse_f64()?;
                    }
                }
                out.push(LocEntry::Matrix(m));
            }
            "2" => {
                let mut pairs = Vec::new();
                loop {
                    let a = c.parse_i32()?;
                    if a == 0 {
                        break;
                    }
                    pairs.push((a, c.parse_i32()?));
                }
                out.push(LocEntry::Compound(pairs));
            }
            tag => {
                return Err(OcctBrepError::Unsupported(
                    format!("location record tag '{tag}'").into(),
                ));
            }
        }
    }
    Ok(out)
}

fn expand_knots(multiplicity_knots: &[(f64, i32)]) -> Vec<f64> {
    let mut knots = Vec::new();
    for &(u, q) in multiplicity_knots {
        for _ in 0..q.max(0) {
            knots.push(u);
        }
    }
    knots
}

fn parse_bspline_knot_line(c: &mut Cursor<'_>, k: usize) -> Result<Vec<(f64, i32)>, OcctBrepError> {
    let mut v = Vec::with_capacity(k);
    for _ in 0..k {
        v.push((c.parse_f64()?, c.parse_i32()?));
    }
    Ok(v)
}

fn parse_curve3(c: &mut Cursor<'_>) -> Result<Curve3, OcctBrepError> {
    let ty = c.next()?;
    match ty {
        "1" => Ok(Curve3::Line(Line3 {
            origin: c.parse_point3()?,
            direction: c.parse_dir3()?,
        })),
        "2" => {
            let center = c.parse_point3()?;
            let n = c.parse_dir3()?;
            let _dx = c.parse_dir3()?;
            let _dy = c.parse_dir3()?;
            Ok(Curve3::Circle(Circle3::new(center, n, c.parse_f64()?)))
        }
        "3" => Ok(Curve3::Ellipse(Ellipse3 {
            center: c.parse_point3()?,
            normal: c.parse_dir3()?,
            major_dir: c.parse_dir3()?,
            major_radius: c.parse_f64()?,
            minor_radius: c.parse_f64()?,
        })),
        "4" => Ok(Curve3::Parabola(Parabola3 {
            vertex: c.parse_point3()?,
            normal: c.parse_dir3()?,
            axis_dir: c.parse_dir3()?,
            focal_param: (2.0 * c.parse_f64()?).max(1e-12),
        })),
        "5" => Ok(Curve3::Hyperbola(Hyperbola3 {
            center: c.parse_point3()?,
            normal: c.parse_dir3()?,
            major_dir: c.parse_dir3()?,
            semi_major: c.parse_f64()?,
            semi_minor: c.parse_f64()?,
        })),
        "6" => {
            let rat = c.parse_i32()? != 0;
            let deg = c.parse_usize()?;
            let mut pts = Vec::with_capacity(deg + 1);
            let mut wts = Vec::with_capacity(deg + 1);
            for _ in 0..=deg {
                pts.push(c.parse_point3()?);
                wts.push(if rat { c.parse_f64()? } else { 1.0 });
            }
            Ok(Curve3::Bezier(BezierCurve3 {
                control_points: pts,
                weights: wts,
            }))
        }
        "7" => {
            let rat = c.parse_i32()? != 0;
            c.expect("0")?;
            let deg = c.parse_usize()?;
            let pole_count = c.parse_usize()?;
            let mk_count = c.parse_usize()?;
            let mut pts = Vec::with_capacity(pole_count);
            let mut wts = Vec::with_capacity(pole_count);
            for _ in 0..pole_count {
                pts.push(c.parse_point3()?);
                wts.push(if rat { c.parse_f64()? } else { 1.0 });
            }
            let mk = parse_bspline_knot_line(c, mk_count)?;
            Ok(Curve3::BSpline(BSplineCurve3 {
                degree: deg,
                knots: expand_knots(&mk),
                control_points: pts,
                weights: wts,
            }))
        }
        "8" => Err(OcctBrepError::Unsupported(
            "trimmed 3D curve (type 8) not mapped to rcad Curve3".into(),
        )),
        "9" => {
            let d = c.parse_f64()?;
            let inner = parse_curve3(c)?;
            Ok(Curve3::Offset(OffsetCurve3 {
                offset_distance: d,
                basis: Box::new(inner),
                offset_dir: DVec3::Z,
            }))
        }
        _ => Err(OcctBrepError::Unsupported(
            format!("3D curve type {ty}").into(),
        )),
    }
}

fn parse_curve2d(c: &mut Cursor<'_>) -> Result<Curve2d, OcctBrepError> {
    let ty = c.next()?;
    match ty {
        "1" => Ok(Curve2d::Line(Line2d {
            origin: c.parse_point2()?,
            direction: c.parse_dir2()?,
        })),
        "2" => {
            let center = c.parse_point2()?;
            let _dx = c.parse_dir2()?;
            let _dy = c.parse_dir2()?;
            Ok(Curve2d::Circle(Circle2d { center, x_dir: glam::DVec2::X, y_dir: glam::DVec2::Y, radius: c.parse_f64()?, }))
        }
        "3" => Ok(Curve2d::Ellipse(Ellipse2d {
            center: c.parse_point2()?,
            major_dir: c.parse_dir2()?,
            major_radius: c.parse_f64()?,
            minor_radius: c.parse_f64()?,
        })),
        "4" | "5" => Err(OcctBrepError::Unsupported(
            format!("2D curve type {ty}").into(),
        )),
        "6" => {
            let rat = c.parse_i32()? != 0;
            let deg = c.parse_usize()?;
            let mut pts = Vec::with_capacity(deg + 1);
            let mut wts = Vec::with_capacity(deg + 1);
            for _ in 0..=deg {
                pts.push(c.parse_point2()?);
                wts.push(if rat { c.parse_f64()? } else { 1.0 });
            }
            Ok(Curve2d::Bezier(BezierCurve2 {
                control_points: pts,
                weights: wts,
            }))
        }
        "7" => {
            let rat = c.parse_i32()? != 0;
            c.expect("0")?;
            let deg = c.parse_usize()?;
            let pole_count = c.parse_usize()?;
            let mk_count = c.parse_usize()?;
            let mut pts = Vec::with_capacity(pole_count);
            let mut wts = Vec::with_capacity(pole_count);
            for _ in 0..pole_count {
                pts.push(c.parse_point2()?);
                wts.push(if rat { c.parse_f64()? } else { 1.0 });
            }
            let mk = parse_bspline_knot_line(c, mk_count)?;
            Ok(Curve2d::BSpline(BSplineCurve2 {
                degree: deg,
                knots: expand_knots(&mk),
                control_points: pts,
                weights: wts,
            }))
        }
        "8" | "9" => Err(OcctBrepError::Unsupported(
            format!("2D curve type {ty}").into(),
        )),
        _ => Err(OcctBrepError::Unsupported(
            format!("2D curve type {ty}").into(),
        )),
    }
}

fn parse_surface(c: &mut Cursor<'_>) -> Result<Surface3, OcctBrepError> {
    let ty = c.next()?;
    match ty {
        "1" => {
            let p = c.parse_point3()?;
            let n = c.parse_dir3()?;
            let du = c.parse_dir3()?;
            let _dv = c.parse_dir3()?;
            Ok(Surface3::Plane(Plane::with_axes(p, n, du)))
        }
        "2" => {
            let p = c.parse_point3()?;
            let dv = c.parse_dir3()?;
            let dx = c.parse_dir3()?;
            let _dy = c.parse_dir3()?;
            Ok(Surface3::Cylinder(CylindricalSurface {
                origin: p,
                axis: dv,
                radius: c.parse_f64()?,
                ref_dir: dx,
            }))
        }
        "3" => {
            let p = c.parse_point3()?;
            let dz = c.parse_dir3()?;
            let _dx = c.parse_dir3()?;
            let _dy = c.parse_dir3()?;
            let r = c.parse_f64()?;
            let phi = c.parse_f64()?;
            Ok(Surface3::Cone(ConicalSurface {
                apex: p,
                axis: dz,
                radius: r,
                half_angle_rad: phi.abs(),
            }))
        }
        "4" => {
            let p = c.parse_point3()?;
            let dz = c.parse_dir3()?;
            let _dx = c.parse_dir3()?;
            let _dy = c.parse_dir3()?;
            Ok(Surface3::Sphere(SphericalSurface {
                center: p,
                axis: dz,
                radius: c.parse_f64()?,
                ref_dir: any_perpendicular(dz),
            }))
        }
        "5" => {
            let p = c.parse_point3()?;
            let dz = c.parse_dir3()?;
            let _dx = c.parse_dir3()?;
            let _dy = c.parse_dir3()?;
            Ok(Surface3::Torus(ToroidalSurface {
                center: p,
                axis: dz,
                major_radius: c.parse_f64()?,
                minor_radius: c.parse_f64()?,
            }))
        }
        "6" => {
            let dir = c.parse_dir3()?;
            let prof = parse_curve3(c)?;
            Ok(Surface3::LinearExtrusion(
                rcad_kernel::geom::LinearExtrusionSurface {
                    profile: Box::new(prof),
                    direction: dir.normalize_or_zero(),
                },
            ))
        }
        "7" => Ok(Surface3::Revolution(rcad_kernel::geom::RevolutionSurface {
            axis_origin: c.parse_point3()?,
            axis_dir: c.parse_dir3()?.normalize_or_zero(),
            profile: Box::new(parse_curve3(c)?),
        })),
        "8" => Err(OcctBrepError::Unsupported(
            "Bezier surface (type 8) not implemented".into(),
        )),
        "9" => {
            let ru = c.parse_i32()? != 0;
            let rv = c.parse_i32()? != 0;
            c.expect("0")?;
            c.expect("0")?;
            let mu = c.parse_usize()?;
            let mv = c.parse_usize()?;
            let nu = c.parse_usize()?;
            let nv = c.parse_usize()?;
            let ku = c.parse_usize()?;
            let kv = c.parse_usize()?;
            let mut grid = vec![vec![DVec3::ZERO; nv]; nu];
            let mut wg = vec![vec![1.0f64; nv]; nu];
            for i in 0..nu {
                for j in 0..nv {
                    grid[i][j] = c.parse_point3()?;
                    if ru || rv {
                        wg[i][j] = c.parse_f64()?;
                    }
                }
            }
            Ok(Surface3::BSpline(BSplineSurface {
                degree_u: mu,
                degree_v: mv,
                knots_u: expand_knots(&parse_bspline_knot_line(c, ku)?),
                knots_v: expand_knots(&parse_bspline_knot_line(c, kv)?),
                control_points: grid,
                weights: wg,
            }))
        }
        "10" => {
            let u0 = c.parse_f64()?;
            let u1 = c.parse_f64()?;
            let v0 = c.parse_f64()?;
            let v1 = c.parse_f64()?;
            Ok(Surface3::Trimmed(TrimmedSurface::new(
                parse_surface(c)?,
                u0,
                u1,
                v0,
                v1,
            )))
        }
        "11" => {
            let d = c.parse_f64()?;
            Ok(Surface3::Offset(OffsetSurface {
                basis: Box::new(parse_surface(c)?),
                offset_distance: d,
            }))
        }
        _ => Err(OcctBrepError::Unsupported(
            format!("surface type {ty}").into(),
        )),
    }
}

fn skip_polygon3d(c: &mut Cursor<'_>, count: usize) -> Result<(), OcctBrepError> {
    for _ in 0..count {
        let nc = c.parse_usize()?;
        let hasp = c.parse_i32()? != 0;
        let _defl = c.parse_f64()?;
        for _ in 0..nc {
            let _ = c.parse_point3()?;
        }
        if hasp {
            for _ in 0..nc {
                let _ = c.parse_f64()?;
            }
        }
    }
    Ok(())
}

fn skip_polygon_on_tris(c: &mut Cursor<'_>, count: usize) -> Result<(), OcctBrepError> {
    for _ in 0..count {
        let m = c.parse_usize()?;
        for _ in 0..m {
            let _ = c.parse_i32()?;
        }
        c.expect("p")?;
        let _d = c.parse_f64()?;
        let pflag = c.parse_i32()? != 0;
        if pflag {
            for _ in 0..m {
                let _ = c.parse_f64()?;
            }
        }
    }
    Ok(())
}

fn parse_triangulation(c: &mut Cursor<'_>) -> Result<Triangulation, OcctBrepError> {
    let node_count = c.parse_usize()?;
    let tri_count = c.parse_usize()?;
    let has_uv = c.parse_i32()? != 0;
    let _defl = c.parse_f64()?;
    let mut nodes = Vec::with_capacity(node_count);
    for _ in 0..node_count {
        nodes.push(c.parse_point3()?);
    }
    if has_uv {
        for _ in 0..node_count {
            let _ = c.parse_f64()?;
            let _ = c.parse_f64()?;
        }
    }
    let mut triangles = Vec::with_capacity(tri_count);
    for _ in 0..tri_count {
        let a = c.parse_usize()?;
        let b = c.parse_usize()?;
        let cc = c.parse_usize()?;
        if a == 0 || b == 0 || cc == 0 {
            return Err(OcctBrepError::Unsupported("bad triangulation index".into()));
        }
        triangles.push([a - 1, b - 1, cc - 1]);
    }
    Ok(Triangulation { nodes, triangles })
}

fn parse_ref_token(tok: &str) -> Result<(char, i32), OcctBrepError> {
    let mut it = tok.chars();
    let o = it
        .next()
        .ok_or_else(|| OcctBrepError::ParseInt(tok.to_string()))?;
    let rest = it.as_str();
    if o == '+' || o == '-' || o == 'i' || o == 'e' {
        let n: i32 = rest
            .parse()
            .map_err(|_| OcctBrepError::ParseInt(tok.to_string()))?;
        return Ok((o, n));
    }
    Err(OcctBrepError::BadToken {
        expected: "+/-/i/e prefixed shape ref".into(),
        got: Some(tok.to_string()),
    })
}

fn parse_subshape_refs_star(c: &mut Cursor<'_>) -> Result<Vec<(char, i32, i32)>, OcctBrepError> {
    let mut out = Vec::new();
    loop {
        let t = c.next()?;
        if t == "*" {
            break;
        }
        let (o, n) = parse_ref_token(t)?;
        let loc = c.parse_i32()?;
        // Non-zero `loc` indexes the Locations table (OCCT composition); ignored until transforms apply.
        out.push((o, n, loc));
    }
    Ok(out)
}

fn skip_edge_representations(
    c: &mut Cursor<'_>,
) -> Result<Option<(usize, f64, f64)>, OcctBrepError> {
    let mut curve3d = None;
    loop {
        let tag = c.next()?;
        if tag == "0" {
            break;
        }
        match tag {
            "1" => {
                let ci = c.parse_usize()?;
                let _loc = c.parse_i32()?;
                if _loc != 0 {
                    return Err(OcctBrepError::Unsupported(
                        "non-identity location on edge 3D curve".into(),
                    ));
                }
                let u0 = c.parse_f64()?;
                let u1 = c.parse_f64()?;
                curve3d = Some((ci, u0, u1));
            }
            "2" => {
                let _c2 = c.parse_usize()?;
                let _sf = c.parse_usize()?;
                let _lc = c.parse_i32()?;
                let _u0 = c.parse_f64()?;
                let _u1 = c.parse_f64()?;
            }
            "3" => {
                let _c2a = c.parse_usize()?;
                let _c2b = c.parse_usize()?;
                let _cont = c.next()?;
                let _sf = c.parse_usize()?;
                let _lc = c.parse_i32()?;
                let _u0 = c.parse_f64()?;
                let _u1 = c.parse_f64()?;
                let _ = c.parse_f64()?;
                let _ = c.parse_f64()?;
                let _ = c.parse_f64()?;
                let _ = c.parse_f64()?;
            }
            "4" => {
                let _cont = c.next()?;
                let _sf1 = c.parse_usize()?;
                let _lc1 = c.parse_i32()?;
                let _sf2 = c.parse_usize()?;
                let _lc2 = c.parse_i32()?;
            }
            "5" => {
                let _p = c.parse_usize()?;
                let _lc = c.parse_i32()?;
            }
            "6" => {
                let _poly_on_tri = c.parse_usize()?;
                let _tri = c.parse_usize()?;
                let _lc = c.parse_i32()?;
            }
            "7" => {
                let _p1 = c.parse_usize()?;
                let _p2 = c.parse_usize()?;
                let _tr = c.parse_usize()?;
                let _lc = c.parse_i32()?;
            }
            _ => {
                return Err(OcctBrepError::Unsupported(
                    format!("edge representation tag {tag}").into(),
                ));
            }
        }
    }
    Ok(curve3d)
}

fn is_shape_flags(s: &str) -> bool {
    s.len() == 7 && s.chars().all(|c| c == '0' || c == '1')
}

fn parse_vertex_block(c: &mut Cursor<'_>) -> Result<TShapeKind, OcctBrepError> {
    let tol = c.parse_f64()?;
    let p = c.parse_point3()?;
    loop {
        match c.peek() {
            Some("0") => {
                c.next()?;
                c.expect("0")?;
                break;
            }
            Some(_) => {
                let _u = c.parse_f64()?;
                let k = c.next()?;
                match k {
                    "1" => {
                        let _cv = c.parse_usize()?;
                    }
                    "2" => {
                        let _c2 = c.parse_usize()?;
                        let _sf = c.parse_usize()?;
                    }
                    "3" => {
                        let _v = c.parse_f64()?;
                        let _sf = c.parse_usize()?;
                    }
                    _ => {
                        return Err(OcctBrepError::Unsupported(
                            format!("vertex representation kind {k}").into(),
                        ));
                    }
                }
            }
            None => return Err(OcctBrepError::UnexpectedEof),
        }
    }
    Ok(TShapeKind::Vertex { tol, p })
}

fn parse_tshape_record(c: &mut Cursor<'_>) -> Result<TShape, OcctBrepError> {
    let kind = c.next()?;
    let tkind = match kind {
        "Ve" => parse_vertex_block(c)?,
        "Ed" => {
            let tol = c.parse_f64()?;
            let _sp = c.parse_i32()?;
            let _sr = c.parse_i32()?;
            let _deg = c.parse_i32()?;
            let curve3d = skip_edge_representations(c)?;
            let fl = c.next()?.to_string();
            if !is_shape_flags(&fl) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags after edge representations".into(),
                    got: Some(fl),
                });
            }
            let refs = parse_subshape_refs_star(c)?;
            if refs.len() != 2 {
                return Err(OcctBrepError::Unsupported(
                    "edge must reference exactly two vertices".into(),
                ));
            }
            return Ok(TShape {
                kind: TShapeKind::Edge {
                    tol,
                    curve3d,
                    v0: refs[0],
                    v1: refs[1],
                },
                flags: fl,
                children: Vec::new(),
            });
        }
        "Wi" => {
            let fl = c.next()?.to_string();
            if !is_shape_flags(&fl) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags".into(),
                    got: Some(fl),
                });
            }
            let edges = parse_subshape_refs_star(c)?;
            return Ok(TShape {
                kind: TShapeKind::Wire { edges },
                flags: fl,
                children: Vec::new(),
            });
        }
        "Fa" => {
            let natural = c.parse_i32()? as u8;
            let tol = c.parse_f64()?;
            let surface = c.parse_usize()?;
            let loc = c.parse_i32()?;
            if loc != 0 {
                return Err(OcctBrepError::Unsupported(
                    "non-zero face location not supported".into(),
                ));
            }
            let tri = if c.peek() == Some("2") {
                c.next()?;
                Some(c.parse_usize()?)
            } else {
                None
            };
            let flags = c.next()?.to_string();
            if !is_shape_flags(&flags) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags".into(),
                    got: Some(flags),
                });
            }
            let children = parse_subshape_refs_star(c)?;
            return Ok(TShape {
                kind: TShapeKind::Face {
                    natural,
                    tol,
                    surface,
                    loc,
                    triangulation: tri,
                },
                flags,
                children,
            });
        }
        "Sh" => {
            let fl = c.next()?.to_string();
            if !is_shape_flags(&fl) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags".into(),
                    got: Some(fl),
                });
            }
            let ch = parse_subshape_refs_star(c)?;
            return Ok(TShape {
                kind: TShapeKind::Shell { faces: ch },
                flags: fl,
                children: Vec::new(),
            });
        }
        "So" => {
            let fl = c.next()?.to_string();
            if !is_shape_flags(&fl) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags".into(),
                    got: Some(fl),
                });
            }
            let ch = parse_subshape_refs_star(c)?;
            return Ok(TShape {
                kind: TShapeKind::Solid { shells: ch },
                flags: fl,
                children: Vec::new(),
            });
        }
        "CS" => {
            let fl = c.next()?.to_string();
            if !is_shape_flags(&fl) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags".into(),
                    got: Some(fl),
                });
            }
            let ch = parse_subshape_refs_star(c)?;
            return Ok(TShape {
                kind: TShapeKind::CompSolid { children: ch },
                flags: fl,
                children: Vec::new(),
            });
        }
        "Co" => {
            let fl = c.next()?.to_string();
            if !is_shape_flags(&fl) {
                return Err(OcctBrepError::BadToken {
                    expected: "7-char shape flags".into(),
                    got: Some(fl),
                });
            }
            let ch = parse_subshape_refs_star(c)?;
            return Ok(TShape {
                kind: TShapeKind::Compound { children: ch },
                flags: fl,
                children: Vec::new(),
            });
        }
        _ => {
            return Err(OcctBrepError::BadToken {
                expected: "Ve|Ed|Wi|Fa|Sh|So|CS|Co".into(),
                got: Some(kind.to_string()),
            });
        }
    };
    let flags = c.next()?.to_string();
    if !is_shape_flags(&flags) {
        return Err(OcctBrepError::BadToken {
            expected: "7-char shape flags".into(),
            got: Some(flags),
        });
    }
    let children = parse_subshape_refs_star(c)?;
    Ok(TShape {
        kind: tkind,
        flags,
        children,
    })
}

fn shape_index(n_shapes: usize, back: i32) -> usize {
    n_shapes - (back.unsigned_abs() as usize)
}

fn face_normal(surfaces: &[Surface3], si: usize) -> DVec3 {
    if si == 0 || si > surfaces.len() {
        return DVec3::Z;
    }
    match &surfaces[si - 1] {
        Surface3::Plane(pl) => pl.normal.normalize_or_zero(),
        _ => DVec3::Z,
    }
}

fn triangle_normal(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
    (b - a).cross(c - a).normalize_or_zero()
}

fn build_brep(
    curves: Vec<Curve3>,
    _curve2ds: Vec<Curve2d>,
    surfaces: Vec<Surface3>,
    _tris: Vec<Triangulation>,
    shapes: Vec<TShape>,
) -> Result<BRep, OcctBrepError> {
    let n = shapes.len();
    if n == 0 {
        return Err(OcctBrepError::EmptyResult("no TShapes".into()));
    }

    let mut brep = BRep::new();

    // shape_refs[i] = ShapeRef for the TShape with parsed index i
    let mut shape_refs: Vec<Option<topods::ShapeRef>> = vec![None; n];

    // Pass 1: Create all vertex TShapes
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Vertex { p, tol } = &s.kind {
            let r = brep.add_tvertex(*p);
            brep.vertex_mut(r).tolerance = *tol;
            shape_refs[i] = Some(r);
        }
    }

    // Pass 2: Create all edge TShapes
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Edge {
            curve3d, v0, v1, tol, ..
        } = &s.kind
        {
            let i0 = shape_index(n, v0.1);
            let i1 = shape_index(n, v1.1);
            let vr0 = shape_refs[i0]
                .ok_or_else(|| OcctBrepError::EmptyResult("edge references missing vertex".into()))?;
            let vr1 = shape_refs[i1]
                .ok_or_else(|| OcctBrepError::EmptyResult("edge references missing vertex".into()))?;

            let (edge_curve, edge_range) = if let Some((ci, u0, u1)) = curve3d {
                if *ci > 0 && *ci <= curves.len() {
                    (Some(curves[*ci - 1].clone()), [*u0, *u1])
                } else {
                    (None, [0.0, 0.0])
                }
            } else {
                (None, [0.0, 0.0])
            };

            let e = brep.add_tedge(edge_curve, vr0, vr1, edge_range);
            brep.edge_mut(e).tolerance = *tol;
            shape_refs[i] = Some(e);
        }
    }

    // Pass 3: Create all wire TShapes
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Wire { edges: wrefs } = &s.kind {
            let mut edge_srs = Vec::with_capacity(wrefs.len());
            for (o, back, _loc) in wrefs {
                let si = shape_index(n, *back);
                let esr = shape_refs[si].ok_or_else(|| {
                    OcctBrepError::EmptyResult("wire references non-edge".into())
                })?;
                // Apply orientation from the reference character
                let orient = match o {
                    '+' => topods::Orientation::Forward,
                    '-' => topods::Orientation::Reversed,
                    'i' => topods::Orientation::Internal,
                    'e' => topods::Orientation::External,
                    _ => topods::Orientation::Forward,
                };
                edge_srs.push(topods::ShapeRef::synthetic_with_orientation(esr.index, orient));
            }
            let w = brep.add_twire(edge_srs);
            // Set flags from the parsed shape flags
            if s.flags.len() == 7 {
                let parsed = u16::from_str_radix(&s.flags, 2).unwrap_or(0);
                brep.wire_mut(w).flags = parsed;
            }
            shape_refs[i] = Some(w);
        }
    }

    // Pass 4: Create all face TShapes
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Face {
            natural,
            tol,
            surface,
            triangulation: _tri,
            ..
        } = &s.kind
        {
            let face_surface = if *surface > 0 && *surface <= surfaces.len() {
                Some(surfaces[*surface - 1].clone())
            } else {
                None
            };

            // Outer wire from first child
            if s.children.is_empty() {
                return Err(OcctBrepError::EmptyResult("face has no wires".into()));
            }
            let wsi = shape_index(n, s.children[0].1);
            let outer_wire_sr = shape_refs[wsi].ok_or_else(|| {
                OcctBrepError::EmptyResult("face outer is not a wire".into())
            })?;

            // Inner wires from remaining children
            let mut inner_wire_srs = Vec::new();
            for ch in s.children.iter().skip(1) {
                let wsi2 = shape_index(n, ch.1);
                let iw_sr = shape_refs[wsi2].ok_or_else(|| {
                    OcctBrepError::EmptyResult("face inner is not a wire".into())
                })?;
                // Apply orientation from the reference
                let orient = match ch.0 {
                    '+' => topods::Orientation::Forward,
                    '-' => topods::Orientation::Reversed,
                    'i' => topods::Orientation::Internal,
                    'e' => topods::Orientation::External,
                    _ => topods::Orientation::Forward,
                };
                inner_wire_srs.push(topods::ShapeRef::synthetic_with_orientation(iw_sr.index, orient));
            }

            let natural_flag = *natural != 0;
            let f = brep.add_tface(
                face_surface,
                outer_wire_sr,
                inner_wire_srs,
                None,   // sample_point
                None,   // uv_domain
                vec![], // internal_vertices
                natural_flag,
            );
            brep.face_mut(f).tolerance = *tol;
            shape_refs[i] = Some(f);
        }
    }

    // Pass 5: Create all shell TShapes
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Shell { faces: frefs } = &s.kind {
            let mut face_srs = Vec::with_capacity(frefs.len());
            for (_o, back, _loc) in frefs {
                let fi = shape_index(n, *back);
                let fsr = shape_refs[fi].ok_or_else(|| {
                    OcctBrepError::EmptyResult("shell references non-face".into())
                })?;
                face_srs.push(fsr);
            }
            let sh = brep.add_tshell(face_srs);
            shape_refs[i] = Some(sh);
        }
    }

    // Pass 6: Create all solid TShapes
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Solid { shells: srefs } = &s.kind {
            let mut shell_srs = Vec::with_capacity(srefs.len());
            for (_o, back, _loc) in srefs {
                let shi = shape_index(n, *back);
                let ssr = shape_refs[shi].ok_or_else(|| {
                    OcctBrepError::EmptyResult("solid references non-shell".into())
                })?;
                shell_srs.push(ssr);
            }
            brep.add_tsolid(shell_srs);
            // No need to store shape_refs for solids — we don't need to reference them further
        }
    }

    // Handle CompSolid and Compound shapes (referenced by other top-level shapes)
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::CompSolid { children: refs } = &s.kind {
            let mut child_srs = Vec::with_capacity(refs.len());
            for (_o, back, _loc) in refs {
                let ci = shape_index(n, *back);
                if let Some(csr) = shape_refs[ci] {
                    child_srs.push(csr);
                }
            }
            brep.add_tcompsolid(child_srs);
        }
    }
    for (i, s) in shapes.iter().enumerate() {
        if let TShapeKind::Compound { children: refs } = &s.kind {
            let mut child_srs = Vec::with_capacity(refs.len());
            for (_o, back, _loc) in refs {
                let ci = shape_index(n, *back);
                if let Some(csr) = shape_refs[ci] {
                    child_srs.push(csr);
                }
            }
            brep.add_tcompound(child_srs);
        }
    }

    if brep.solid_count() == 0 {
        return Err(OcctBrepError::EmptyResult("no solid found in BREP".into()));
    }

    Ok(brep)
}

/// Reader for OCCT ASCII `.brep` files.
pub struct OcctBrepReader;

impl OcctBrepReader {
    /// Read an ASCII OCCT `.brep` file from disk.
    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, OcctBrepError> {
        let s =
            std::fs::read_to_string(path.as_ref()).map_err(|e| OcctBrepError::Io(e.to_string()))?;
        Self::parse_string(&s)
    }

    /// Parse ASCII OCCT BREP from a string (UTF-8).
    pub fn parse_string(content: &str) -> Result<BRep, OcctBrepError> {
        if content.as_bytes().contains(&0) {
            return Err(OcctBrepError::InvalidHeader(
                "binary BREP is not supported".into(),
            ));
        }
        let t = tokenize(content);
        let mut c = Cursor::new(&t);
        if c.next()? != "DBRep_DrawableShape" {
            return Err(OcctBrepError::InvalidHeader(
                "expected DBRep_DrawableShape".into(),
            ));
        }
        // Version is one logical line but tokenizes to multiple words (e.g. "CASCADE Topology V1, (c) Matra-Datavision").
        let mut ver = String::new();
        loop {
            let t = c.peek().ok_or(OcctBrepError::InvalidHeader(
                "expected CASCADE Topology line before Locations".into(),
            ))?;
            if t == "Locations" {
                break;
            }
            if !ver.is_empty() {
                ver.push(' ');
            }
            ver.push_str(c.next()?);
        }
        if !ver.contains("CASCADE") || !ver.contains("Topology") {
            return Err(OcctBrepError::InvalidHeader(
                "expected CASCADE Topology V1/V2 line".into(),
            ));
        }
        let _locs = parse_locations(&mut c)?;
        c.expect("Curve2ds")?;
        let n2 = c.parse_usize()?;
        let mut curve2ds = Vec::with_capacity(n2);
        for _ in 0..n2 {
            curve2ds.push(parse_curve2d(&mut c)?);
        }
        c.expect("Curves")?;
        let nc = c.parse_usize()?;
        let mut curves = Vec::with_capacity(nc);
        for _ in 0..nc {
            curves.push(parse_curve3(&mut c)?);
        }
        c.expect("Polygon3D")?;
        let np3 = c.parse_usize()?;
        skip_polygon3d(&mut c, np3)?;
        c.expect("PolygonOnTriangulations")?;
        let npt = c.parse_usize()?;
        skip_polygon_on_tris(&mut c, npt)?;
        c.expect("Surfaces")?;
        let ns = c.parse_usize()?;
        let mut surfaces = Vec::with_capacity(ns);
        for _ in 0..ns {
            surfaces.push(parse_surface(&mut c)?);
        }
        c.expect("Triangulations")?;
        let nt = c.parse_usize()?;
        let mut tris = Vec::with_capacity(nt);
        for _ in 0..nt {
            tris.push(parse_triangulation(&mut c)?);
        }
        c.expect("TShapes")?;
        let nsh = c.parse_usize()?;
        let mut shapes = Vec::with_capacity(nsh);
        for _ in 0..nsh {
            shapes.push(parse_tshape_record(&mut c)?);
        }
        let _ = c.next()?;
        let _ = c.next()?;
        build_brep(curves, curve2ds, surfaces, tris, shapes)
    }
}

#[cfg(test)]
mod tests {
    use super::OcctBrepReader;

    #[test]
    fn rejects_bad_header() {
        assert!(OcctBrepReader::parse_string("not a brep").is_err());
    }

    #[test]
    fn rejects_binary_marker() {
        let mut s = String::from("DBRep_DrawableShape\n");
        s.push('\0');
        assert!(OcctBrepReader::parse_string(&s).is_err());
    }
}
