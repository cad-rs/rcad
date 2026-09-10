//! OCCT BRepOffsetAPI_ThruSections (TKOffset) — loft between section wires.
//!
//! OCCT BRepOffsetAPI_ThruSections builds a shell/solid interpolating the
//! given section wires (BRepOffsetAPI_ThruSections.cxx, via BRepFill).  This
//! port implements the RULED loft used when exactly two sections are given or
//! myIsRuled is true (BRepOffsetAPI_ThruSections::Build L514-523 ->
//! CreateRuled L542-702 -> BRepFill_Generator::Perform L613-1186):
//!
//!   - the working sections are paired edge-by-edge (BRepFill_CompatibleWires;
//!     this port requires equal edge counts and pairs by index — the ported
//!     GTests construct compatible wires),
//!   - each edge pair generates one ruled face: the V-degree-1 BSpline
//!     surface whose two V-poles are the profiled section curves
//!     (GeomFill_Generator::Perform L32-81; GeomFill_Profiler::Perform
//!     L166-246 unifies degree, parameter range and knot vectors),
//!   - for isSolid the extremities are closed with planar faces fitted from
//!     the end wires (MakeSolid L191-254 + PerformPlan L119-157) and the
//!     solid is oriented via BRepClass3d_SolidClassifier::PerformInfinitePoint.
//!
//! The smoothed loft (CreateSmoothed L706-1116, GeomFill_AppSurf) is not
//! ported.
//!
//! OCCT file refs: BRepOffsetAPI_ThruSections.cxx, BRepFill_Generator.cxx,
//! BRepFill_CompatibleWires.cxx, GeomFill_Generator.cxx,
//! GeomFill_Profiler.cxx.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{
    BSplineCurve3, BSplineSurface, Circle3, Curve2d, Curve3, Line2d, Line3, Surface3,
};
use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topo::topods::{self, BRep, Orientation, TShape};
use std::f64::consts::PI;

use crate::topalgo::brep_class3d::solid_classifier::SolidClassifier;

/// OCCT BRepOffsetAPI_ThruSections — section-wire loft builder.
#[derive(Debug, Default)]
pub struct ThruSections {
    /// OCCT isSolid: whether to build a solid (vs an open shell).
    pub is_solid: bool,
    /// OCCT ruled: ruled (vs smooth) interpolation.
    pub ruled: bool,
    /// OCCT pres3d: 3D tolerance for the loft.
    pub pres3d: f64,
    wires: Vec<Shape>,
    done: bool,
    shape: Option<Shape>,
}

impl ThruSections {
    /// OCCT BRepOffsetAPI_ThruSections(isSolid, ruled, pres3d).
    pub fn new(is_solid: bool, ruled: bool, pres3d: f64) -> Self {
        ThruSections {
            is_solid,
            ruled,
            pres3d,
            wires: Vec::new(),
            done: false,
            shape: None,
        }
    }

    /// OCCT BRepOffsetAPI_ThruSections::AddWire(wire).
    pub fn add_wire(&mut self, wire: Shape) {
        self.wires.push(wire);
    }

    /// OCCT BRepOffsetAPI_ThruSections::Build (BRepOffsetAPI_ThruSections.cxx
    /// L339-538).  The section compatibility pass (BRepFill_CompatibleWires,
    /// L380-509) is covered for the ported cases by the index pairing in
    /// create_ruled; with two sections or myIsRuled the ruled shell is built
    /// (L514-518 CreateRuled), otherwise the smoothed path is not ported.
    pub fn build(&mut self, brep: &mut BRep) {
        self.done = false;
        self.shape = None;
        if self.wires.len() < 2 {
            return;
        }
        // BRepOffsetAPI_ThruSections::Build L514-523: only the ruled path.
        if self.wires.len() != 2 && !self.ruled {
            return;
        }
        let shell = match create_ruled(brep, &self.wires) {
            Some(s) => s,
            None => return,
        };
        let shape = if self.is_solid {
            // MakeSolid (L191-254): close the extremities with planar faces
            // and orient the solid.
            match make_solid(brep, &shell, &self.wires) {
                Some(s) => s,
                None => return,
            }
        } else {
            shell
        };
        self.shape = Some(shape);
        self.done = true;
    }

    /// OCCT BRepOffsetAPI_ThruSections::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BRepOffsetAPI_ThruSections::Shape().
    pub fn shape(&self) -> Option<Shape> {
        self.shape.clone()
    }
}

// =============================================================================
// BRepFill_Generator (BRepFill_Generator.cxx L613-1186) — ruled faces
// =============================================================================

/// The edges of a wire in wire order (BRepTools_WireExplorer).
fn wire_edges(brep: &BRep, wire: &Shape) -> Vec<Shape> {
    match wire.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// True when the wire is closed (the wire flag, or first == last vertex).
fn wire_is_closed(brep: &BRep, wire: &Shape) -> bool {
    let edges = wire_edges(brep, wire);
    if edges.is_empty() {
        return false;
    }
    if let TShape::Wire(wd) = wire.data.as_ref() {
        if wd.flags & topods::tshape_flags::CLOSED != 0 {
            return true;
        }
    }
    // OCCT BRepFill_CompatibleWires L942-955: a wire is also closed when its
    // first and last vertices are the same.
    let (_, l) = edge_vertices(brep, &edges[edges.len() - 1]);
    let (f, _) = edge_vertices(brep, &edges[0]);
    l.is_same(&f)
}

/// (traversal start, traversal end) vertices of an edge in its wire — the
/// orientation of the edge wrapper decides (BRepFill_Generator.cxx L736-751).
fn edge_vertices(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    let ed = brep.edge(e.clone());
    if e.orientation == Orientation::Reversed {
        (ed.last.clone(), ed.first.clone())
    } else {
        (ed.first.clone(), ed.last.clone())
    }
}

/// The first and last vertices of a wire in its traversal direction
/// (TopExp::Vertices(wire), BRepFill_CompatibleWires L2498-2500).
fn wire_endpoints(brep: &BRep, edges: &[Shape]) -> (DVec3, DVec3) {
    let (s, _) = edge_vertices(brep, &edges[0]);
    let (_, e) = edge_vertices(brep, &edges[edges.len() - 1]);
    (brep.vertex(s.clone()).point, brep.vertex(e.clone()).point)
}

/// True when the wire's points are collinear (PlaneOfWire fails — the
/// straight-segment branch of SearchOrigin, L2531-2540).
fn wire_is_straight(brep: &BRep, edges: &[Shape]) -> bool {
    let mut pts: Vec<DVec3> = Vec::new();
    for e in edges {
        let ed = brep.edge(e.clone());
        let push = |p: DVec3, pts: &mut Vec<DVec3>| {
            if !pts.iter().any(|q| (q - p).length_squared() < 1e-18) {
                pts.push(p);
            }
        };
        push(brep.vertex(ed.first.clone()).point, &mut pts);
        push(brep.vertex(ed.last.clone()).point, &mut pts);
    }
    if pts.len() < 3 {
        return true;
    }
    let d0 = pts[1] - pts[0];
    let len0 = d0.length();
    if len0 < 1e-30 {
        return true;
    }
    let dir = d0 / len0;
    pts.iter().skip(2).all(|p| {
        let d = (*p - pts[0]).cross(dir).length();
        d < 1e-6 * (p - pts[0]).length().max(1.0)
    })
}

/// The 3D curve of an edge in world coordinates, reversed together with its
/// range for a REVERSED edge (BRepFill_Generator.cxx L787-849).
fn edge_curve_oriented(brep: &BRep, e: &Shape) -> Option<(Curve3, [f64; 2])> {
    let ed = brep.edge(e.clone());
    let c = ed.curve.clone()?;
    let mut r = ed.range;
    if e.orientation == Orientation::Reversed {
        r = [-r[1], -r[0]];
        Some((curve_reversed(&c), r))
    } else {
        Some((c, r))
    }
}

/// OCCT Geom_Curve::Reverse — the reversed parameterization.
fn curve_reversed(c: &Curve3) -> Curve3 {
    match c {
        Curve3::Line(l) => Curve3::Line(Line3::new(l.origin + l.direction, -l.direction)),
        Curve3::Circle(cir) => {
            let mut out = *cir;
            out.y_dir = -out.y_dir;
            Curve3::Circle(out)
        }
        Curve3::BSpline(b) => {
            let mut b = b.clone();
            b.control_points.reverse();
            b.weights.reverse();
            // Mirror the interior knots about the range center and reverse
            // their order: evaluating the reversed curve at -t over the
            // mirrored range reproduces the original geometry.
            let (k0, k1) = (b.knots[b.degree], b.knots[b.knots.len() - b.degree - 1]);
            b.knots = b.knots.iter().rev().map(|k| (k1 + k0) - k).collect();
            Curve3::BSpline(b)
        }
        _ => c.clone(),
    }
}

/// BRepFill_Generator::Perform for the two-section case: one ruled face per
/// corresponding edge pair of the two working wires.
fn create_ruled(brep: &mut BRep, wires: &[Shape]) -> Option<Shape> {
    let w1 = wires.first()?.clone();
    let w2 = wires.last()?.clone();
    let mut e1s = wire_edges(brep, &w1);
    let mut e2s = wire_edges(brep, &w2);
    if e1s.is_empty() || e1s.len() != e2s.len() {
        return None;
    }
    // BRepFill_CompatibleWires::SearchOrigin (L2459-2628): open sections are
    // reorganized so consecutive wires run in the same direction — the second
    // wire is reversed when the angle between the section direction vectors
    // is >= PI/2 (or, for straight segments, when the crossed distance sum is
    // smaller).  Closed wires are handled by ComputeOrigin (index pairing is
    // sufficient for the ported cases).
    let open1 = !wire_is_closed(brep, &w1);
    let open2 = !wire_is_closed(brep, &w2);
    if open1 && open2 {
        let (p1o, p2o) = wire_endpoints(brep, &e1s);
        let (p1, p2) = wire_endpoints(brep, &e2s);
        let straight = wire_is_straight(brep, &e1s) || wire_is_straight(brep, &e2s);
        let parcours = if straight {
            let dist1 = (p1o - p1).length() + (p2o - p2).length();
            let dist2 = (p1o - p2).length() + (p2o - p1).length();
            dist2 >= dist1
        } else {
            let v0 = p2o - p1o;
            let v = p2 - p1;
            v0.angle_between(v) < PI / 2.0
        };
        if !parcours {
            e2s = e2s
                .iter()
                .rev()
                .map(|e| Shape {
                    orientation: if e.orientation == Orientation::Reversed {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    },
                    ..e.clone()
                })
                .collect();
        }
    }
    let n = e1s.len();
    let closed = wire_is_closed(brep, &w1) && wire_is_closed(brep, &w2);

    // BRepFill_Generator.cxx L673-761: vertex roles of every edge pair.
    let mut v1f = Vec::with_capacity(n);
    let mut v1l = Vec::with_capacity(n);
    let mut v2f = Vec::with_capacity(n);
    let mut v2l = Vec::with_capacity(n);
    for (e1, e2) in e1s.iter().zip(e2s.iter()) {
        let (a, b) = edge_vertices(brep, e1);
        v1f.push(a);
        v1l.push(b);
        let (a, b) = edge_vertices(brep, e2);
        v2f.push(a);
        v2l.push(b);
    }

    // The longitudinal edges (L873-971): the left edge of face i connects
    // V1f_i -> V2f_i (the surface iso U=f1).  For closed wires the right
    // edge of face i is the left edge of face (i+1)%n (Map-based sharing,
    // L884-971); for open wires each face owns both edges.
    let mut left: Vec<Shape> = Vec::with_capacity(n);
    let mut right: Vec<Shape> = Vec::with_capacity(n);
    let rev = |sr: Shape| Shape {
        orientation: Orientation::Reversed,
        ..sr
    };
    for i in 0..n {
        let p1 = brep.vertex(v1f[i].clone()).point;
        let p2 = brep.vertex(v2f[i].clone()).point;
        let dir = p2 - p1;
        let len = dir.length();
        if len < 1e-30 {
            return None;
        }
        left.push(brep.add_tedge(
            Some(Curve3::Line(Line3::new(p1, dir / len))),
            v1f[i].clone(),
            rev(v2f[i].clone()),
            [0.0, len],
        ));
    }
    if closed {
        for i in 0..n {
            right.push(left[(i + 1) % n].clone());
        }
    } else {
        for i in 0..n {
            let p1 = brep.vertex(v1l[i].clone()).point;
            let p2 = brep.vertex(v2l[i].clone()).point;
            let dir = p2 - p1;
            let len = dir.length();
            if len < 1e-30 {
                return None;
            }
            right.push(brep.add_tedge(
                Some(Curve3::Line(Line3::new(p1, dir / len))),
                v1l[i].clone(),
                rev(v2l[i].clone()),
                [0.0, len],
            ));
        }
    }

    // Faces (L851-1104): one per edge pair, boundary
    // [Edge1, right, Edge2.Reversed(), left.Reversed] (L1092-1101).
    let mut faces: Vec<Shape> = Vec::with_capacity(n);
    for i in 0..n {
        let e1 = &e1s[i];
        let e2 = &e2s[i];
        let (c1, r1) = edge_curve_oriented(brep, e1)?;
        let (c2, r2) = edge_curve_oriented(brep, e2)?;
        // GeomFill_Generator (L851-857): the ruled surface between the two
        // section curves.
        let surface = ruled_surface(&c1, r1, &c2, r2)?;
        let wire = brep.add_twire(vec![
            e1.clone(),
            right[i].clone(),
            rev(e2.clone()),
            rev(left[i].clone()),
        ]);
        let face = brep.add_tface(
            Some(surface),
            wire,
            vec![],
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            vec![],
            false,
        );
        faces.push(face);
    }

    // Pcurves (L1005-1101): straight lines in the face UV domain [0,1]^2.
    // edge_mut_inplace (not edge_mut) so the Arc identity survives: the
    // section and longitudinal edges are already referenced by the wires,
    // and a clone-on-write would leave the wire edges without their pcurves
    // (BRep_Builder::UpdateEdge semantics, see topods.rs).
    for i in 0..n {
        let face = &faces[i];
        let fkey = (face.ptr_id(), face.location);
        // Edge1: (t, f2=0) — the U axis at V = 0.
        brep.edge_mut_inplace(e1s[i].clone())
            .pcurves
            .insert(fkey, (Curve2d::Line(Line2d::new(DVec3::ZERO.truncate(), DVec3::X.truncate())), 0.0, 1.0));
        // Edge2: (t, l2=1).
        brep.edge_mut_inplace(e2s[i].clone())
            .pcurves
            .insert(fkey, (Curve2d::Line(Line2d::new(DVec3::new(0.0, 1.0, 0.0).truncate(), DVec3::X.truncate())), 0.0, 1.0));
        // left edge: (f1=0, t).
        brep.edge_mut_inplace(left[i].clone())
            .pcurves
            .insert(fkey, (Curve2d::Line(Line2d::new(DVec3::ZERO.truncate(), DVec3::Y.truncate())), 0.0, 1.0));
        // right edge: (l1=1, t).
        brep.edge_mut_inplace(right[i].clone())
            .pcurves
            .insert(fkey, (Curve2d::Line(Line2d::new(DVec3::new(1.0, 0.0, 0.0).truncate(), DVec3::Y.truncate())), 0.0, 1.0));
    }

    Some(brep.add_tshell(faces))
}

// =============================================================================
// GeomFill_Generator (GeomFill_Generator.cxx L32-81) + GeomFill_Profiler
// (GeomFill_Profiler.cxx L166-246) — the ruled surface between two curves
// =============================================================================

/// OCCT GeomFill_Profiler::AddCurve (L125-162) + GeomConvert: convert an edge
/// 3D curve to its (rational) BSpline form over [0,1].  Lines become
/// degree-1, circular arcs exact degree-2 rational, B-splines pass through.
fn curve_to_bspline(c: &Curve3, range: [f64; 2]) -> Option<BSplineCurve3> {
    match c {
        Curve3::Line(l) => {
            let p0 = l.origin;
            let p1 = l.origin + l.direction * (range[1] - range[0]);
            Some(BSplineCurve3 {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![p0, p1],
                weights: vec![1.0, 1.0],
                is_periodic: false,
            })
        }
        Curve3::Circle(cir) => circle_arc_bspline(cir, range),
        Curve3::BSpline(b) => {
            let mut b = b.clone();
            // Reparametrize to [0,1] (the caller unifies both curves).
            reparametrize(&mut b, 0.0, 1.0);
            Some(b)
        }
        Curve3::Bezier(bz) => {
            let d = bz.control_points.len().checked_sub(1)?;
            let mut knots = Vec::with_capacity(2 * (d + 1));
            for _ in 0..=d {
                knots.push(0.0);
            }
            for _ in 0..=d {
                knots.push(1.0);
            }
            Some(BSplineCurve3 {
                degree: d,
                knots,
                control_points: bz.control_points.clone(),
                weights: bz.weights.clone(),
                is_periodic: false,
            })
        }
        _ => None,
    }
}

/// Exact rational degree-2 BSpline of a circular arc (segments of at most
/// PI/2, weights [1, cos(seg/2), 1], interior knot multiplicity 2).
fn circle_arc_bspline(cir: &Circle3, range: [f64; 2]) -> Option<BSplineCurve3> {
    let [t1, t2] = range;
    let delta = t2 - t1;
    if delta <= 0.0 || delta > 2.0 * PI + 1e-9 {
        return None;
    }
    let nseg = (delta / (PI / 2.0)).ceil().max(1.0) as usize;
    let seg = delta / nseg as f64;
    let (o, r, x, y) = (cir.center, cir.radius, cir.x_dir, cir.y_dir);
    let pt = |a: f64| o + r * (a.cos() * x + a.sin() * y);
    let mut poles: Vec<DVec3> = Vec::with_capacity(2 * nseg + 1);
    let mut weights: Vec<f64> = Vec::with_capacity(2 * nseg + 1);
    let mut knots: Vec<f64> = Vec::with_capacity(2 * nseg + 4);
    for s in 0..nseg {
        let a = t1 + s as f64 * seg;
        let m = a + seg / 2.0;
        let w = (seg / 2.0).cos();
        if s == 0 {
            poles.push(pt(a));
            weights.push(1.0);
        }
        poles.push(o + r * (m.cos() * x + m.sin() * y) / w);
        weights.push(w);
        poles.push(pt(a + seg));
        weights.push(1.0);
    }
    knots.push(0.0);
    knots.push(0.0);
    knots.push(0.0);
    for s in 1..nseg {
        let u = s as f64 / nseg as f64;
        knots.push(u);
        knots.push(u);
    }
    knots.push(1.0);
    knots.push(1.0);
    knots.push(1.0);
    Some(BSplineCurve3 {
        degree: 2,
        knots,
        control_points: poles,
        weights,
        is_periodic: false,
    })
}

/// OCCT BSplCLib::Reparametrize — affine rescale of the knot vector.
fn reparametrize(c: &mut BSplineCurve3, t0: f64, t1: f64) {
    let d = c.degree;
    if c.knots.len() <= 2 * d {
        return;
    }
    let k0 = c.knots[d];
    let k1 = c.knots[c.knots.len() - d - 1];
    let span = k1 - k0;
    if span.abs() < 1e-300 {
        return;
    }
    for k in c.knots.iter_mut() {
        *k = t0 + (*k - k0) * (t1 - t0) / span;
    }
}

/// GeomFill_Profiler::UnifyByInsertingAllKnots (L29-51): insert the interior
/// knots of each curve into the other so both share one knot vector.
fn unify_knots(a: &mut BSplineCurve3, b: &mut BSplineCurve3) {
    let ia = interior_knots(a);
    let ib = interior_knots(b);
    for (u, m) in &ib {
        insert_knots_h(a, *u, *m);
    }
    for (u, m) in &ia {
        insert_knots_h(b, *u, *m);
    }
}

/// Interior (knot value, multiplicity) pairs of a clamped BSpline.
fn interior_knots(c: &BSplineCurve3) -> Vec<(f64, usize)> {
    let mut out: Vec<(f64, usize)> = Vec::new();
    for (i, &k) in c.knots.iter().enumerate() {
        if i < c.degree + 1 || i >= c.knots.len() - c.degree - 1 {
            continue;
        }
        match out.last_mut() {
            Some((last, m)) if (*last - k).abs() < 1e-15 => *m += 1,
            _ => out.push((k, 1)),
        }
    }
    out
}

/// Boehm knot insertion on homogeneous coordinates (P&T A5.1) — raises the
/// multiplicity of the interior knot `u` to `m` (max degree).
fn insert_knots_h(c: &mut BSplineCurve3, u: f64, m: usize) {
    let d = c.degree;
    if m == 0 {
        return;
    }
    let cur = c
        .knots
        .iter()
        .filter(|&&k| (k - u).abs() < 1e-12)
        .count();
    if cur >= m || cur >= d {
        return;
    }
    for _ in cur..m {
        if !insert_knot_one(c, u) {
            break;
        }
    }
}

/// Insert a single copy of the interior knot `u` (Boehm, P&T A5.1); returns
/// false when the insertion is not interior (clamped boundary).
fn insert_knot_one(c: &mut BSplineCurve3, u: f64) -> bool {
    let d = c.degree;
    let n = c.control_points.len();
    let knots = &c.knots; // len = n + d + 1
    // 0-based span index: the LAST knot <= u (u interior, so 0 < k < n+d+1).
    let mut k = 0usize;
    for (i, &tk) in knots.iter().enumerate() {
        if tk <= u + 1e-12 {
            k = i;
        } else {
            break;
        }
    }
    // k must satisfy d <= k <= n (1-based span k+1 in [d+1, n+1]).
    if k < d {
        k = d;
    }
    if k > n {
        k = n;
    }
    if u <= knots[d] + 1e-12 || u >= knots[n + d] - 1e-12 {
        return false;
    }
    // Homogeneous control points.
    let h: Vec<(DVec3, f64)> = (0..n)
        .map(|i| (c.control_points[i] * c.weights[i], c.weights[i]))
        .collect();
    let mut new_h: Vec<(DVec3, f64)> = Vec::with_capacity(n + 1);
    // 1-based affected pole range [k-d+2, k+1] (P&T A5.1 with the 1-based
    // span index k+1).
    let lo = k + 2 - d;
    let hi = k + 1;
    for i1 in 1..=(n + 1) {
        // 1-based new pole index.
        let new_pole = if i1 < lo {
            h[i1 - 1].clone()
        } else if i1 > hi {
            h[i1 - 2].clone()
        } else {
            // Q_i1 = alpha P_i1 + (1-alpha) P_{i1-1},
            // alpha = (u - U[i1]) / (U[i1+d] - U[i1]) — 1-based U.
            let ui = knots[i1 - 1];
            let uid = knots[i1 - 1 + d];
            let denom = uid - ui;
            let alpha = if denom.abs() < 1e-300 { 1.0 } else { (u - ui) / denom };
            let (pi, wi) = &h[i1 - 1];
            let (pim, wim) = &h[i1 - 2];
            (*pi * alpha + *pim * (1.0 - alpha), *wi * alpha + *wim * (1.0 - alpha))
        };
        new_h.push(new_pole);
    }
    let mut new_knots: Vec<f64> = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(u);
    new_knots.extend_from_slice(&knots[k + 1..]);
    c.knots = new_knots;
    c.control_points = new_h.iter().map(|(p, w)| *p / *w).collect();
    c.weights = new_h.iter().map(|(_, w)| *w).collect();
    true
}

/// GeomFill_Generator::Perform (L32-81): the V-degree-1 BSpline surface whose
/// V-poles are the two unified section curves — the ruled surface
/// S(u,v) = (1-v) C1(u) + v C2(u).
fn ruled_surface(c1: &Curve3, r1: [f64; 2], c2: &Curve3, r2: [f64; 2]) -> Option<Surface3> {
    let mut b1 = curve_to_bspline(c1, r1)?;
    let mut b2 = curve_to_bspline(c2, r2)?;
    // GeomFill_Profiler::Perform L192-215: common degree, then reparameterize
    // every curve to the same range.
    let degree = b1.degree.max(b2.degree);
    if b1.degree < degree || b2.degree < degree {
        // increase_degree preserves poles only for non-rational curves; the
        // ported sections never mix rational arcs with a higher degree.
        let rational = |b: &BSplineCurve3| b.weights.iter().any(|w| (*w - 1.0).abs() > 1e-12);
        if rational(&b1) || rational(&b2) {
            return None;
        }
        if b1.degree < degree {
            b1.increase_degree(degree);
        }
        if b2.degree < degree {
            b2.increase_degree(degree);
        }
    }
    reparametrize(&mut b1, 0.0, 1.0);
    reparametrize(&mut b2, 0.0, 1.0);
    // UnifyByInsertingAllKnots (L223-243).
    unify_knots(&mut b1, &mut b2);
    let nu = b1.control_points.len();
    if nu != b2.control_points.len() || nu == 0 {
        return None;
    }
    let control_points: Vec<Vec<DVec3>> = (0..nu)
        .map(|i| vec![b1.control_points[i], b2.control_points[i]])
        .collect();
    let weights: Vec<Vec<f64>> = (0..nu)
        .map(|i| vec![b1.weights[i], b2.weights[i]])
        .collect();
    Some(Surface3::BSpline(BSplineSurface {
        degree_u: degree,
        degree_v: 1,
        knots_u: b1.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        weights,
    }))
}

// =============================================================================
// MakeSolid (BRepOffsetAPI_ThruSections.cxx L191-254) + PerformPlan (L119-157)
// =============================================================================

/// BRepOffsetAPI_ThruSections::MakeSolid: close the shell extremities with
/// planar faces from the end wires, then orient the solid via the infinite
/// point classification (L237-253).
fn make_solid(brep: &mut BRep, shell: &Shape, wires: &[Shape]) -> Option<Shape> {
    let w1 = wires.first()?;
    let w2 = wires.last()?;
    let mut faces = match shell.data.as_ref() {
        TShape::Shell(sd) => sd.faces.clone(),
        _ => return None,
    };
    // PerformPlan (L206-234): the end caps, oriented opposite to the band.
    if let Some(cap) = perform_plan(brep, w1) {
        if is_same_oriented(brep, &cap, &faces) {
            faces.push(Shape {
                orientation: Orientation::Reversed,
                ..cap
            });
        } else {
            faces.push(cap);
        }
    }
    if let Some(cap) = perform_plan(brep, w2) {
        if is_same_oriented(brep, &cap, &faces) {
            faces.push(Shape {
                orientation: Orientation::Reversed,
                ..cap
            });
        } else {
            faces.push(cap);
        }
    }
    let shell = brep.add_tshell(faces);
    let mut solid = brep.add_tsolid(vec![shell.clone()]);
    // L242-250: BRepClass3d_SolidClassifier PerformInfinitePoint — reverse
    // the shell when the infinite point is inside.
    let mut clas3d = SolidClassifier::from_shape(&solid);
    clas3d.perform_infinite_point(rcad_kernel::core::precision::CONFUSION);
    if clas3d.state() == 3 {
        // IN
        solid = brep.add_tsolid(vec![Shape {
            orientation: Orientation::Reversed,
            ..shell
        }]);
    }
    Some(solid)
}

/// BRepOffsetAPI_ThruSections::IsSameOriented (L164-187): the cap must
/// traverse its shared edge opposite to the adjacent shell face.  Returns
/// true when the cap's first edge runs OPPOSITE to the shell face's edge —
/// the cap is then kept; otherwise it is reversed by the caller.
fn is_same_oriented(brep: &BRep, cap: &Shape, shell_faces: &[Shape]) -> bool {
    let Some(first) = first_wire_edge_of_face(brep, cap) else {
        return true;
    };
    for f in shell_faces {
        for e in face_wire_edges(brep, f) {
            if e.is_same(&first) {
                // OCCT L185-186: Or1 (cap) != Or2 (shell) — the two faces
                // traverse the edge in opposite directions.
                return first.orientation != e.orientation;
            }
        }
    }
    true
}

fn first_wire_edge_of_face(brep: &BRep, face: &Shape) -> Option<Shape> {
    match face.data.as_ref() {
        TShape::Face(fd) => {
            let wd = brep.wire(fd.outer_wire.clone());
            let e = wd.edges.first()?.clone();
            Some(Shape {
                orientation: face.orientation.compose(e.orientation),
                ..e
            })
        }
        _ => None,
    }
}

fn face_wire_edges(brep: &BRep, face: &Shape) -> Vec<Shape> {
    match face.data.as_ref() {
        TShape::Face(fd) => brep.wire(fd.outer_wire.clone()).edges.clone(),
        _ => Vec::new(),
    }
}

/// BRepOffsetAPI_ThruSections::PerformPlan (L119-157): a planar face on the
/// wire (BRepBuilderAPI_FindPlane + BRepBuilderAPI_MakeFace).  Falls back to
/// rcad's make_face_from_wire_brep for polygon wires; single-edge closed
/// profiles (B-spline sections) fit the plane from curve samples.
fn perform_plan(brep: &mut BRep, wire: &Shape) -> Option<Shape> {
    use rcad_kernel::geom::CurveEval;
    let edges = wire_edges(brep, wire);
    if edges.is_empty() {
        return None;
    }
    // Collect distinct points: edge vertices plus curve samples (the closed
    // B-spline profile has a single edge whose endpoints coincide).
    let mut pts: Vec<DVec3> = Vec::new();
    let mut push = |p: DVec3, pts: &mut Vec<DVec3>| {
        if !pts.iter().any(|q| (q - p).length_squared() < 1e-18) {
            pts.push(p);
        }
    };
    for e in &edges {
        let ed = brep.edge(e.clone());
        push(brep.vertex(ed.first.clone()).point, &mut pts);
        push(brep.vertex(ed.last.clone()).point, &mut pts);
        if let Some(c) = &ed.curve {
            let [a, b] = ed.range;
            for s in 1..=8 {
                let t = a + (b - a) * (s as f64) / 8.0;
                push(c.point_at(t), &mut pts);
            }
        }
    }
    if pts.len() < 3 {
        return None;
    }
    let p0 = pts[0];
    let mut plane = None;
    'outer: for i in 1..pts.len() {
        for j in (i + 1)..pts.len() {
            let n = (pts[i] - p0).cross(pts[j] - p0);
            if n.length_squared() > 1e-24 {
                plane = Some(rcad_kernel::geom::Plane::new(p0, n));
                break 'outer;
            }
        }
    }
    let plane = plane?;
    let face = brep.add_tface(
        Some(Surface3::Plane(plane)),
        wire.clone(),
        vec![],
        None,
        None,
        vec![],
        false,
    );
    // BRepLib::SameParameter: planar 2D curves on the cap face.
    let fkey = (face.ptr_id(), face.location);
    let to_uv = |p: DVec3| {
        DVec3::new(
            (p - plane.origin).dot(plane.u_dir),
            (p - plane.origin).dot(plane.v_dir),
            0.0,
        )
    };
    for e in &edges {
        let ed = brep.edge(e.clone());
        let uv_a = to_uv(brep.vertex(ed.first.clone()).point);
        let uv_b = to_uv(brep.vertex(ed.last.clone()).point);
        let d = uv_b - uv_a;
        let pc = match &ed.curve {
            // Project the 3D curve control points for curved edges (the
            // closed B-spline profile) so the boundary integral follows the
            // edge; straight edges keep the two-point line.
            Some(Curve3::BSpline(b)) => {
                let cpts: Vec<DVec2> = b.control_points.iter().map(|&p| to_uv(p).truncate()).collect();
                let mut weights = b.weights.clone();
                if weights.is_empty() {
                    weights = vec![1.0; cpts.len()];
                }
                Curve2d::BSpline(rcad_kernel::geom::BSplineCurve2 {
                    degree: b.degree,
                    knots: b.knots.clone(),
                    control_points: cpts,
                    weights,
                })
            }
            _ => Curve2d::Line(Line2d::new(
                uv_a.truncate(),
                if d.length_squared() > 1e-30 {
                    d.truncate().normalize()
                } else {
                    DVec2::X
                },
            )),
        };
        let [ta, tb] = match &ed.curve {
            Some(Curve3::BSpline(_)) => [0.0, 1.0],
            _ => [0.0, (uv_b - uv_a).truncate().length()],
        };
        brep.edge_mut_inplace(e.clone())
            .pcurves
            .insert(fkey, (pc, ta, tb));
    }
    Some(face)
}
