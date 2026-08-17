// OCCT BRepClass3d_BndBoxTree (BRepClass3d_BndBoxTree.hxx / .cxx)
// BVH tree selectors for point and line queries during solid classification.
//
// OCCT builds an NCollection_UBTree of the solid's vertices/edges bounding
// boxes and runs BRepClass3d_BndBoxTreeSelectorPoint/Line over it
// (BRepClass3d_SClassifier::Perform L214-230/L288-372). The selectors' Accept
// predicates (BRepClass3d_BndBoxTree.cxx L39-73/L78-135) perform the actual
// point-vs-edge (Extrema_ExtPC), vertex distance, line-vs-edge (Extrema_ExtCC)
// and line-vs-vertex (Extrema_ExtPElC) checks. rcad does not carry a UBTree
// translation; the tree traversal is replaced by a linear scan over the
// collected vertices/edges with the same Reject/Accept predicates (the UBTree
// only accelerates the scan, the selection result is identical).

use rcad_kernel::base::extrema::closest_point_on_curve;
use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::precision::CONFUSION;
use glam::DVec3;

/// OCCT BRepClass3d_BndBoxTreeSelectorPoint — selects edges/vertices whose
/// tolerance region contains the query point (Accept, BndBoxTree.cxx L39-73).
pub struct BndBoxTreeSelectorPoint {
    point: DVec3,
    found: bool,
    // rcad: linear scan over the shape list (edges first, then vertices),
    // mirroring OCCT's aMapEV order (TopExp::MapShapes(aE, myMapEV) adds the
    // edge then its vertices; BndBoxTree.cxx Accept switches on the shape type).
    edges: Vec<Curve3>,
    edge_tols: Vec<f64>,
    edge_ranges: Vec<[f64; 2]>,
    vertices: Vec<DVec3>,
    vert_tols: Vec<f64>,
}

impl BndBoxTreeSelectorPoint {
    pub fn new(
        edges: Vec<Curve3>,
        edge_tols: Vec<f64>,
        edge_ranges: Vec<[f64; 2]>,
        vertices: Vec<DVec3>,
        vert_tols: Vec<f64>,
    ) -> Self {
        BndBoxTreeSelectorPoint {
            point: DVec3::ZERO,
            found: false,
            edges,
            edge_tols,
            edge_ranges,
            vertices,
            vert_tols,
        }
    }

    /// OCCT: SetCurrentPoint(P) (BndBoxTree.hxx L42).
    pub fn set_current_point(&mut self, p: DVec3) {
        self.point = p;
        self.found = false;
    }

    /// OCCT: the UBTree Select() result — true when Accept returned true.
    pub fn found(&self) -> bool {
        self.found
    }

    pub fn point(&self) -> DVec3 {
        self.point
    }

    /// OCCT: Select over the tree — run the Accept predicate on every element.
    /// myStop in OCCT's Accept exits the traversal early; rcad scans all and
    /// returns as soon as one element accepts (same result).
    pub fn select(&mut self) -> usize {
        self.found = false;
        // OCCT Accept L43-72: edges first.
        let edges = self.edges.clone();
        let edge_tols = self.edge_tols.clone();
        let edge_ranges = self.edge_ranges.clone();
        for (i, curve) in edges.iter().enumerate() {
            if self.accept_edge(curve, edge_tols[i], edge_ranges[i]) {
                self.found = true;
                return 1;
            }
        }
        // OCCT Accept L73-90: vertices.
        let vertices = self.vertices.clone();
        let vert_tols = self.vert_tols.clone();
        for (i, v) in vertices.iter().enumerate() {
            let t = vert_tols[i];
            if (v - self.point).length_squared() < t * t {
                self.found = true;
                return 1;
            }
        }
        0
    }

    /// OCCT BndBoxTree.cxx L43-72 (EDGE branch): Extrema_ExtPC(point, edge)
    /// with the edge tolerance — any extremum closer than the tolerance wins.
    fn accept_edge(&self, curve: &Curve3, tol: f64, range: [f64; 2]) -> bool {
        // Extrema_ExtPC(myP, C, f, l) — the projection restricted to the edge
        // range (closest_point_on_curve_range wraps periodic curves, matching
        // GeomAPI_ProjectPointOnCurve::Init(C, f, l)).
        let proj = rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
            curve,
            self.point,
            range[0],
            range[1],
            64,
        );
        proj.distance < tol
    }
}

/// OCCT BRepClass3d_BndBoxTreeSelectorLine — selects edges/vertices whose
/// tolerance region intersects the query line (Accept, BndBoxTree.cxx L78-135).
pub struct BndBoxTreeSelectorLine {
    line_origin: DVec3,
    line_dir: DVec3,
    max_param: f64,
    is_valid: bool,
    // Edge-Line interferences (BndBoxTree.cxx L96-123): edge idx, param on the
    // edge (P1.Parameter), param on the line (P2.Parameter).
    edge_params: Vec<(usize, f64, f64)>,
    // Vertex-Line interferences (L124-135): vertex idx, param on the line.
    vert_params: Vec<(usize, f64)>,
    edges: Vec<Curve3>,
    edge_tols: Vec<f64>,
    edge_ranges: Vec<[f64; 2]>,
    vertices: Vec<DVec3>,
    vert_tols: Vec<f64>,
}

impl BndBoxTreeSelectorLine {
    pub fn new(
        edges: Vec<Curve3>,
        edge_tols: Vec<f64>,
        edge_ranges: Vec<[f64; 2]>,
        vertices: Vec<DVec3>,
        vert_tols: Vec<f64>,
    ) -> Self {
        BndBoxTreeSelectorLine {
            line_origin: DVec3::ZERO,
            line_dir: DVec3::X,
            max_param: 0.0,
            is_valid: true,
            edge_params: Vec::new(),
            vert_params: Vec::new(),
            edges,
            edge_tols,
            edge_ranges,
            vertices,
            vert_tols,
        }
    }

    /// OCCT: SetCurrentLine(L, MaxParam) (BndBoxTree.hxx L90-94) — loads the
    /// line into myLC with the range [-PConfusion, MaxParam].
    pub fn set_current_line(&mut self, origin: DVec3, dir: DVec3, max_param: f64) {
        self.line_origin = origin;
        self.line_dir = dir;
        self.max_param = max_param;
        self.is_valid = true;
    }

    pub fn clear_results(&mut self) {
        self.edge_params.clear();
        self.vert_params.clear();
        self.is_valid = true;
    }

    pub fn is_valid(&self) -> bool {
        self.is_valid
    }

    pub fn edge_params(&self) -> &[(usize, f64, f64)] {
        &self.edge_params
    }
    pub fn vert_params(&self) -> &[(usize, f64)] {
        &self.vert_params
    }
    pub fn nb_edge_params(&self) -> usize {
        self.edge_params.len()
    }
    pub fn nb_vert_params(&self) -> usize {
        self.vert_params.len()
    }

    /// OCCT: Select over the tree — run the Accept predicate on every element.
    pub fn select(&mut self) -> usize {
        self.edge_params.clear();
        self.vert_params.clear();
        self.is_valid = true;
        let mut n = 0usize;
        let edges = self.edges.clone();
        let edge_tols = self.edge_tols.clone();
        let edge_ranges = self.edge_ranges.clone();
        for (i, curve) in edges.iter().enumerate() {
            if self.accept_edge(i, curve, edge_tols[i], edge_ranges[i]) {
                n += 1;
            }
        }
        let vertices = self.vertices.clone();
        let vert_tols = self.vert_tols.clone();
        for (i, v) in vertices.iter().enumerate() {
            if self.accept_vertex(i, *v, vert_tols[i]) {
                n += 1;
            }
        }
        n
    }

    /// OCCT BndBoxTree.cxx L96-123 (EDGE branch): Extrema_ExtCC(edge, line)
    /// over both ranges. A parallel line makes the selector invalid (the
    /// tangent case cannot be used for classification). Interferences closer
    /// than the edge tolerance are recorded with the parameters on both curves.
    fn accept_edge(&mut self, idx: usize, curve: &Curve3, tol: f64, range: [f64; 2]) -> bool {
        // OCCT Extrema_ExtCC::IsParallel (ExtCC.cxx): true only for the
        // line-line branch (ExtElC sets myIsPar); the line-conic branches
        // leave it false. rcad: line_line_extrema reports parallelism by
        // collapsing to a single constant-distance solution — detect it the
        // same way (ExtElC line-line L268-357).
        if let Curve3::Line(el) = curve {
            let a_d1 = self.line_dir.normalize_or_zero();
            let a_d2 = el.direction.normalize_or_zero();
            if a_d1.length_squared() > 0.5 && a_d2.length_squared() > 0.5 {
                let a_sq_sin_a = 1.0 - a_d1.dot(a_d2) * a_d1.dot(a_d2);
                if a_sq_sin_a < 1e-30 {
                    // OCCT L104-107: IsParallel -> the tangent case is invalid.
                    self.is_valid = false;
                    return false;
                }
            }
        }
        // OCCT Extrema_ExtCC(C, myLC, f, l, myLC.FirstParameter(),
        // myLC.LastParameter()) — the line range is [-PConfusion, MaxParam]
        // (SetCurrentLine L91-93). The first argument is the gp_Lin myLC
        // (Adaptor3d_Curve of the line); ext_cc_line_conic takes the Line3.
        let l0 = -CONFUSION;
        let l1 = self.max_param;
        let ext = rcad_kernel::base::extrema::ext_cc_line_conic(
            &rcad_kernel::geom::Line3 {
                origin: self.line_origin,
                direction: self.line_dir,
            },
            l0,
            l1,
            curve,
            range[0],
            range[1],
        );
        let mut inside = false;
        for (d, p1, p2) in &ext.interior {
            if *d < tol {
                self.edge_params.push((idx, *p1, *p2));
                inside = true;
            }
        }
        inside
    }

    /// OCCT BndBoxTree.cxx L124-135 (VERTEX branch): Extrema_ExtPElC(vertex,
    /// line) — the closest point on the infinite line within the vertex
    /// tolerance; the line parameter is recorded.
    fn accept_vertex(&mut self, idx: usize, v: DVec3, tol: f64) -> bool {
        let dir_sq = self.line_dir.dot(self.line_dir);
        if dir_sq < 1e-30 {
            return false;
        }
        let t = (v - self.line_origin).dot(self.line_dir) / dir_sq;
        let q = self.line_origin + t * self.line_dir;
        if (q - v).length_squared() < tol * tol {
            self.vert_params.push((idx, t));
            return true;
        }
        false
    }
}

/// Helper: bounding box of a vertex (point box, no gap) — used by the BVH
/// Reject predicate (OCCT BRepBndLib::Add on a vertex, Bnd_Box with the
/// vertex tolerance as the gap; the tree traversal Reject is
/// Bnd_Box::IsOut(Pnt/Line) with the accumulated gap).
pub fn vertex_box(p: DVec3, tol: f64) -> BndBox {
    let mut b = BndBox::from_point(p);
    b.set_gap(tol);
    b
}

/// Helper: bounding box of an edge curve (sampled curve box, gap = edge
/// tolerance) — semantic equivalent of OCCT BRepBndLib::Add(edge).
pub fn edge_box(curve: &Curve3, range: [f64; 2], tol: f64) -> BndBox {
    let mut b = BndBox::new();
    let t1 = range[0];
    let t2 = range[1];
    if t1.is_finite() && t2.is_finite() {
        for k in 0..=16 {
            let t = t1 + (t2 - t1) * (k as f64) / 16.0;
            b.add_point(curve.point_at(t));
        }
    } else if t1.is_finite() {
        for k in 0..=16 {
            let t = t1 + k as f64;
            b.add_point(curve.point_at(t));
        }
    } else if t2.is_finite() {
        for k in 0..=16 {
            let t = t2 - k as f64;
            b.add_point(curve.point_at(t));
        }
    }
    b.set_gap(tol);
    b
}
