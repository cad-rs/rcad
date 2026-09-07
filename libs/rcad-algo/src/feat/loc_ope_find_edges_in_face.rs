// OCCT LocOpe_FindEdgesInFace.hxx L89-113 + LocOpe_FindEdgesInFace.cxx
// L516-667 + LocOpe_FindEdgesInFace.lxx L688-723 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_FindEdgesInFace.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_FindEdgesInFace.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_FindEdgesInFace.lxx
//
// OCCT inheritance chain: none (standalone value class).
//
// Architecture differences (referenced from the affected functions):
// 1. BRep_Tool::Curve/Surface with TopLoc_Location — rcad standalone feat
//    shapes carry geometry on the TShape with identity locations; the
//    OCCT location transform branch applies trivially (same reduction as
//    loc_ope_find_edges.rs).
// 2. gp_Pln::Contains / gp_Ax1::IsParallel / gp_Ax1::IsCoaxial /
//    gp_Ax3::IsCoplanar — the gp package has no standalone rcad module yet;
//    the predicates are re-hosted below as pure-math helpers with the
//    standard gp formulas (gp_Pln.cxx / gp_Ax1.cxx / gp_Ax3.cxx).
// 3. ElSLib::Parameters(Pln/Cylinder, P, U, V) — re-hosted as pure-math
//    analytic helpers (ElSLib.cxx); geom_lib::parameters_surface is the
//    numeric solver and is NOT used (OCCT is analytic here).
// 4. BRepTopAdaptor_TopolTool (TPT) — the 2D point-in-face classifier needs
//    the face topology through BRepAdaptor_Surface, which rcad only carries
//    for DS-anchored faces (topalgo::brep_top_adaptor is DS/pool-bound).
//    Until that lands, the local Tpt stand-in reports TopAbs_OUT for every
//    classify (the conservative missing-dependency behaviour, same pattern
//    as the 3a IsDone()==false stubs) — only the shared-edge branch (partage
//    d'edge) of Set keeps collecting. GAP: wire topol_tool_brep in.
// 5. NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> is a HashMap
//    keyed by (TShape ptr, Location) — orientation ignored, the same key
//    scheme as feat::brep_feat_builder::OcctShapeMap.
//
// first consumer: BRepFeat_Form family (3b) — BRepFeat_MakePipe /
// MakeRevol / MakeLinearForm use LocOpe_FindEdgesInFace to select the
// profile edges lying on the support face.

use crate::feat::brep_feat_builder::explorer;
use glam::DVec2;
use glam::DVec3;
use rcad_kernel::geom::{ Curve3, CurveEval, Surface3 };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{ ShapeType, State, TShape };
use std::collections::HashMap;

/// OCCT gp_Pln::Contains(gp_Lin, LinTol, AngTol) (gp_Pln.cxx) — the plane
/// contains the line when the line location lies in the plane (within
/// LinTol) and the line direction is parallel to the plane (within AngTol).
fn gp_pln_contains_line(
    pl_origin: DVec3,
    pl_normal: DVec3,
    li_origin: DVec3,
    li_dir: DVec3,
    lin_tol: f64,
    ang_tol: f64,
) -> bool {
    // OCCT gp_Pln.cxx Contains: distance of the point to the plane.
    let dist = (li_origin - pl_origin).dot(pl_normal);
    if dist.abs() > lin_tol {
        return false;
    }
    // OCCT gp_Pln.cxx: direction parallel to the plane.
    li_dir.dot(pl_normal).abs() <= ang_tol
}

/// OCCT gp_Ax1::IsParallel(Other, AngTol) (gp_Ax1.cxx) — directions parallel
/// (or anti-parallel) within AngTol.
fn gp_ax1_is_parallel(d1: DVec3, d2: DVec3, ang_tol: f64) -> bool {
    d1.cross(d2).length() <= ang_tol
}

/// OCCT gp_Lin::Distance(P) (gp_Lin.cxx) — distance of a point to the line.
fn gp_lin_distance(li_origin: DVec3, li_dir: DVec3, p: DVec3) -> f64 {
    let d = p - li_origin;
    (d - li_dir * d.dot(li_dir)).length()
}

/// OCCT gp_Ax1::IsCoaxial(Other, AngTol, LinTol) (gp_Ax1.cxx) — parallel
/// directions and the other location on the first axis (within LinTol).
fn gp_ax1_is_coaxial(p1: DVec3, d1: DVec3, p2: DVec3, d2: DVec3, ang_tol: f64, lin_tol: f64) -> bool {
    if !gp_ax1_is_parallel(d1, d2, ang_tol) {
        return false;
    }
    let d = p2 - p1;
    let n = d - d1 * d.dot(d1);
    n.length() <= lin_tol
}

/// OCCT gp_Ax3::IsCoplanar(Other, LinTol, AngTol) (gp_Ax3.cxx) — parallel
/// normals and coplanar locations.
fn gp_ax3_is_coplanar(
    n1: DVec3,
    p1: DVec3,
    n2: DVec3,
    p2: DVec3,
    lin_tol: f64,
    ang_tol: f64,
) -> bool {
    if !gp_ax1_is_parallel(n1, n2, ang_tol) {
        return false;
    }
    // OCCT gp_Ax3.cxx: the other location lies in this plane (and vice versa).
    (p2 - p1).dot(n1).abs() <= lin_tol && (p2 - p1).dot(n2).abs() <= lin_tol
}

/// OCCT ElSLib::Parameters(gp_Pln, P, U, V) (ElSLib.cxx) — planar projection
/// onto the plane axes: U = (P-Loc).XDir, V = (P-Loc).YDir.
fn elslib_parameters_plane(pl: &rcad_kernel::geom::Plane, p: DVec3) -> (f64, f64) {
    let d = p - pl.origin;
    (d.dot(pl.u_dir), d.dot(pl.v_dir))
}

/// OCCT ElSLib::Parameters(gp_Cylinder, P, U, V) (ElSLib.cxx) —
/// U = atan2((P-Loc).YDir, (P-Loc).XDir), V = (P-Loc).AxisDir.
fn elslib_parameters_cylinder(cy: &rcad_kernel::geom::CylindricalSurface, p: DVec3) -> (f64, f64) {
    let d = p - cy.origin;
    let x_dir = cy.ref_dir;
    let y_dir = cy.axis.cross(x_dir);
    (d.dot(y_dir).atan2(d.dot(x_dir)), d.dot(cy.axis))
}

/// OCCT BRep_Tool::Surface(F) with the Geom_RectangularTrimmedSurface basis
/// strip (cxx L536-542).
fn brep_tool_surface_basis(face: &Shape) -> Option<Surface3> {
    match face.data.as_ref() {
        TShape::Face(fd) => {
            let s = fd.surface.as_ref()?.clone();
            Some(match s {
                // OCCT cxx L538-542: Geom_RectangularTrimmedSurface -> basis.
                Surface3::Trimmed(ts) => (*ts.basis).clone(),
                _ => s,
            })
        }
        _ => None,
    }
}

/// OCCT BRep_Tool::Curve(edg, Loc, f, l) (architecture difference #1).
fn brep_tool_curve(edg: &Shape) -> Option<(Curve3, f64, f64)> {
    match edg.data.as_ref() {
        TShape::Edge(ed) => {
            let c = ed.curve.as_ref()?;
            Some((c.clone(), ed.range[0], ed.range[1]))
        }
        _ => None,
    }
}

/// OCCT C->Value(u) — curve evaluation through the rcad CurveEval vehicle.
fn curve_value(c: &Curve3, u: f64) -> DVec3 {
    c.point_at(u)
}

/// OCCT BRepTopAdaptor_TopolTool stand-in (architecture difference #4) —
/// every classify reports TopAbs_OUT until the standalone-face classifier
/// lands.
struct Tpt;

impl Tpt {
    fn new() -> Self {
        Tpt
    }
    /// OCCT BRepTopAdaptor_TopolTool::Classify(P, Tol) — missing-dependency
    /// stub returning TopAbs_OUT (see architecture difference #4).
    fn classify(&mut self, _p: DVec2, _tol: f64) -> State {
        State::Out
    }
}

/// OCCT LocOpe_FindEdgesInFace (LocOpe_FindEdgesInFace.hxx L89-113).
pub struct LocOpeFindEdgesInFace {
    my_shape: Shape,      // OCCT: myShape (TopoDS_Shape)
    my_face: Shape,       // OCCT: myFace (TopoDS_Face)
    my_list: Vec<Shape>,  // OCCT: myList (NCollection_List<TopoDS_Shape>)
    my_it: Option<usize>, // OCCT: myIt (list iterator; None = uninitialized)
}

impl Default for LocOpeFindEdgesInFace {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeFindEdgesInFace {
    /// OCCT LocOpe_FindEdgesInFace::LocOpe_FindEdgesInFace() (lxx L688).
    pub fn new() -> Self {
        LocOpeFindEdgesInFace {
            my_shape: Shape::null(),
            my_face: Shape::null(),
            my_list: Vec::new(),
            my_it: None,
        }
    }

    /// OCCT LocOpe_FindEdgesInFace::LocOpe_FindEdgesInFace(S, F)
    /// (lxx L692-695) — Rust has no overloading: the value constructor
    /// carries the `_shape_face` suffix.
    pub fn new_shape_face(the_s: &Shape, the_f: &Shape) -> Self {
        let mut res = LocOpeFindEdgesInFace::new();
        res.set(the_s, the_f);
        res
    }

    /// OCCT LocOpe_FindEdgesInFace::Set(Sh, F) (cxx L516-667).
    pub fn set(&mut self, the_sh: &Shape, the_f: &Shape) {
        self.my_shape = the_sh.clone();
        self.my_face = the_f.clone();
        self.my_list.clear();

        // OCCT cxx L522: NCollection_Map<TopoDS_Shape,...> M.
        let mut m: HashMap<(u64, u32), ()> = HashMap::new();
        // OCCT cxx L533-534.
        let tol = rcad_kernel::precision::CONFUSION;
        let tolang = rcad_kernel::precision::ANGULAR;

        // OCCT cxx L536-542: S = BRep_Tool::Surface(F); trimmed -> basis.
        let Some(s) = brep_tool_surface_basis(the_f) else {
            return;
        };
        // OCCT cxx L544-547: only Plane / CylindricalSurface — "pour le moment".
        if !matches!(s, Surface3::Plane(_) | Surface3::Cylinder(_)) {
            return;
        }

        // OCCT cxx L549-550: BRepAdaptor_Surface + BRepTopAdaptor_TopolTool.
        // Architecture difference #4: the classifier is not wired for
        // standalone faces; Tpt::classify reports TopAbs_OUT (conservative).
        let mut tpt = Tpt::new();

        // OCCT cxx L552:
        for edg in explorer(&self.my_shape, ShapeType::Edge, ShapeType::Shape) {
            // OCCT cxx L554.
            let mut to_add = false;
            // OCCT cxx L556-559: M.Contains(edg) -> continue.
            if m.contains_key(&(edg.ptr_id(), edg.location)) {
                continue;
            }

            // OCCT cxx L561-567: is the edge shared with myFace?
            let mut shared = false;
            for expf in explorer(&self.my_face, ShapeType::Edge, ShapeType::Shape) {
                if expf.ptr_id() == edg.ptr_id() && expf.location == edg.location {
                    shared = true;
                    break;
                }
            }
            if shared {
                // OCCT cxx L568-573: "partage d'edge".
                m.insert((edg.ptr_id(), edg.location), ());
                self.my_list.push(edg.clone());
                continue;
            }

            // OCCT cxx L575-582: C = BRep_Tool::Curve + Transformed; trimmed.
            let Some((mut c, f, l)) = brep_tool_curve(&edg) else {
                continue;
            };
            // OCCT cxx L578-582: trimmed -> basis.
            if matches!(c, Curve3::Trimmed(_)) {
                c = crate::feat::loc_ope_find_edges::basis_curve(&c);
            }
            // OCCT cxx L583-586: only Line / Circle — "pour le moment".
            if !matches!(c, Curve3::Line(_) | Curve3::Circle(_)) {
                continue;
            }
            // OCCT cxx L587-594: pl / cy extraction (cx guarded per branch —
            // OCCT declares pl/cy up front and assigns per surface kind).
            let pl = if let Surface3::Plane(p) = &s { Some(*p) } else { None };
            let cy = if let Surface3::Cylinder(cyl) = &s { Some(*cyl) } else { None };

            // OCCT cxx L596-614: Tc == Geom_Line branch.
            if let Curve3::Line(li) = &c {
                if let Some(pl) = &pl {
                    if gp_pln_contains_line(pl.origin, pl.normal, li.origin, li.direction, tol, tolang) {
                        to_add = true;
                    }
                } else {
                    // OCCT cxx L606-613: Ts == Geom_CylindricalSurface.
                    if let Some(cy) = &cy {
                        if gp_ax1_is_parallel(cy.axis, li.direction, tolang)
                            && ((gp_lin_distance(li.origin, li.direction, cy.origin) - cy.radius).abs() < tol)
                        {
                            to_add = true;
                        }
                    }
                }
            } else if let Curve3::Circle(ci) = &c {
                // OCCT cxx L615-633: Tc == Geom_Circle branch.
                if let Some(pl) = &pl {
                    if gp_ax3_is_coplanar(pl.normal, pl.origin, ci.normal, ci.center, tol, tolang) {
                        to_add = true;
                    }
                } else if let Some(cy) = &cy {
                    if (cy.radius - ci.radius).abs() < tol
                        && gp_ax1_is_coaxial(cy.origin, cy.axis, ci.center, ci.normal, tolang, tol)
                    {
                        to_add = true;
                    }
                }
            }

            if to_add {
                // OCCT cxx L636-659: "On classifie 3 points."
                let mut p = [DVec3::ZERO; 3];
                p[0] = curve_value(&c, f);
                p[1] = curve_value(&c, l);
                p[2] = curve_value(&c, (f + l) / 2.0);
                let mut i = 0;
                while i < 3 {
                    // OCCT cxx L647-654: ElSLib::Parameters(pl/cy, p[i], U, V).
                    let (u, v) = if let Some(pl) = &pl {
                        elslib_parameters_plane(pl, p[i])
                    } else {
                        elslib_parameters_cylinder(cy.as_ref().expect("cylinder surface"), p[i])
                    };
                    // OCCT cxx L655: TPT.Classify(gp_Pnt2d(U,V),
                    // Precision::Confusion()) == TopAbs_OUT -> break.
                    if tpt.classify(DVec2::new(u, v), rcad_kernel::precision::CONFUSION) == State::Out
                    {
                        break;
                    }
                    i += 1;
                }
                // OCCT cxx L660-664: if (i >= 3).
                if i >= 3 {
                    m.insert((edg.ptr_id(), edg.location), ());
                    self.my_list.push(edg.clone());
                }
            }
        }
    }

    /// OCCT LocOpe_FindEdgesInFace::Init() (lxx L699-702).
    pub fn init(&mut self) {
        self.my_it = if self.my_list.is_empty() { None } else { Some(0) };
    }

    /// OCCT LocOpe_FindEdgesInFace::More() (lxx L706-709).
    pub fn more(&self) -> bool {
        self.my_it.map_or(false, |i| i < self.my_list.len())
    }

    /// OCCT LocOpe_FindEdgesInFace::Edge() (lxx L713-716) — TopoDS::Edge
    /// cast; the list only ever receives edges (cxx L571/L663).
    pub fn edge(&self) -> Shape {
        self.my_list[self.my_it.expect("More() == true")].clone()
    }

    /// OCCT LocOpe_FindEdgesInFace::Next() (lxx L720-723).
    pub fn next(&mut self) {
        self.my_it = self.my_it.map(|i| i + 1);
    }
}
